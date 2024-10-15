#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <mpi.h>
#include <time.h>

/* AUTHOR : Charles Bouillaguet <charles.bouillaguet@lip6.fr>
   USAGE  : compile with -lm (and why not -O3)
            redirect the standard output to a text file
            gcc heatsink.c -O3 -lm -o heatsink
            ./heatsink > steady_state.txt
            then run the indicated python script for graphical rendering

   DISCLAIMER : this code does not claim to an absolute realism.
                this code could be obviously improved, but has been written so as
                to make as clear as possible the physics principle of the simulation.
*/

/* one can change the matter of the heatsink, its size, the power of the CPU, etc. */
#define ALUMINIUM
#define NORMAL /* MEDIUM is faster, and FAST is even faster (for debugging) */
#define DUMP_STEADY_STATE
const double L = 0.15;            /* length (x) of the heatsink (m) */
const double l = 0.12;            /* height (y) of the heatsink (m) */
const double E = 0.008;           /* width (z) of the heatsink (m) */
const double watercooling_T = 20; /* temperature of the fluid for water-cooling, (°C) */
const double CPU_TDP = 280;       /* power dissipated by the CPU (W) */

/* dl: "spatial step" for simulation (m) */
/* dt: "time step" for simulation (s) */
#ifdef FAST
double dl = 0.004;
double dt = 0.004;
#endif

#ifdef MEDIUM
double dl = 0.002;
double dt = 0.002;
#endif

#ifdef NORMAL
double dl = 0.001;
double dt = 0.001;
#endif

#ifdef CHALLENGE
double dl = 0.0001;
double dt = 0.00001;
#endif

/* sink_heat_capacity: specific heat capacity of the heatsink (J / kg / K) */
/* sink_density: density of the heatsink (kg / m^3) */
/* sink_conductivity: thermal conductivity of the heatsink (W / m / K) */
/* euros_per_kg: price of the matter by kilogram */
#ifdef ALUMINIUM
double sink_heat_capacity = 897;
double sink_density = 2710;
double sink_conductivity = 237;
double euros_per_kg = 1.594;
#endif

#ifdef COPPER
double sink_heat_capacity = 385;
double sink_density = 8960;
double sink_conductivity = 390;
double euros_per_kg = 5.469;
#endif

#ifdef GOLD
double sink_heat_capacity = 128;
double sink_density = 19300;
double sink_conductivity = 317;
double euros_per_kg = 47000;
#endif

#ifdef IRON
double sink_heat_capacity = 444;
double sink_density = 7860;
double sink_conductivity = 80;
double euros_per_kg = 0.083;
#endif

const double Stefan_Boltzmann = 5.6703e-8;   /* (W / m^2 / K^4), radiation of black body */
const double heat_transfer_coefficient = 10; /* coefficient of thermal convection (W / m^2 / K) */
double CPU_surface;

/*
 * Return True if the CPU is in contact with the heatsink at the point (x,y).
 * This describes an AMD EPYC "Rome".
 */
static inline bool CPU_shape(double x, double y)
{
    x -= (L - 0.0754) / 2;
    y -= (l - 0.0585) / 2;
    bool small_y_ok = (y > 0.015 && y < 0.025) || (y > 0.0337 && y < 0.0437);
    bool small_x_ok = (x > 0.0113 && x < 0.0186) || (x > 0.0193 && x < 0.0266) || (x > 0.0485 && x < 0.0558) || (x > 0.0566 && x < 0.0639);
    bool big_ok = (x > 0.03 && x < 0.045 && y > 0.0155 && y < 0.0435);
    return big_ok || (small_x_ok && small_y_ok);
}

/* returns the total area of the surface of contact between CPU and heatsink (in m^2) */
double CPU_contact_surface()
{
    double S = 0;
    for (double x = dl / 2; x < L; x += dl)
        for (double y = dl / 2; y < l; y += dl)
            if (CPU_shape(x, y))
                S += dl * dl;
    return S;
}

/* Returns the new temperature of the cell (i, j, k). For this, there is an access to neighbor
 * cells (left, right, top, bottom, front, back), except if (i, j, k) is on the external surface. */
static inline double update_temperature(const double *T, int u, int n, int m, int o, int i, int j, int k)
{
    /* quantity of thermal energy that must be brought to a cell to make it heat up by 1°C */
    const double cell_heat_capacity = sink_heat_capacity * sink_density * dl * dl * dl; /* J.K */
    const double dl2 = dl * dl;
    double thermal_flux = 0;

    if (i > 0)
        thermal_flux += (T[u - m * o] - T[u]) * sink_conductivity * dl; /* neighbor x-1 */
    else
    {
        thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[u], 4);
        thermal_flux -= heat_transfer_coefficient * dl2 * (T[u] - watercooling_T);
    }

    if (i < n - 1)
        thermal_flux += (T[u + m * o] - T[u]) * sink_conductivity * dl; /* neighbor x+1 */
    else
    {
        thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[u], 4);
        thermal_flux -= heat_transfer_coefficient * dl2 * (T[u] - watercooling_T);
    }

    if (j > 0)
        thermal_flux += (T[u - o] - T[u]) * sink_conductivity * dl; /* neighbor y-1 */
    else
    {
        /* Bottom cell: does it receive it from the CPU ? */
        if (CPU_shape(i * dl, k * dl))
            thermal_flux += CPU_TDP / CPU_surface * dl2;
        else
        {
            thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[u], 4);
            thermal_flux -= heat_transfer_coefficient * dl2 * (T[u] - watercooling_T);
        }
    }

    if (j < m - 1)
        thermal_flux += (T[u + o] - T[u]) * sink_conductivity * dl; /* neighbor y+1 */
    else
    {
        thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[u], 4);
        thermal_flux -= heat_transfer_coefficient * dl2 * (T[u] - watercooling_T);
    }

    if (k > 0)
        thermal_flux += (T[u - 1] - T[u]) * sink_conductivity * dl; /* neighbor z-1 */
    else
    {
        thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[u], 4);
        thermal_flux -= heat_transfer_coefficient * dl2 * (T[u] - watercooling_T);
    }

    if (k < o - 1)
        thermal_flux += (T[u + 1] - T[u]) * sink_conductivity * dl; /* neighbor z+1 */
    else
    {
        thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[u], 4);
        thermal_flux -= heat_transfer_coefficient * dl2 * (T[u] - watercooling_T);
    }

    /* adjust temperature depending on the heat flux */
    return T[u] + thermal_flux * dt / cell_heat_capacity;
}

int main()
{
    // Openning the right files for each difficulty
    FILE *fichier = NULL;
    char *fileNames[4] = {"Fast1d.txt", "Medium1d.txt", "Normal1d.txt", "Challenge1d.txt"};
    double indices[4] = {0.004, 0.002, 0.001, 0.0001};

    int i = 0;
    for (i = 0; i <= 3; i++)
    {

        if (dl == indices[i])
        {
            fichier = fopen(fileNames[i], "a");
            break;
        }
    }

    if (fichier == NULL)
    {
        fprintf(stderr, "Impossible d'ouvrir le fichier.\n");
        return 1; // Code d'erreur
    }
    // Initialisation of the time function
    clock_t starting, ending;
    double cpu_time_used;
    starting = clock();

    // Initialisation of the temperature function
    CPU_surface = CPU_contact_surface();
    double V = L * l * E;
    int n = ceil(L / dl);
    int m = ceil(E / dl);
    int o = ceil(l / dl);

    // Everything related to MPI now
    MPI_Request trash;
    int rang, p;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rang);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    int average_size = ceil((double)n / p);
    int start = rang * average_size;
    int end;
    if (average_size * (rang + 1) < n)
    {
        end = average_size * (rang + 1);
    }
    else
    {
        end = n;
    }

    int size = end - start;
    int index = start * m * o;
    int np = ceil((double)n / average_size);

    MPI_Group world_group;
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);
    // MPI Process Group Handling
    MPI_Group new_group;
    MPI_Comm newworld;

    // Check if the total number of processes (p) is greater than the desired number of processes (np)
    if (p > np)
    {
        // Define a range to exclude unnecessary ranks from the world group
        int ranges[3] = {np, p - 1, 1};
        // Create a new MPI group
        MPI_Group_range_excl(world_group, 1, &ranges, &new_group);
        // Create a new MPI communicator based on the new group
        MPI_Comm_create(MPI_COMM_WORLD, new_group, &newworld);
    }
    else
    {
        // If the total number of processes is not greater than np, use the entire world group
        MPI_Comm_create(MPI_COMM_WORLD, world_group, &newworld);
    }
    if (newworld == MPI_COMM_NULL)
    {
        // Terminate MPI and exit the program
        MPI_Finalize();
        exit(0);
    }
    // Adjust the value of p to np
    p = np;

    fprintf(stderr, "HEATSINK\n");
    fprintf(stderr, "\tDimension (cm) [x,y,z] = %.1f x %.1f x %.1f\n", 100 * L, 100 * E, 100 * l);
    fprintf(stderr, "\tVolume = %.1f cm^3\n", V * 1e6);
    fprintf(stderr, "\tWeight = %.2f kg\n", V * sink_density);
    fprintf(stderr, "\tPrice = %.2f €\n", V * sink_density * euros_per_kg);
    fprintf(stderr, "SIMULATION\n");
    fprintf(stderr, "\tGrid (x,y,z) = %d x %d x %d (%.1fMo)\n", n, m, o, 7.6293e-06 * n * m * o);
    fprintf(stderr, "\tdt = %gs\n", dt);
    fprintf(stderr, "CPU\n");
    fprintf(stderr, "\tPower = %.0fW\n", CPU_TDP);
    fprintf(stderr, "\tArea = %.1f cm^2\n", CPU_surface * 10000);

    /* temperature of each cell, in degree Kelvin. */
    double *T;
    double *R;
    if (rang == 0)
    {
        T = (double *)malloc((average_size * p + 1) * m * o * sizeof(*T));
        R = (double *)malloc((average_size * p + 1) * m * o * sizeof(*R));
    }
    else
    {
        T = (double *)malloc((average_size + 2) * m * o * sizeof(*T));
        R = (double *)malloc((average_size + 2) * m * o * sizeof(*R));
    }

    T += m * o;
    R += m * o;
    if (T == NULL || R == NULL)
    {
        perror("Allocation problem with T or R");
        exit(1);
    }

    /* initially the heatsink is at the temperature of the water-cooling fluid */
    for (int u = -m * o; u < (size + 1) * m * o; u++)
        R[u] = T[u] = watercooling_T + 273.15;

    /* let's go! we switch the CPU on and launch the simulation until it reaches a stationary state. */
    double t = 0;
    int n_steps = 0;
    int convergence = 0;

    while (convergence == 0)
    {
        // if t is true, MPI communication is performed
        if (t)
        {
            // Send data to the right neighbor
            if (end < n)
                MPI_Isend(T + (size - 1) * m * o, m * o, MPI_DOUBLE, rang + 1, t, newworld, &trash);

            // Send data to the left neighbor
            if (start)
                MPI_Isend(T, m * o, MPI_DOUBLE, rang - 1, t, newworld, &trash);

            // Receive data from the right neighbor
            if (end < n)
                MPI_Recv(T + (size)*m * o, m * o, MPI_DOUBLE, rang + 1, t, newworld, MPI_STATUS_IGNORE);

            // Receive data from the left neighbor
            if (start)
                MPI_Recv(T - m * o, m * o, MPI_DOUBLE, rang - 1, t, newworld, MPI_STATUS_IGNORE);
        }
        // Update temp
        for (int i = start; i < end; i++)
        {
            int w = i * m * o - index;
            for (int j = 0; j < m; j++)
            {
                int v = w + j * o;
                for (int k = 1; k < o; k++)
                {
                    int u = v + k;
                    // Call the update_temperature function for each cell and the result is stored in R
                    R[u] = update_temperature(T, u, n, m, o, i, j, k);
                }
            }
        }
        /* each second, we test the convergence, and print a short progress report */
        if (n_steps % ((int)(1 / dt)) == 0)
        {
            double delta_T = 0;
            for (int u = 0; u < size * m * o; u++)
            {
                delta_T += (R[u] - T[u]) * (R[u] - T[u]);
            }
            // Sum and square root the squared differences globally using MPI
            MPI_Allreduce(&delta_T, &delta_T, 1, MPI_DOUBLE, MPI_SUM, newworld);
            delta_T = sqrt(delta_T) / dt;

            // Print progress report
            if (rang == 0)
                fprintf(stderr, "t = %.3fs ; convergence = %g\n", t, delta_T);

            // Check for convergence criterion
            if (delta_T < 0.1)
                convergence = 1;
        }

        // New temp in R
        double *temp = R;
        R = T;
        T = temp;

        t += dt;
        n_steps += 1;
    }
    MPI_Gather(T, m * o * average_size, MPI_DOUBLE,
               T, m * o * average_size, MPI_DOUBLE, 0, newworld);
#ifdef DUMP_STEADY_STATE
    if (!rang)
    {
        printf("###### STEADY STATE; t = %.1f\n", t);
        for (int k = 0; k < o; k++)
        { // z
            printf("# z = %g\n", k * dl);
            for (int j = 0; j < m; j++)
            { // y
                for (int i = 0; i < n; i++)
                { // x
                    printf("%.1f ", T[i * o * m + j * o + k] - 273.15);
                }
                printf("\n");
            }
        }
        printf("\n");
        ending = clock();
        cpu_time_used = ((double)(ending - starting)) / CLOCKS_PER_SEC;
        printf("Le temps d'execution de la fonction est %f secondes avec %d CPU\n", cpu_time_used, p);

        // Recording the time for n CPU
        fprintf(fichier, "%f %d", cpu_time_used, p);
        fputc('\n', fichier);
        fclose(fichier);
        fprintf(stderr, "For graphical rendering: python3 rendu_picture_steady.py [filename.txt] %d %d %d\n", n, m, o);
#endif
    }
    MPI_Finalize();
    exit(EXIT_SUCCESS);
}
