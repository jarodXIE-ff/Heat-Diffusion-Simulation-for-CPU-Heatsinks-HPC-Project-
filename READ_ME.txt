I'LL BRIEFLY EXPLAIN HOW TO USE THE CODE

First of all, you need to have access to the Grid5000 servers, and with your login and password, connect with an SSH connection.
Then send the folder with the command scp -r PROJET_PPAR USERNAME@access.grid5000.fr:nancy/ (you need to be in the right directory in your terminal).

Then type make in the terminal to compile the various .c files 

Allocate nodes to run the parallelized code with oarsub -l core=100 -I for example.

Run the code with mpiexec -hostfile $OAR_NODE_FILE -n 75 HeatSink1d > test1d.txt (75 here taken abitrairement)

To retrieve the average times of these simulations: python3 Average_Times.py Normal1d.txt (if you've set the NORMAL mode, for example)

For visual rendering: python3 rendu_picture_steady.py test1d.txt n m o (n,m,o are given)

To display a curve showing the current time mode in relation to the number of nodes: python3 courbe.py Normal1d.txt (if you've set the NORMAL mode, for example) 


THX AND GOOD LUCK :)

(I've left the image rendering with Normal mode and the curve for Normal mode as the calculation takes a long time)