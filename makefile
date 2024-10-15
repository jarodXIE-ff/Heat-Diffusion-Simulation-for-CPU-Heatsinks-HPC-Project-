CC = mpicc
CFLAGS = -Wall -O3
LIBS = -lm

all: HeatSink1d SequentialCode

SequentialCode: SequentialCode.c
	$(CC) $(CFLAGS) SequentialCode.c -o SequentialCode $(LIBS)

HeatSink1d: HeatSink1d.c
	$(CC) $(CFLAGS) HeatSink1d.c -o HeatSink1d $(LIBS)
clean:
	rm -f HeatSink1d