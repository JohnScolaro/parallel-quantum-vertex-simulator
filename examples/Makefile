PROJ = helloflops2
SRC = $(PROJ).c
OBJ = $(PROJ).o
CC = icc
FLAGS = -qopenmp -static
LINKS = -lpthread -lm -ldl
OPT = -qopt-report=3 -qopt-report-phase=vec -O3

all:
	$(CC) $(FLAGS) $(OPT) -mmic $(SRC) -o $(PROJ) $(LINKS)

test:
	icc -mkl -O3 -mmic -openmp -L opt/intel/lib/mic -Wno-unknown-pragmas -std=c99 -vec-report3 helloflops2.c -o helloflops2test -liomp5

clean:
	rm -f $(PROJ)
	rm -f *.o
	rm -f *.optrpt
