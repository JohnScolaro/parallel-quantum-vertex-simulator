PROJ = helloflops2
SRC = $(PROJ).c
OBJ = $(PROJ).o
CC = icc
FLAGS = -qopenmp -static
LINKS = -lpthread -lm -ldl
OPT = -qopt-report=3 -qopt-report-phase=vec -O3

all:
	$(CC) $(FLAGS) $(OPT) -mmic $(SRC) -o $(PROJ) $(LINKS)

clean:
	rm -f $(PROJ)
	rm -f *.o
	rm -f *.optrpt
