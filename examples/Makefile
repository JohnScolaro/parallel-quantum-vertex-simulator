PROJ = helloflops2
SRC = $(PROJ).c
OBJ = $(PROJ).o
CC = icc
FLAGS = -qopenmp -static
LINKS = -lpthread -lm -ldl
OPT = -qopt-report=3 -qopt-report-phase=vec -O3
INTEL_LINKS = -Wl,--start-group ${MKLROOT}/lib/mic/libmkl_intel_ilp64.a ${MKLROOT}/lib/mic/libmkl_core.a ${MKLROOT}/lib/mic/libmkl_intel_thread.a -Wl,--end-group -lpthread -lm -ldl
COMPILER_OPTIONS = -DMKL_ILP64 -qopenmp -I${MKLROOT}/include -mmic

all:
	$(CC) $(FLAGS) $(OPT) -mmic $(SRC) -o $(PROJ) $(LINKS)

test:
	icc -mkl -O3 -mmic -qopenmp -L opt/intel/lib/mic -Wno-unknown-pragmas -std=c99 -qopt-report=3 -qopt-report-phase=vec helloflops2.c -o helloflops2test -liomp5 -static

test2:
	icc ${COMPILER_OPTIONS} helloflops2.c -o helloflopstest2 ${INTEL_LINKS} -static


clean:
	rm -f $(PROJ)
	rm -f *.o
	rm -f *.optrpt
