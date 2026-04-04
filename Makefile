CC = gcc
FC = gfortran
CFLAGS = -O3 -fopenmp
FFLAGS = -O3 -fopenmp

all: ring_detector

ring_detector: src/main.c src/ring_engine.f90
	$(FC) $(FFLAGS) -c src/ring_engine.f90 -o ring_engine.o
	$(CC) $(CFLAGS) -c src/main.c -o main.o
	$(FC) $(FFLAGS) main.o ring_engine.o -o ring_detector

lib: src/ring_engine.f90
	$(FC) -shared -fPIC $(FFLAGS) src/ring_engine.f90 -o ringdetect/libringengine.so

clean:
	rm -f *.o *.mod ring_detector ringdetect/libringengine.so
