CC = gcc
CFLAGS = -O3
GSLFLAGS = -lgsl -lgslcblas -lm
COMMON_INCLUDES = definitions.h includes.h
MATH_INCLUDES = math_includes.h
COMMON_OBJECTS = neuron.o general.o
UNTHREADED_OBJECTS = $(COMMON_OBJECTS) main_unthreaded.o ode_unthreaded.o
THREADED_OBJECTS = $(COMMON_OBJECTS) ode_threaded.o main_threaded.o
THREADING = -DTHREADED -pthread

all : threaded unthreaded

threaded: $(COMMON_OBJECTS) $(THREADED_OBJECTS) Makefile
	$(CC) $(CFLAGS) $(GSLFLAGS) $(THREADED_OBJECTS) $(THREADING) -o threaded

unthreaded: $(COMMON_OBJECTS) $(UNTHREADED_OBJECTS) Makefile
	$(CC) $(CFLAGS) $(GSLFLAGS) $(UNTHREADED_OBJECTS) -o unthreaded

neuron.o: $(COMMON_INCLUDES) neuron.h neuron.c
	$(CC) $(CFLAGS) -c neuron.c

general.o: $(COMMON_INCLUDES) general.h general.c
	$(CC) $(CFLAGS) -c general.c

ode_unthreaded.o: $(COMMON_INCLUDES) $(MATH_INCLUDES) ode_unthreaded.h ode_unthreaded.c
	$(CC) $(CFLAGS) $(GSLFLAGS) -c ode_unthreaded.c -o ode_unthreaded.o

ode_threaded.o: $(COMMON_INCLUDES) $(MATH_INCLUDES) ode_threaded.h ode_threaded.c
	$(CC) $(CFLAGS) $(GSLFLAGS) $(THREADING) -c ode_threaded.c

main_unthreaded.o: $(COMMON_INCLUDES) general.h neuron.h ode_unthreaded.h main.c
	$(CC) $(CFLAGS) $(GSLFLAGS) -c main.c -o main_unthreaded.o

main_threaded.o: $(COMMON_INCLUDES) general.h neuron.h ode_threaded.h main.c
	$(CC) $(CFLAGS) $(GSLFLAGS) $(THREADING) -DTHREADED -c main.c -o main_threaded.o

clean:
	rm *.o threaded unthreaded