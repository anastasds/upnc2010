CC = gcc
CFLAGS = -Wall -O0 -g
GSLFLAGS = -lgsl -lgslcblas -lm
COMMON_INCLUDES = definitions.h includes.h
MATH_INCLUDES = math_includes.h
COMMON_OBJECTS = neuron.o general.o main.o currents.o
UNTHREADED_OBJECTS = $(COMMON_OBJECTS) ode_unthreaded.o
THREADED_OBJECTS = $(COMMON_OBJECTS) ode_threaded.o

all : threaded unthreaded

threaded: $(COMMON_OBJECTS) $(THREADED_OBJECTS) Makefile
	$(CC) $(CFLAGS) $(GSLFLAGS) $(THREADED_OBJECTS) -pthread -o threaded

unthreaded: $(COMMON_OBJECTS) $(UNTHREADED_OBJECTS) Makefile
	$(CC) $(CFLAGS) $(GSLFLAGS) $(UNTHREADED_OBJECTS) -o unthreaded

neuron.o: $(COMMON_INCLUDES) neuron.h neuron.c Makefile
	$(CC) $(CFLAGS) -c neuron.c

general.o: $(COMMON_INCLUDES) general.h general.c Makefile
	$(CC) $(CFLAGS) -c general.c

ode_unthreaded.o: $(COMMON_INCLUDES) $(MATH_INCLUDES) currents.c ode.h ode.c Makefile
	$(CC) $(CFLAGS) $(GSLFLAGS) -c ode.c -o ode_unthreaded.o

ode_threaded.o: $(COMMON_INCLUDES) $(MATH_INCLUDES) currents.c ode.h ode.c Makefile
	$(CC) $(CFLAGS) $(GSLFLAGS) -pthread -DTHREADED -c ode.c -o ode_threaded.o

currents.o: $(COMMON_INCLUDES) $(MATH_INCLUDES) neuron.h currents.h currents.c Makefile
	$(CC) $(CFLAGS) $(GSLFLAGS) -c currents.c -o currents.o

main.o: $(COMMON_INCLUDES) general.h neuron.h ode.h main.c Makefile
	$(CC) $(CFLAGS) $(GSLFLAGS) $(THREADING) -c main.c

clean:
	-rm *.o threaded unthreaded