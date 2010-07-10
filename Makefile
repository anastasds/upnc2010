CC = gcc
CFLAGS = -O3
GSLFLAGS = -lgsl -lgslcblas -lm

all : proj1

proj1: neuron.o general.o ode.o main.o
	$(CC) $(CFLAGS) $(GSLFLAGS) neuron.o general.o ode.o main.o -o proj1

neuron.o: neuron.c
	$(CC) $(CFLAGS) -c neuron.c

general.o: general.c
	$(CC) $(CFLAGS) -c general.c

ode.o: ode.c
	$(CC) $(CFLAGS) $(GSLFLAGS) -c ode.c

main.o: main.c
	$(CC) $(CFLAGS) $(GSLFLAGS) -c main.c
clean:
	rm *.o proj1