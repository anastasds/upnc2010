#ifndef ODE_H
#define ODE_H

#ifdef THREADED
  #define NUM_THREADS 2
int hh_ode_threaded(double t, const double y[], double f[], void *params);
#endif

void * ode_update_neurons_threaded(void * thread_params);
void ode_update_neurons(struct network * network, long start, long num, const double * y, double * f, double t);
int ode_run(struct network * network, double t, double t1, double step_size, double error);
int hh_ode (double t, const double y[], double f[], void *params);
long get_ode_system_dimension(struct network * network);
void output_data(struct network * network, double t, const double * y);

#endif
