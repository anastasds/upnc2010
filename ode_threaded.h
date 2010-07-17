#ifndef ODE_THREADED_H
#define ODE_THREADED_H

#define NUM_THREADS 2

void * ode_update_neurons_thread(void * thread_params);
int hh_ode_threaded(double t, const double y[], double f[], void *params);
int ode_run_threaded(struct network * network, double t, double t1, double step_size, double error);

#endif
