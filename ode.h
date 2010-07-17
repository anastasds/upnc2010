#ifndef ODE_H
#define ODE_H

int hh_ode (double t, const double y[], double f[], void *params);
int ode_run(struct network * network, double t, double t1, double step_size, double error);

#endif
