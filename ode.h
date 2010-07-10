#ifndef DEFINITIONS
#define DEFINITIONS
#include "definitions.h"
#endif

#ifndef INCLUDES
#define INCLUDES
#include "includes.h"
#endif

#ifndef ODE_H
#define ODE_H

int ode_test(struct neuron_state * neuron_state, struct neuron_params * neuron_params);
int hh_ode_jac (double t, const double y[], double *dfdy, double dfdt[], void *params);
int hh_ode (double t, const double y[], double f[], void *params);

#endif
