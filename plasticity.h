#ifndef PLASTICITY_H
#define PLASTICITY_H

double evolve_P(struct network * network, const double * y, long offset);
double evolve_V(struct network * network, const double * y, long offset);
double evolve_A(struct network * network, const double * y, long offset);
double evolve_B(struct network * network, const double * y, long offset);
double evolve_D(struct network * network, const double * y, long offset);
double evolve_W(struct network * network, const double * y, long offset);
double p_sigma(struct network * network, double x);
double a_sigma(struct network * network, double x);
double v_sigma(struct network * network, double x);
double d_sigma(struct network * network, double x);
double b_sigma(struct network * network, double x);

#endif
