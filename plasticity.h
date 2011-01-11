#ifndef PLASTICITY_H
#define PLASTICITY_H

double evolve_P(struct network * network, long i, long j, long k, const double * y);
double evolve_V(struct network * network, long i, long j, long k, const double * y);
double evolve_A(struct network * network, long i, long j, long k, const double * y);
double evolve_B(struct network * network, long i, long j, long k, const double * y);
double evolve_D(struct network * network, long i, long j, long k, const double * y);
double evolve_W(struct network * network, long i, long j, long k, const double * y);
double p_sigma(struct network * network, double x, long i, long j, long k);
double a_sigma(struct network * network, double x, long i, long j, long k);
double v_sigma(struct network * network, double x, long i, long j, long k);
double d_sigma(struct network * network, double x, long i, long j, long k);
double b_sigma(struct network * network, double x, long i, long j, long k);

#endif
