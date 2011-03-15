#ifndef CURRENTS_H
#define CURRENTS_H

double Na_current(struct network * network, long num_neuron, long num_compartment, double * f, const double * y, double t);
double Kdr_current(struct network * network, long num_neuron, long num_compartment, double * f, const double * y, double t);
double A_current(struct network * network, long num_neuron, long num_compartment, double * f, const double * y, double t);
double KCa_current(struct network * network, long num_neuron, long num_compartment, double * f, const double * y, double t);
double CaT_current(struct network * network, long num_neuron, long num_compartment, double * f, const double * y, double t);
double CaT_f(double z);
double CaL_current(struct network * network, long num_neuron, long num_compartment, double * f, const double * y, double t);
double L_current(struct network * network, long num_neuron, long num_compartment, double * f, const double * y, double t);
double NMDA_current(struct network * network, long num_neuron, long num_compartment, double * f, const double * y, double t, double i_CaL);
double AMPA_current(struct network * network, long num_neuron, long num_compartment, double * f, const double * y, double t);
double presynaptic_activity(struct network * network, long num_neuron, long num_compartment, long num_link,  double * f, const double * y, double t);
double square_wave_f_pre(struct network * network, long num_neuron, long num_compartment, double * f, const double * y, double t);
double dayan_abbott_f_pre(struct network * network, long num_neuron, long num_compartment, double * f, const double * y, double t);
double destexhe_transmitter_release(struct network * network, long num_neuron, long num_compartment, long num_link, double * f, const double * y, double t);
long calc_y_array_offset(struct network * network, long num_neuron, long num_compartment);
double eff(double z);
double eff2(double z);
double eff3(double z);
double heav(double x);
double max(double a, double b);
double min(double a, double b);
#endif
