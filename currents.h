#ifndef CURRENTS_H
#define CURRENTS_H

double Na_current(struct network * network, long num_neuron, long num_compartment, double * f, const double * y, double t);
double Kdr_current(struct network * network, long num_neuron, long num_compartment, double * f, const double * y, double t);
double A_current(struct network * network, long num_neuron, long num_compartment, double * f, const double * y, double t);
double KCa_current(struct network * network, long num_neuron, long num_compartment, double * f, const double * y, double t);
double CaT_current(struct network * network, long num_neuron, long num_compartment, double * f, const double * y, double t);
double L_current(struct network * network, long num_neuron, long num_compartment, double * f, const double * y, double t);
double NMDA_current(struct network * network, long num_neuron, long num_compartment, double * f, const double * y, double t, double i_CaT);
double AMPA_current(struct network * network, long num_neuron, long num_compartment, double * f, const double * y, double t);
double presynaptic_activity(struct network * network, long num_neuron, long num_compartment);
double eff(double z);
double eff2(double z);
double eff3(double z);
double heav(double x);
#endif
