#ifndef CURRENTS_H
#define CURRENTS_H

double Na_current(struct network * network, long num_neuron, long num_compartment, double * f, const double * y);
double Kdr_current(struct network * network, long num_neuron, long num_compartment, double * f, const double * y);
double A_current(struct network * network, long num_neuron, long num_compartment, double * f, const double * y);
double KCa_current(struct network * network, long num_neuron, long num_compartment, double * f, const double * y);
double CaT_current(struct network * network, long num_neuron, long num_compartment, double * f, const double * y);
double L_current(struct network * network, long num_neuron, long num_compartment, double * f, const double * y);
double NMDA_current(struct network * network, long num_neuron, long num_compartment, double * f, const double * y);
double AMPA_current(struct network * network, long num_neuron, long num_compartment, double * f, const double * y);

#endif
