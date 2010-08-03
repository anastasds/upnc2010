#ifndef STIMULATE_H
#define STIMULATE_H

struct stimuli * init_stimuli_struct();
void init_stimuli(struct network * network, char * filename);
void stimulate(struct stimuli * stimuli, long num_neuron, long num_compartment, double start, double end, double current);
void destroy_stimuli(struct stimuli * stimuli);
double apply_stimulus(struct network * network, long num_neuron, long num_compartment, double t);
struct stimulus * find_stimulus(struct stimuli * stimuli, long num_neuron, long num_compartment);

#endif
