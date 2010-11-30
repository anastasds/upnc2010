#ifndef STIMULATE_H
#define STIMULATE_H

struct stimuli * init_stimuli_struct();
void init_stimuli(struct network * network, char * filename);
//void stimulate(struct stimuli * stimuli, double start, double end, int direct, double current);
void append_stimulus(struct stimuli * stimuli, struct stimulus * stimulus);
void stimulate(struct network * network, long num_neuron, long num_compartment, double start, double end, int direct, double current);
void destroy_stimuli(struct stimuli * stimuli);
struct stimulus * apply_stimulus(struct network * network, long num_neuron, long num_compartment, double t);


#endif
