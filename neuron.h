#ifndef NEURON_H
#define NEURON_H

struct neuron {
  int num_links;
  struct neuron_state * state;
  struct neuron_params * params;
  struct neuron_link ** links;
};

struct neuron_state {
  int num_params;
  char ** names;
  float * values;
};

struct neuron_params {
  int num_params;
  char ** names;
  float * values;
};

struct neuron_ode_params {
  struct neuron_state * state;
  struct neuron_params * params;
};

struct network {
  long size;
  struct neuron ** neurons;
};

struct neuron_link {
  long to;
  float weight;
  float conduction_time;
};

struct link_queue {
  struct link_node * head;
  struct link_node * tail;
};

struct link_node {
  long from;
  long to;
  float weight;
  float ctime;
  struct link_node * next;
};

struct network * create_network(char * filename);
struct neuron * create_neuron();
struct neuron_params * init_neuron_params(char * filename);
struct neuron_state * init_init_neuron_state(char * filename);
void init_network_states(struct network * network, struct neuron_state * init_neuron_state);
struct neuron_state * copy_neuron_state(struct neuron_state * init_neuron_state);
void assoc_network_params(struct network * network, struct neuron_params * params);
void destroy_network_params(struct neuron_params * params);
void destroy_init_neuron_state(struct neuron_state * init_neuron_state);
void destroy_network(struct network * network);
void cleanup(struct network * network, struct neuron_state * init_neuron_state, struct neuron_params * params);
FILE * open_file_to_section(char * filename, char * section);
void remove_newline(char * line);
void link_neurons(struct network * network, char * filename);
struct neuron_link * create_link(long to, float weight, float ctime);
void queue_link(struct link_queue * link_queue, long from, long to, float weight, float ctime);
void create_queued_links(struct network * network, struct link_queue * link_queue);
void output_state(struct network * network, struct neuron_state * state, struct neuron_params * params, char * filename);
void write_to_file(FILE * fp, char * line);
void init_nondefault_states(struct network * network, char * filename);

#endif