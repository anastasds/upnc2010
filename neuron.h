#ifndef NEURON_H
#define NEURON_H

struct thread_params {
  struct network * network;
  const double * y;
  double * f;
  double t;
  long start;
  long num;
};

struct buffer {
  int size;
  float * values;
};

struct neuron {
  struct neuron_compartment ** compartments;
  struct neuron_params * params;
};

struct neuron_compartment {
  struct neuron * neuron;
  struct compartment_state * state;
  struct stimuli * stimuli;
  struct buffer * buffer;
  int num_links;
  struct neuron_link ** links;
  long spike_count;
  int flag;
};

struct compartment_state {
  int num_params;
  char ** names;
  double * values;
};

struct init_compartment_states {
  long num_states;
  struct compartment_state ** states;
};

struct neuron_params {
  int num_params;
  char ** names;
  double * values;
};

struct network {
  long size;
  long compartments;
  long num_discontinuities;
  long passed_discontinuities;
  struct neuron ** neurons;
  struct stimuli * stimuli;
  double * discontinuities;
};

struct neuron_link {
  long from;
  long from_compartment;
  long to;
  long to_compartment;
  double weight;
  double conduction_time;
};

struct link_queue {
  struct link_node * head;
  struct link_node * tail;
};

struct link_node {
  long from;
  long from_compartment;
  long to;
  long to_compartment;
  double weight;
  double ctime;
  struct link_node * next;
};

struct stimulus {
  double from;
  double to;
  double current;
  int direct;
  struct stimulus * next;
};

struct stimuli {
  struct stimulus * head;
  struct stimulus * tail;
  long length;
};

struct network * create_network(char * filename);
void create_neuron_compartments(struct neuron * neuron, long n);
double * get_runtime(char * filename);
struct neuron * create_neuron();
struct neuron_params * init_neuron_params(char * filename);
struct compartment_state * init_init_compartment_state(long num_compartment, char * filename);
struct init_compartment_states * init_init_compartment_states(struct network * network, char * filename);
void init_network_states(struct network * network, struct init_compartment_states * init_compartment_states);
struct compartment_state * copy_compartment_state(struct compartment_state * init_compartment_state);
void assoc_network_params(struct network * network, struct neuron_params * params);
void destroy_network_params(struct neuron_params * params);
void destroy_init_compartment_states(struct init_compartment_states * init_compartment_states);
void destroy_network(struct network * network);
void cleanup(struct network * network, struct init_compartment_states * init_compartment_states, struct neuron_params * params);
FILE * open_file_to_section(char * filename, char * section);
char * remove_newline(char * line);
void link_neurons(struct network * network, char * filename);
struct neuron_link * create_link(long from, long from_compartment, long to, long to_compartment, double weight, double ctime);
void queue_link(struct link_queue * link_queue, long from, long from_compartment, long to, long to_compartment, double weight, double ctime);
void create_queued_links(struct network * network, struct link_queue * link_queue);
void output_state(struct network * network, struct init_compartment_states * states, struct neuron_params * params, char * filename);
void write_to_file(FILE * fp, char * line);
void init_nondefault_states(struct network * network, char * filename);
void print_network(struct network * network);
struct buffer * allocate_buffer(int size);
void print_buffer(struct buffer * buffer);

#endif
