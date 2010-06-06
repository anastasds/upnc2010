struct neuron {
  struct neuron_state * state;
  struct neuron_params * params;
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

struct network {
  long size;
  struct neuron ** neurons;
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
