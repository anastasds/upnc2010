#include "definitions.h"
#include "includes.h"
#include "math_includes.h"

int main(int argc, char** argv)
{
  char * input_filename = get_input_filename(argc, argv);
  struct network * network = create_network(input_filename);
  struct neuron_params * params = init_neuron_params(input_filename);
  struct neuron_state * init_neuron_state = init_init_neuron_state(input_filename);
  
  init_network_states(network, init_neuron_state);
  init_nondefault_states(network, input_filename);
  assoc_network_params(network, params);
  link_neurons(network, input_filename);

  ode_test(network->neurons[0]->state,params);

  if(argc > 2) output_state(network, init_neuron_state, params, argv[2]);
  else if(argc < 2) free(input_filename);
  cleanup(network, init_neuron_state, params);

  return 0;
}