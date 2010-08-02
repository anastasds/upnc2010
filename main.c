#include "definitions.h"
#include "includes.h"
#include "math_includes.h"
#include "general.h"
#include "neuron.h"
#include "ode.h"

int main(int argc, char** argv)
{
  char * input_filename = get_input_filename(argc, argv);
  struct network * network = create_network(input_filename);
  struct neuron_params * params = init_neuron_params(input_filename);
  struct init_compartment_states * init_compartment_states = init_init_compartment_states(network, input_filename);
  
  init_network_states(network, init_compartment_states);
  init_nondefault_states(network, input_filename);
  assoc_network_params(network, params);
  link_neurons(network, input_filename);

  ode_run(network, 0, 50.0, 1.0e-6, 1.0e-6);

  //print_network(network);

  if(argc > 2) output_state(network, init_compartment_states, params, argv[2]);
  else if(argc < 2) free(input_filename);
  cleanup(network, init_compartment_states, params);

  return 0;
}
