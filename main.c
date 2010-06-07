#define MAX_LINE_LEN 80
#define DEBUG 0
#define DEFAULT_CONFIG_FILE "template.txt"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "neuron.h"
#include "neuron.c"

int main(int argc, char** argv)
{
  char * input_filename;

  if(argc < 2)
    {
      input_filename = malloc((strlen(DEFAULT_CONFIG_FILE) + 1) * sizeof(char));
      strcpy(input_filename, DEFAULT_CONFIG_FILE);
    }
  else
    input_filename = argv[1];

  struct network * network = create_network(input_filename);
  struct neuron_params * params = init_neuron_params(input_filename);
  struct neuron_state * init_neuron_state = init_init_neuron_state(input_filename);
  
  init_network_states(network, init_neuron_state);
  assoc_network_params(network, params);
  link_neurons(network, input_filename);

  if(argc > 2)
    output_state(network, init_neuron_state, params, argv[2]);

  cleanup(network, init_neuron_state, params);
  if(argc < 2)
    free(input_filename);

  return 0;
}
