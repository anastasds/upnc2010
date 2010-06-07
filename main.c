#define NETWORK_SIZE 1000
#define MAX_LINE_LEN 80
#define DEBUG 2

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "neuron.h"
#include "neuron.c"

int main(int argc, char** argv)
{
  char filename[] = "template.txt";
  struct network * network = create_network(filename);
  struct neuron_params * params = init_neuron_params(filename);
  struct neuron_state * init_neuron_state = init_init_neuron_state(filename);
  
  init_network_states(network, init_neuron_state);
  assoc_network_params(network, params);

  link_neurons(network, filename);

  cleanup(network, init_neuron_state, params);

  return 0;
}
