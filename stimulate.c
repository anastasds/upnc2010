#include "definitions.h"
#include "includes.h"
#include "math_includes.h"
#include "neuron.h"
#include "stimulate.h"

void stimulate(struct stimuli * stimuli, double start, double end, int direct, double current)
{
  struct stimulus * new_stimulus = (struct stimulus *)malloc(sizeof(struct stimulus));
  new_stimulus->from = start;
  new_stimulus->to = end;
  new_stimulus->direct = direct;
  new_stimulus->current = current;
  new_stimulus->next = NULL;

  if(stimuli->head == NULL)
    {
      stimuli->head = new_stimulus;
      stimuli->tail = new_stimulus;
    }
  else
    {
      stimuli->tail->next = new_stimulus;
      stimuli->tail = new_stimulus;
    }
}

void init_stimuli(struct network * network, char * filename)
{
  long num_neuron, num_compartment;
  double start, end, current;
  int direct;
  char line[MAX_LINE_LEN];
  FILE * fp = open_file_to_section(filename,"@STIMULI");
  
  fgets(line, MAX_LINE_LEN, fp);
  remove_newline(line);

  while(strcmp(line,"") != 0 && strcmp(line,"0") != 0)
    {
      sscanf(line, "%ld %ld %d %lf %lf %lf", &num_neuron, &num_compartment, &direct, &start, &end, &current);
      if(num_neuron < 0 || num_neuron > network->size - 1)
	{
	  printf("init_stimuli(): config file says to stimulate neuron that doesn't exist (#%ld, network size is %ld\n",num_neuron,network->size);
	  exit(-1);
	}

      if(num_compartment < 0 || num_compartment > network->compartments - 1)
	{
	  printf("init_stimuli(): config file says to stimulate compartment that doesn't exist (#%ld, neuron has %ld compartments\n",num_compartment,network->compartments);
	  exit(-1);
	}

      if(network->neurons[num_neuron]->compartments[num_compartment]->stimuli == NULL)
	network->neurons[num_neuron]->compartments[num_compartment]->stimuli = init_stimuli_struct();
      
      stimulate(network->neurons[num_neuron]->compartments[num_compartment]->stimuli, start, end, direct, current);

      fgets(line, MAX_LINE_LEN, fp);
      remove_newline(line);
    }

  fclose(fp);
}

struct stimuli * init_stimuli_struct()
{
  struct stimuli * stimuli = (struct stimuli *)malloc(sizeof(struct stimuli));
  stimuli->head = NULL;
  stimuli->tail = NULL;
  return stimuli;
}

struct stimulus * apply_stimulus(struct network * network, long num_neuron, long num_compartment, double t)
{
  if(network->neurons[num_neuron]->compartments[num_compartment]->stimuli == NULL)
    return NULL;
  struct stimulus * stimulus = network->neurons[num_neuron]->compartments[num_compartment]->stimuli->head;
  if(stimulus == NULL)
    return NULL;
  while(stimulus != NULL)
    {
      if(t > stimulus->from && t < stimulus->to)
	return stimulus;
      stimulus = stimulus->next;
    }
  return NULL;
}

void destroy_stimuli(struct stimuli * stimuli)
{
  if(stimuli == NULL)
    return;
  struct stimulus * tmp;
  while((tmp = stimuli->head) != NULL)
    {
      stimuli->head = tmp->next;
      free(tmp);
    }
  free(stimuli);
}
