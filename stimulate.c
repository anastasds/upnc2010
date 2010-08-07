#include "definitions.h"
#include "includes.h"
#include "neuron.h"
#include "stimulate.h"

void stimulate(struct stimuli * stimuli, long num_neuron, long num_compartment, double start, double end, double current)
{
  struct stimulus * new_stimulus = (struct stimulus *)malloc(sizeof(struct stimulus));
  new_stimulus->from = start;
  new_stimulus->to = end;
  new_stimulus->neuron = num_neuron;
  new_stimulus->compartment = num_compartment;
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
  struct stimuli * stimuli = init_stimuli_struct();
  long num_neuron, num_compartment;
  double start, end, current;
  char line[MAX_LINE_LEN];
  FILE * fp = open_file_to_section(filename,"@STIMULI");
  
  fgets(line, MAX_LINE_LEN, fp);
  remove_newline(line);

  while(strcmp(line,"") != 0 && strcmp(line,"0") != 0)
    {
      sscanf(line, "%ld %ld %lf %lf %lf", &num_neuron, &num_compartment, &start, &end, &current);
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

      stimulate(stimuli, num_neuron, num_compartment, start, end, current);
      network->neurons[num_neuron]->compartments[num_compartment]->stimulated = TRUE;

      fgets(line, MAX_LINE_LEN, fp);
      remove_newline(line);
    }

  fclose(fp);

  network->stimuli = stimuli;
}

struct stimuli * init_stimuli_struct()
{
  struct stimuli * stimuli = (struct stimuli *)malloc(sizeof(struct stimuli));
  stimuli->head = NULL;
  stimuli->tail = NULL;
  return stimuli;
}

double apply_stimulus(struct network * network, long num_neuron, long num_compartment, double t)
{
  struct stimulus * stimulus = network->stimuli->head;
  while(1)
    {
      stimulus = find_stimulus(stimulus, num_neuron, num_compartment);
      if(stimulus == NULL)
	return 0.0;
      if(t > stimulus->from && t < stimulus->to)
	return stimulus->current;
      stimulus = stimulus->next;
    }

}

struct stimulus * find_stimulus(struct stimulus * stimulus, long num_neuron, long num_compartment)
{
  while(stimulus != NULL && (stimulus->neuron != num_neuron || stimulus->compartment != num_compartment))
    {
      stimulus = stimulus->next;
    }
  return stimulus;
}

void destroy_stimuli(struct stimuli * stimuli)
{
  struct stimulus * tmp;
  while((tmp = stimuli->head) != NULL)
    {
      stimuli->head = tmp->next;
      free(tmp);
    }
  free(stimuli);
}
