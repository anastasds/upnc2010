#include "definitions.h"
#include "includes.h"
#include "math_includes.h"
#include "neuron.h"
#include "stimulate.h"

void stimulate(struct network * network, long num_neuron, long num_compartment, double start, double end, int direct, double current)
{
  if(start == end) return;
  struct stimulus * new_stimulus = (struct stimulus *)malloc(sizeof(struct stimulus));
  new_stimulus->from = start;
  new_stimulus->to = end;
  new_stimulus->direct = direct;
  new_stimulus->current = current;
  new_stimulus->next = NULL;

  append_stimulus(network->neurons[num_neuron]->compartments[num_compartment]->stimuli, new_stimulus);

  new_stimulus = (struct stimulus *)malloc(sizeof(struct stimulus));
  new_stimulus->from = start;
  new_stimulus->to = end;
  new_stimulus->direct = direct;
  new_stimulus->current = current;
  new_stimulus->next = NULL;

  append_stimulus(network->stimuli, new_stimulus);
}

void identify_discontinuities(struct network * network)
{
  long i = 0, uniques = 1, j = 0;
  if(network->stimuli == NULL || network->stimuli->length == 0) return;
  double *discontinuities,* times = (double *)malloc(network->stimuli->length * 2 * sizeof(double));
  struct stimulus * stim = network->stimuli->head;
  while(stim != NULL)
    {
      times[i] = stim->from;
      times[i+1] = stim->to;
      if(times[i] == 0.0) times[i] = times[i+1];
      stim = stim->next;
      i += 2;
    }
  qsort((void*)times,network->stimuli->length * 2,sizeof(double),compare_doubles);
  for(i = 0; i < network->stimuli->length * 2 - 1; i++)
      if(times[i] != times[i+1])
	uniques+=1;

  discontinuities = (double *)malloc(uniques * sizeof(double));
  discontinuities[0] = times[0];
  for(i = 0; j < uniques-1; i++)
    {
      if(times[i] != discontinuities[j])
	discontinuities[++j] = times[i];
    }
  free(times);
  network->discontinuities = discontinuities;
  network->num_discontinuities = uniques;
  network->passed_discontinuities = 0;

}

int compare_doubles(const void * a, const void * b)
{
  return (*((double*)a) <= *((double*)b)) ? -1 : 1;
}

void append_stimulus(struct stimuli * stimuli, struct stimulus * stimulus)
{
  stimuli->length++;
  if(stimuli->head == NULL)
    {
      stimuli->head = stimulus;
      stimuli->tail = stimulus;
    }
  else
    {
      stimuli->tail->next = stimulus;
      stimuli->tail = stimulus;
    }
}

void prepare_tetanus(struct network * network, char * filename)
{
  long num_neuron, num_compartment;
  double start, length, current, periodicity, i;
  int direct;
  char line[MAX_LINE_LEN];
  FILE * fp = open_file_to_section(filename,"@TETANUS");

  fgets(line, MAX_LINE_LEN, fp);
  remove_newline(line);

  while(strcmp(line,"") != 0 && strcmp(line,"0") != 0)
    {
      sscanf(line, "%ld %ld %d %lf %lf %lf %lf", &num_neuron, &num_compartment, &direct, &start, &length, &periodicity, &current);
      if(network->neurons[num_neuron]->compartments[num_compartment]->stimuli == NULL)
	network->neurons[num_neuron]->compartments[num_compartment]->stimuli = init_stimuli_struct();
      
      if(network->stimuli == NULL)
	network->stimuli = init_stimuli_struct();

      for(i = start; i < network->runtime[1]; i += periodicity)
	{
	  stimulate(network, num_neuron, num_compartment, i, i + length, direct, current);
	}
      fgets(line, MAX_LINE_LEN, fp);
      remove_newline(line);
    }
  fclose(fp);
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
      
      if(network->stimuli == NULL)
	network->stimuli = init_stimuli_struct();

      stimulate(network, num_neuron, num_compartment, start, end, direct, current);

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
  stimuli->length = 0;
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
