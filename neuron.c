#include "definitions.h"
#include "includes.h"
#include "neuron.h"

// initializes struct network of size specified in input file;
// creates each neuron in the network using create_neuron();
// returns pointer to network
struct network * create_network(char * filename)
{
  long i;
  struct network * network = malloc(sizeof(struct network));
  FILE * fp = open_file_to_section(filename, "@NETWORK_SIZE");

  fscanf(fp, "%ld %ld", &(network->size),&(network->compartments));
  fclose(fp);

  if(DEBUG > 0)
    printf("* network->size set to %ld\n",network->size);

  network->neurons = malloc(network->size*sizeof(struct neuron *));
  for(i = 0; i < network->size; i++)
    {
      network->neurons[i] = create_neuron();
      create_neuron_compartments(network->neurons[i],network->compartments);
    }
  return network;
}

void create_neuron_compartments(struct neuron * neuron, long n)
{
  long i;
  neuron->compartments = (struct neuron_compartment **)malloc(n * sizeof(struct neuron_compartment *));
  for(i = 0; i < n; i++)
    {
      neuron->compartments[i] = (struct neuron_compartment *) malloc(sizeof(struct neuron_compartment));
      neuron->compartments[i]->neuron = neuron;
      neuron->compartments[i]->state = NULL;
    }
}

// malloc()'s a new struct neuron;
// initializes elements to 0 or NULL;
// returns pointer;
struct neuron * create_neuron()
{
  struct neuron * new_neuron = malloc(sizeof(struct neuron));
  new_neuron->num_links = 0;
  //new_neuron->state = NULL;
  new_neuron->params = NULL;
  new_neuron->links = NULL;
  return new_neuron;
}

// opens config file to @PARAMS;
// malloc()s struct network_params, populates it from config file;
// returns pointer to network parameters struct
struct neuron_params * init_neuron_params(char * filename)
{
  int i = 0;
  struct neuron_params * params = malloc(sizeof(struct neuron_params));
  FILE * fp = open_file_to_section(filename, "@PARAMS");
  char line[MAX_LINE_LEN];

  fscanf(fp, "%d",&(params->num_params));
  if(DEBUG > 0)
    printf("* params->num_params set to %d\n", params->num_params);
  params->names = malloc(params->num_params * sizeof(char *));
  params->values = malloc(params->num_params * sizeof(float));
  
  for(i = 0; i < params->num_params; i++)
    {
      fscanf(fp, "%s", line);
      fscanf(fp, "%f", &params->values[i]);
      params->names[i] = malloc((strlen(line) + 1) * sizeof(char));
      strcpy(params->names[i], line);
      if(DEBUG > 0)
	  printf("* params->names[%d]: %s = %f\n", i, params->names[i], params->values[i]);
    }
  fclose(fp);
  if(DEBUG > 0)
    printf("* init_neuron_params successful\n");
  return params;
}

// malloc()s a neuron_state with state variables set to values
// defined to be default in @DEFAULT_STATE section of config file;
// returns pointer to neuron_state()
struct compartment_state * init_init_compartment_state(char * filename)
{
  int i = 0;
  struct compartment_state * state = malloc(sizeof(struct compartment_state));
  FILE * fp = open_file_to_section(filename, "@DEFAULT_STATE");
  char line[MAX_LINE_LEN];

  fscanf(fp, "%d",&(state->num_params));
  if(DEBUG > 0)
    printf("* init_compartment_state->num_params set to %d\n", state->num_params);
  state->names = malloc(state->num_params * sizeof(char *));
  state->values = malloc(state->num_params * sizeof(float));
  
  for(i = 0; i < state->num_params; i++)
    {
      fscanf(fp, "%s", line);
      fscanf(fp, "%f", &state->values[i]);
      state->names[i] = malloc((strlen(line) + 1) * sizeof(char));
      strcpy(state->names[i], line);
      if(DEBUG > 0)
	  printf("* state->names[%d]: %s = %f\n", i, state->names[i], state->values[i]);
    }
  fclose(fp);
  if(DEBUG > 0)
    printf("* init_init_compartment_state successful\n");
  return state;
}

// copies default neuron states as defined in init_init_neuron_state()
// to each neurons in the network
void init_network_states(struct network * network, struct compartment_state * init_compartment_state)
{
  long i, j;
  for(i = 0; i < network->size; i++)
    for(j = 0; j < network->compartments; j++)
      network->neurons[i]->compartments[j]->state = copy_compartment_state(init_compartment_state);
}

// parses @INIT_STATES to update network neurons to initial state variable values that
// differ from the default values set in @DEFAULT_STATE of config file
void init_nondefault_states(struct network * network, char * filename)
{
  FILE * fp;
  int num_params, num_read;
  long i, n, c;
  float param;
  char *tmp, *line = malloc(MAX_LINE_LEN * sizeof(char));;
  tmp = line;

  if(network->size == 0)
    {
      printf("init_nondefault_states(): network has size 0\n");
      exit(-1);
    }
  if(network->neurons == NULL)
    {
      printf("init_nondefault_states() called with null network->neurons despite nonzero network size\n");
      exit(-1);
    }
  if(network->neurons[0] == NULL)
    {
      printf("init_nondefault_states() called with null network->neurons[0] despite nonzero network size\n");
      exit(-1);
    }
  if(network->neurons[0]->compartments[0]->state == NULL)
    {
      printf("init_nondefault_states() called before init_network_states, apparently\n");
      exit(-1);
    }

  num_params = network->neurons[0]->compartments[0]->state->num_params;

  fp = open_file_to_section(filename, "@INIT_STATES");
  fgets(line, MAX_LINE_LEN, fp);
  remove_newline(line);
  while(strcmp(line,"") != 0 && strcmp(line,"0") != 0)
    {
      sscanf(line, "%ld %ld%n", &n, &c, &num_read);
      if(n < 0 || n > network->size)
	{
	  printf("@INIT_STATES tried to set neuron %ld's state despite network size of %ld\n",n, network->size);
	  exit(-1);
	}
      if(c < 0 || c > network->compartments)
	{
	  printf("@INIT_STATES tried to set compartment %ld of neuron %ld despite compartment number of %ld\n", c, n, network->compartments);
	  exit(-1);
	}
      if(network->neurons[n]->compartments[c]->state == NULL)
	{
	  printf("init_nondefault_states(): neuron %ld compartment %ld has NULL state\n", n, c);
	  exit(-1);
	}
      if(network->neurons[n]->compartments[c]->state->values == NULL)
	{
	  printf("init_nondefault_states(): neuron %ld compartment %ld has NULL state->values\n", n, c);
	  exit(-1);
	}
      line += num_read;
      for(i = 0; i < num_params; i++)
	{
	  sscanf(line, "%f%n", &param, &num_read);
	  network->neurons[n]->compartments[c]->state->values[i] = param;
	  line += num_read;
	}
      line = tmp;
      fgets(line, MAX_LINE_LEN, fp);
      remove_newline(line);
    }
  fclose(fp);
  free(line);
}

// malloc()s a new struct neuron_state;
// copies all struct elements from given neuron_state to the new one;
// returns pointer to new neuron_state
struct compartment_state * copy_compartment_state(struct compartment_state * init_compartment_state)
{
  int i = 0;
  struct compartment_state * new_state = malloc(sizeof(struct compartment_state));

  new_state->num_params = init_compartment_state->num_params;
  new_state->names = malloc(new_state->num_params * sizeof(char *));
  new_state->values = malloc(new_state->num_params * sizeof(float));
  for(i = 0; i < new_state->num_params; i++)
    {
      new_state->names[i] = malloc((strlen(init_compartment_state->names[i]) + 1) * sizeof(char));
      strcpy(new_state->names[i],init_compartment_state->names[i]);
      new_state->values[i] = init_compartment_state->values[i];
    }
  return new_state;
}

// links all network neurons to a single neuron_params struct;
// this means that neuron_params values are constant throughout the network
void assoc_network_params(struct network * network, struct neuron_params * params)
{
  long i;
  for(i = 0; i < network->size; i++)
    network->neurons[i]->params = params;
}

// cleanup function
void destroy_network_params(struct neuron_params * params)
{
  int i;
  for(i = 0; i < params->num_params; i++)
    free(params->names[i]);
  free(params->names);
  free(params->values);
  free(params);
}

// cleanup function
void destroy_init_compartment_state(struct compartment_state * init_compartment_state)
{
  int i;
  for(i = 0; i < init_compartment_state->num_params; i++)
    free(init_compartment_state->names[i]);
  free(init_compartment_state->names);
  free(init_compartment_state->values);
  free(init_compartment_state);
}

// cleanup function
void destroy_network(struct network * network)
{
  long i, j, k;
  for(i = 0; i < network->size; i++)
    {
      for(k = 0; k < network->compartments; k++)
	{
	  for(j = 0; j < network->neurons[i]->compartments[k]->state->num_params; j++)
	    free(network->neurons[i]->compartments[k]->state->names[j]);
	  free(network->neurons[i]->compartments[k]->state->names);
	  free(network->neurons[i]->compartments[k]->state->values);
	  free(network->neurons[i]->compartments[k]->state);
	  free(network->neurons[i]->compartments[k]);
	}
      free(network->neurons[i]->compartments);
      if(network->neurons[i]->num_links > 0)
	{
	  for(j = 0; j < network->neurons[i]->num_links; j++)
	    free(network->neurons[i]->links[j]);
	  free(network->neurons[i]->links);
	}
      free(network->neurons[i]);
    }
  free(network->neurons);
  free(network);
}

// cleanup function
void cleanup(struct network * network, struct compartment_state * init_compartment_state, struct neuron_params * params)
{
  destroy_network_params(params);
  destroy_init_compartment_state(init_compartment_state);
  destroy_network(network);
}

// opens given file, moves file pointer to the line following the
// string specified in char*section, and returns file pointer
FILE * open_file_to_section(char * filename, char * section)
{
  FILE * fp;
  char * line = malloc(MAX_LINE_LEN * sizeof(char));

  if((fp = fopen(filename, "rt")) == NULL)
    {
      printf("Error opening file %s.\n", filename);
      exit(-1);
    }

  if(DEBUG > 1)
    printf("* open_file_to_section(%s, %s)\n",filename, section);
  
  while(fgets(line, MAX_LINE_LEN - 1, fp) != NULL)
    {
      remove_newline(line);
      if(DEBUG > 1)
	printf("* looking for \"%s\", found \"%s\"\n", section, line);
      if(strcmp(line, section) == 0)
	{
	  free(line);
	  if(DEBUG > 1)
	    printf("* open_file_to_section successful\n");
	  return fp;
	}
    }
  
  free(line);
  fclose(fp);
  printf("Error finding section %s in file %s.\n", section, filename);
  exit(-1);
}

// strips trailing \n from given string if it in offset MAX_LINE_LEN or less;
// char *s that are given to this function generally contain data read from a file,
// so they are generally allocated of max size MAX_LINE_LEN, so this should be fine
void remove_newline(char * line)
{
  int i = 0;
  while(*(line + i++) != '\n');
  if(i-1 <= MAX_LINE_LEN)
    line[i-1] = '\0';
}

// parses config file's @LINKS section to create links between neurons
void link_neurons(struct network * network, char * filename)
{
  int random_weights, random_ctimes;
  long i, j, from, from_compartment, to, to_compartment;
  float weight, ctime, p_connected, r, default_weight, default_ctime;
  char line[MAX_LINE_LEN];

  if(DEBUG > 0)
    printf("* link_neurons\n");

  FILE * fp = open_file_to_section(filename, "@LINKS");
  
  if(fgets(line, MAX_LINE_LEN, fp) == NULL)
    {
      printf("Error reading first line after @LINKS.\n");
      exit(-1);
    }

  remove_newline(line);
  if(strcmp(line, "total") == 0)
    {
      if(DEBUG > 0)
	printf("* totally connected network selected\n");
      if(fgets(line, MAX_LINE_LEN, fp) == NULL)
	{
	  printf("Error reading second line after @LINKS.\n");
	  exit(-1);
	}
      
      sscanf(line, "%f %f", &weight, &ctime);
      if(DEBUG > 0)
	printf("* link weights %f, conductance times %f\n", weight, ctime);

      for(i = 0; i < network->size; i++)
	{
	  network->neurons[i]->num_links = network->size - 1;
	  network->neurons[i]->links = malloc((network->size - 1) * sizeof(struct neuron_link *));
	  for(j = 0; j < i; j++)
	    network->neurons[i]->links[j] = create_link(i, network->compartments - 1, j, 0, weight, ctime);

	  for(j = i+1; j < network->size; j++)
	    network->neurons[i]->links[j-1] = create_link(i, network->compartments - 1, j, 0, weight, ctime);
	}
    }
  else if(strcmp(line,"random") == 0)
    {
      srand(time(NULL));
      fscanf(fp, "%f %d %f %d %f", &p_connected, &random_weights, &default_weight, &random_ctimes, &default_ctime);

      struct link_queue * link_queue = malloc(sizeof(struct link_queue));
      link_queue->head = NULL;
      link_queue->tail = NULL;

      for(i = 0; i < network->size; i++)
	{
	  for(j = 0; j < network->size; j++)
	    {
	      if(i == j)
		continue;

	      if(random_weights == 0)
		weight = default_weight;
	      else
		weight = ((float)(rand() % 1000)) / 1000.0;

	      if(random_ctimes == 0)
		ctime = default_ctime;
	      else
		ctime = ((float)(rand() % 1000)) / 1000.0;

	      r = ((float)(rand() % 1000)) / 1000.0;
	      if(r < p_connected)
		{
		  network->neurons[i]->num_links++;
		  queue_link(link_queue, i, network->compartments - 1, j, 0, weight, ctime);
		}
	    }
	}
      create_queued_links(network, link_queue);
    }
  else if(strcmp(line,"defined") == 0)
    {
      if(DEBUG > 0)
	printf("* specifically defined network selected\n");

      while(fgets(line, MAX_LINE_LEN, fp) != NULL && strlen(line) > 0)
	{
	  sscanf(line, "%ld %ld", &from, &from_compartment);
	  if(DEBUG > 1)
	    printf("found link from neuron %ld compartment %ld\n",from, from_compartment);
	  network->neurons[from]->num_links++;
	}
      fclose(fp);

      int * num_created = malloc(network->size * sizeof(int));
      memset(num_created, 0, network->size * sizeof(int));
      for(i = 0; i < network->size; i++)
	{
	  if(network->neurons[i]->num_links > 0)
	    network->neurons[i]->links = (struct neuron_link **)malloc(network->neurons[i]->num_links * sizeof(struct neuron_link *));
	}

      fp = open_file_to_section(filename, "@LINKS");
      fgets(line, MAX_LINE_LEN, fp);
      while(fgets(line, MAX_LINE_LEN, fp) != NULL)
	{
	  sscanf(line, "%ld %ld %ld %ld %f %f", &from, &from_compartment, &to, &to_compartment, &weight, &ctime);
	  network->neurons[from]->links[num_created[from]++] = create_link(from, from_compartment, to, to_compartment, weight, ctime);
	}

      for(i = 0; i < network->size; i++)
	{
	  if(network->neurons[i]->num_links != num_created[i])
	    {
	      printf("Error occurred while creating network links.\n");
	      exit(-1);
	    }
	}
      free(num_created);
    }
  else
    {
      printf("Unknown network architecture chosen: %s\n",line);
      exit(-1);
    }
  fclose(fp);
}

// maloc()s a struct neuron_link;
// populates it with given data;
// returns pointer to new struct
struct neuron_link * create_link(long from, long from_compartment, long to, long to_compartment, float weight, float ctime)
{
  struct neuron_link * link = malloc(sizeof(struct neuron_link));
  link->from = from;
  link->from_compartment = from_compartment;
  link->to = to;
  link->to_compartment = to_compartment;
  link->weight = weight;
  link->conduction_time = ctime;

  if(DEBUG > 1)
    printf("* link: neuron %ld compartment %ld to neuron %ld compartment %ld, wieght %f, ctime %f\n", from, from_compartment, to, to_compartment, weight, ctime);

  return link;
}

// neuron structs have a "links" pointer and a num_links param, but the num_links cannot be
// determined and used to allocate the links struct without having read the entire @LINKS
// section; therefore, a linked list of links to be created is made as @LINKS is parsed,
// and this list is later gone through using create_queued_links
void queue_link(struct link_queue * link_queue, long from, long from_compartment, long to, long to_compartment, float weight, float ctime)
{
  struct link_node * node = malloc(sizeof(struct link_node));
  node->from = from;
  node->from_compartment = from_compartment;
  node->to = to;
  node->to_compartment = to_compartment;
  node->weight = weight;
  node->ctime = ctime;
  node->next = NULL;

  if(link_queue->tail == NULL)
    {
      link_queue->head = node;
      link_queue->tail = node;
    }
  else
    {
      link_queue->tail->next = node;
      link_queue->tail = node;
    }  
}

// goes through the linked list of links to create and actually makes them, i.e.
// updates the neuron structs appropriately
void create_queued_links(struct network * network, struct link_queue * link_queue)
{
  struct link_node * node, *old_node;
  node = link_queue->head;
  long i;

  int * num_created = malloc(network->size * sizeof(int));
  memset(num_created, 0, network->size * sizeof(int));

  for(i = 0; i < network->size; i++)
    {
      if(network->neurons[i]->num_links > 0)
	{
	  network->neurons[i]->links = malloc(network->neurons[i]->num_links * sizeof(struct neuron_link *));
	}
    }

  while(node != NULL)
    {
      network->neurons[node->from]->links[num_created[node->from]++] = create_link(node->from, node->from_compartment, node->to, node->to_compartment, node->weight, node->ctime);
      old_node = node;
      node = node->next;
      free(old_node);
      link_queue->head = node;
    }

  for(i = 0; i < network->size; i++)
    {
      if(network->neurons[i]->num_links != num_created[i])
	{
	  printf("Error creating links.\n");
	  exit(-1);
	}
    }

  free(link_queue);
  free(num_created);
}

// takes the current state of the network and outputs it into a file (with filename given) formatted
// as a config file so that we can analyze the data or use it to progress the simulation later
void output_state(struct network * network, struct compartment_state * state, struct neuron_params * params, char * filename)
{
  FILE * fp;
  long i,j,k;
  char tmp[20];
  char line[MAX_LINE_LEN];

  if((fp = fopen(filename, "wt")) == NULL)
    {
      printf("Could not open output file for writing.\n");
      exit(-1);
    }

  sprintf(line, "@NETWORK_SIZE\n%ld %ld\n\n", network->size, network->compartments);
  write_to_file(fp, line);

  sprintf(line, "@PARAMS\n%d\n", params->num_params);
  write_to_file(fp, line);

  for(i = 0; i < params->num_params; i++)
    {
      sprintf(line,"%s %f\n", params->names[i], params->values[i]);
      write_to_file(fp, line);
    }

  sprintf(line, "\n@DEFAULT_STATE\n%d\n", state->num_params);
  write_to_file(fp, line);

  for(i = 0; i < state->num_params; i++)
    {
      sprintf(line,"%s %f\n", state->names[i], state->values[i]);
      write_to_file(fp, line);
    }

  sprintf(line, "\n@INIT_STATES\n");
  write_to_file(fp, line);

  for(i = 0; i < network->size; i++)
    {
      for(k = 0; k < network->compartments; k++)
	{
	  sprintf(line, "%ld %ld", i, k);
	  for(j = 0; j < state->num_params; j++)
	    {
	      sprintf(tmp, " %f", network->neurons[i]->compartments[k]->state->values[j]);
	      strcat(line, tmp);
	    }
	  strcat(line, "\n");
	  write_to_file(fp,line);
	}
    }

  sprintf(line, "\n@LINKS\ndefined\n");
  write_to_file(fp, line);

  for(i = 0; i < network->size; i++)
    {
      for(j = 0; j < network->neurons[i]->num_links; j++)
	{
	  sprintf(line, "%ld %ld %ld %ld %ld %f %f\n", i, network->neurons[i]->links[j]->from, network->neurons[i]->links[j]->from_compartment, network->neurons[i]->links[j]->to, network->neurons[i]->links[j]->to_compartment, network->neurons[i]->links[j]->weight, network->neurons[i]->links[j]->conduction_time);
	  write_to_file(fp, line);
	}
    }
  fclose(fp);

}

// self-explanatory
void write_to_file(FILE * fp, char * line)
{
  fwrite(line, sizeof(char), strlen(line), fp);
}
