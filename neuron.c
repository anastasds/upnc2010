#include "definitions.h"
#include "includes.h"
#include "neuron.h"
#include "stimulate.h"

// initializes struct network of size specified in input file;
// creates each neuron in the network using create_neuron();
// returns pointer to network
struct network * create_network(char * filename)
{
  long i;
  struct network * network = (struct network *)malloc(sizeof(struct network));
  FILE * fp = open_file_to_section(filename, "@NETWORK_SIZE");

  fscanf(fp, "%ld %ld", &(network->size),&(network->compartments));
  fclose(fp);

  if(DEBUG > 0)
    printf("* network->size set to %ld, %ld compartments each\n",network->size, network->compartments);

  network->neurons = (struct neuron **)malloc(network->size*sizeof(struct neuron *));
  network->stimuli = NULL;
  network->num_discontinuities = 0;
  network->discontinuities = NULL;
  network->passed_discontinuities = 0;

  for(i = 0; i < network->size; i++)
    {
      network->neurons[i] = create_neuron();
      create_neuron_compartments(network->neurons[i],network->compartments);
    }

  return network;
}

double * get_runtime(char * filename)
{
  double * runtime = (double *)malloc(2*sizeof(double));
  FILE * fp = open_file_to_section(filename,"@RUNTIME");

  fscanf(fp,"%lf %lf", &runtime[0], &runtime[1]);
  fclose(fp);
  return runtime;
}

void create_neuron_compartments(struct neuron * neuron, long n)
{
  long i;

  if(n < 1)
    {
      printf("create_neuron_compartments(): number of compartments per neuron was %ld, must be poitive\n",n);
      exit(-1);
    }

  neuron->compartments = (struct neuron_compartment **)malloc(n * sizeof(struct neuron_compartment *));
  for(i = 0; i < n; i++)
    {
      neuron->compartments[i] = (struct neuron_compartment *) malloc(sizeof(struct neuron_compartment));
      neuron->compartments[i]->neuron = neuron;
      neuron->compartments[i]->buffer = allocate_buffer(BUFFER_SIZE);
      neuron->compartments[i]->state = NULL;
      neuron->compartments[i]->links = NULL;
      neuron->compartments[i]->num_links = 0;
      neuron->compartments[i]->stimuli = NULL;
      neuron->compartments[i]->spike_count = 0;
      neuron->compartments[i]->flag = -1;
    }
}

struct buffer * allocate_buffer(int size)
{
  if(size == 0)
    return NULL;

  int i;
  struct buffer * buf = (struct buffer *) malloc(sizeof(struct buffer));
  buf->size = size;
  buf->values = (float *) malloc(size * sizeof(float));
  for(i = 0; i < size; i++)
    buf->values[i] = 0.0;
  return buf;
}

void print_buffer(struct buffer * buffer)
{
  if(buffer == NULL)
    return;

  int i;
  for(i = 0; i < buffer->size; i++)
    printf("%lf ",buffer->values[i]);
}

// malloc()'s a new struct neuron;
// initializes elements to 0 or NULL;
// returns pointer;
struct neuron * create_neuron()
{
  struct neuron * new_neuron = (struct neuron *)malloc(sizeof(struct neuron));
  new_neuron->compartments = NULL;
  new_neuron->params = NULL;
  return new_neuron;
}

// opens config file to @PARAMS;
// malloc()s struct network_params, populates it from config file;
// returns pointer to network parameters struct
struct neuron_params * init_neuron_params(char * filename)
{
  int i = 0;
  struct neuron_params * params = (struct neuron_params *)malloc(sizeof(struct neuron_params));
  FILE * fp = open_file_to_section(filename, "@PARAMS");
  char line[MAX_LINE_LEN];

  fscanf(fp, "%d",&(params->num_params));
  if(DEBUG > 0)
    printf("* params->num_params set to %d\n", params->num_params);
  params->names = (char **)malloc(params->num_params * sizeof(char *));
  params->values = (double *)malloc(params->num_params * sizeof(double));
  
  for(i = 0; i < params->num_params; i++)
    {
      fscanf(fp, "%s", line);
      fscanf(fp, "%lf", &params->values[i]);
      params->names[i] = (char *)malloc((strlen(line) + 1) * sizeof(char));
      strcpy(params->names[i], line);
      if(DEBUG > 0)
	  printf("* params->names[%d]: %s = %lf\n", i, params->names[i], params->values[i]);
    }
  fclose(fp);
  if(DEBUG > 0)
    printf("* init_neuron_params successful\n");
  return params;
}

// malloc()s a neuron_state with state variables set to values
// defined to be default in @DEFAULT_STATE section of config file;
// returns pointer to neuron_state()
struct compartment_state * init_init_compartment_state(long num_compartment, char * filename)
{
  int i = 0;
  struct compartment_state * state = (struct compartment_state *)malloc(sizeof(struct compartment_state));
  FILE * fp;
  char line[MAX_LINE_LEN];

  sprintf(line, "@DEFAULT_STATE %ld", num_compartment);
  fp = open_file_to_section(filename, line);
  fscanf(fp, "%d",&(state->num_params));
  if(DEBUG > 0)
    printf("* init_compartment_state->num_params set to %d\n", state->num_params);
  state->names = (char **)malloc(state->num_params * sizeof(char *));
  state->values = (double *)malloc(state->num_params * sizeof(double));
  
  for(i = 0; i < state->num_params; i++)
    {
      fscanf(fp, "%s", line);
      fscanf(fp, "%lf", &state->values[i]);
      state->names[i] = (char *)malloc((strlen(line) + 1) * sizeof(char));
      strcpy(state->names[i], line);
      if(DEBUG > 0)
	  printf("* state->names[%d]: %s = %lf\n", i, state->names[i], state->values[i]);
    }
  fclose(fp);
  if(DEBUG > 0)
    printf("* init_init_compartment_state successful\n");
  return state;
}

struct init_compartment_states * init_init_compartment_states(struct network * network, char * filename)
{
  long i;
  struct init_compartment_states * new_states = (struct init_compartment_states *)malloc(sizeof(struct init_compartment_states));

  new_states->num_states = network->compartments;
  new_states->states = (struct compartment_state **)malloc(network->compartments * sizeof(struct compartment_state));
  for(i = 0; i < network->compartments; i++)
    new_states->states[i] = init_init_compartment_state(i, filename);

  return new_states;
}

// copies default neuron states as defined in init_init_neuron_state()
// to each neurons in the network
void init_network_states(struct network * network, struct init_compartment_states * init_compartment_states)
{
  long i, j;
  for(i = 0; i < network->size; i++)
    for(j = 0; j < network->compartments; j++)
      network->neurons[i]->compartments[j]->state = copy_compartment_state(init_compartment_states->states[j]);
}

// parses @INIT_STATES to update network neurons to initial state variable values that
// differ from the default values set in @DEFAULT_STATE of config file
void init_nondefault_states(struct network * network, char * filename)
{
  FILE * fp;
  int num_params, num_read;
  long i, n, c;
  double param;
  char *tmp, *line = (char *)malloc(MAX_LINE_LEN * sizeof(char));;
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
	  sscanf(line, "%lf%n", &param, &num_read);
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
  struct compartment_state * new_state = (struct compartment_state *)malloc(sizeof(struct compartment_state));

  new_state->num_params = init_compartment_state->num_params;
  new_state->names = (char **)malloc(new_state->num_params * sizeof(char *));
  new_state->values = (double *)malloc(new_state->num_params * sizeof(double));
  for(i = 0; i < new_state->num_params; i++)
    {
      new_state->names[i] = (char *)malloc((strlen(init_compartment_state->names[i]) + 1) * sizeof(char));
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
void destroy_init_compartment_states(struct init_compartment_states * init_compartment_states)
{
  long i,j;
  for(i = 0; i < init_compartment_states->num_states; i++)
    {
      for(j = 0; j < init_compartment_states->states[i]->num_params; j++)
	{
	  free(init_compartment_states->states[i]->names[j]);
	}
      free(init_compartment_states->states[i]->names);
      free(init_compartment_states->states[i]->values);
      free(init_compartment_states->states[i]);
    }
  free(init_compartment_states->states);
  free(init_compartment_states);
}

// cleanup function
void destroy_network(struct network * network)
{
  long i, j, k, m;

  for(i = 0; i < network->size; i++)
    {
      for(k = 0; k < network->compartments; k++)
	{
	  for(j = 0; j < network->neurons[i]->compartments[k]->state->num_params; j++)
	    free(network->neurons[i]->compartments[k]->state->names[j]);
	  free(network->neurons[i]->compartments[k]->state->names);
	  free(network->neurons[i]->compartments[k]->state->values);
	  free(network->neurons[i]->compartments[k]->state);
	  if(network->neurons[i]->compartments[k]->buffer != NULL && network->neurons[i]->compartments[k]->buffer->values != NULL)
	    {
	      free(network->neurons[i]->compartments[k]->buffer->values);
	      free(network->neurons[i]->compartments[k]->buffer);
	    }
	  destroy_stimuli(network->neurons[i]->compartments[k]->stimuli);
	  if(network->neurons[i]->compartments[k]->num_links > 0)
	    {
	      for(m = 0; m < network->neurons[i]->compartments[k]->num_links; m++)
		free(network->neurons[i]->compartments[k]->links[m]);
	      free(network->neurons[i]->compartments[k]->links);
	    }
	  free(network->neurons[i]->compartments[k]);
	}
      free(network->neurons[i]->compartments);
      free(network->neurons[i]);
    }
  free(network->neurons);
  destroy_stimuli(network->stimuli);
  if(network->num_discontinuities != 0)
    free(network->discontinuities);
  free(network->runtime);
  free(network);
}

// cleanup function
void cleanup(struct network * network, struct init_compartment_states * init_compartment_states, struct neuron_params * params)
{
  destroy_network_params(params);
  destroy_init_compartment_states(init_compartment_states);
  destroy_network(network);
}

// opens given file, moves file pointer to the line following the
// string specified in char*section, and returns file pointer
FILE * open_file_to_section(char * filename, char * section)
{
  FILE * fp;
  char * line = (char *)malloc(MAX_LINE_LEN * sizeof(char));

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
char * remove_newline(char * line)
{
  int i = 0;
  while(*(line + i++) != '\n');
  if(i-1 <= MAX_LINE_LEN)
    line[i-1] = '\0';
  return line;
}

// parses config file's @LINKS section to create links between neurons
void link_neurons(struct network * network, char * filename)
{
  int random_weights, random_ctimes;
  long i, j, from, from_compartment, to, to_compartment, links_found = 0;
  double weight, ctime, p_connected, r, default_weight, default_ctime;
  char line[MAX_LINE_LEN];
  int * num_created;

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
      
      sscanf(line, "%lf %lf", &weight, &ctime);
      if(DEBUG > 0)
	printf("* link weights %lf, conductance times %lf\n", weight, ctime);

      for(i = 0; i < network->size; i++)
	{
	  // for now, we're assuming that all presynaptic activity is received on compartment 0,
	  // and that all synaptic output comes from compartment 1
	  network->neurons[i]->compartments[0]->num_links = network->size - 1;
	  network->neurons[i]->compartments[0]->links = (struct neuron_link **)malloc((network->size - 1) * sizeof(struct neuron_link *));
	  for(j = 0; j < i; j++)
	    network->neurons[i]->compartments[0]->links[j] = create_link(i, network->compartments - 1, j, 0, weight, ctime);
	  
	  for(j = i+1; j < network->size; j++)
	    network->neurons[i]->compartments[0]->links[j-1] = create_link(i, network->compartments - 1, j, 0, weight, ctime);
	}
    }
  else if(strcmp(line,"random") == 0)
    {
      srand(time(NULL));
      fscanf(fp, "%lf %d %lf %d %lf", &p_connected, &random_weights, &default_weight, &random_ctimes, &default_ctime);

      struct link_queue * link_queue = (struct link_queue *)malloc(sizeof(struct link_queue));
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
		weight = ((double)(rand() % 1000)) / 1000.0;

	      if(random_ctimes == 0)
		ctime = default_ctime;
	      else
		ctime = ((double)(rand() % 1000)) / 1000.0;

	      r = ((double)(rand() % 1000)) / 1000.0;
	      if(r < p_connected)
		{
		  network->neurons[j]->compartments[0]->num_links++;
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

      while(fgets(line, MAX_LINE_LEN, fp) != NULL && strlen(line) > 1)
	{
	  sscanf(line, "%ld %ld %ld %ld", &from, &from_compartment, &to, &to_compartment);
	  if(DEBUG > 1)
	    printf("found link to neuron %ld compartment %ld from neuron %ld compartment %ld\n",to, to_compartment, from, from_compartment);
	  network->neurons[to]->compartments[to_compartment]->num_links++;
	  links_found++;
	}

      if(links_found > 0)
	{
	  num_created = (int *)malloc(network->size * network->compartments * sizeof(int));
	  memset(num_created, 0, network->size * network->compartments * sizeof(int));
	  for(i = 0; i < network->size; i++)
	    for(j = 0; j < network->compartments; j++)
	      {
		if(network->neurons[i]->compartments[j]->num_links > 0)
		  network->neurons[i]->compartments[j]->links = (struct neuron_link **)malloc(network->neurons[i]->compartments[j]->num_links * sizeof(struct neuron_link *));
	      }

	  fclose(fp);	 
	  fp = open_file_to_section(filename, "@LINKS");
	  fgets(line, MAX_LINE_LEN, fp);
	  while(fgets(line, MAX_LINE_LEN, fp) != NULL && strlen(line) > 1)
	    {
	      sscanf(line, "%ld %ld %ld %ld %lf %lf", &from, &from_compartment, &to, &to_compartment, &weight, &ctime);
	      network->neurons[to]->compartments[to_compartment]->links[num_created[network->compartments * to + to_compartment]++] = create_link(from, from_compartment, to, to_compartment, weight, ctime);
	    }
	  
	  for(i = 0; i < network->size; i++)
	    for(j = 0; j < network->compartments; j++)
	      if(network->neurons[i]->compartments[j]->num_links != num_created[network->compartments * i + j])
		{
		  printf("Error occurred while creating network links.\n");
		  exit(-1);
		}
	  free(num_created);
	}
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
struct neuron_link * create_link(long from, long from_compartment, long to, long to_compartment, double weight, double ctime)
{
  struct neuron_link * link = (struct neuron_link *)malloc(sizeof(struct neuron_link));
  link->from = from;
  link->from_compartment = from_compartment;
  link->to = to;
  link->to_compartment = to_compartment;
  link->weight = weight;
  link->conduction_time = ctime;

  if(DEBUG > 1)
    printf("* link: neuron %ld compartment %ld to neuron %ld compartment %ld, wieght %lf, ctime %lf\n", from, from_compartment, to, to_compartment, weight, ctime);

  return link;
}

// neuron structs have a "links" pointer and a num_links param, but the num_links cannot be
// determined and used to allocate the links struct without having read the entire @LINKS
// section; therefore, a linked list of links to be created is made as @LINKS is parsed,
// and this list is later gone through using create_queued_links
void queue_link(struct link_queue * link_queue, long from, long from_compartment, long to, long to_compartment, double weight, double ctime)
{
  struct link_node * node = (struct link_node *)malloc(sizeof(struct link_node));
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
  long i,j;

  int * num_created = (int *)malloc(network->size * sizeof(int));
  memset(num_created, 0, network->size * sizeof(int));

  for(i = 0; i < network->size; i++)
    {
      for(j = 0; j < network->compartments; j++)
	{
	  if(network->neurons[i]->compartments[j]->num_links > 0)
	    {
	      network->neurons[i]->compartments[j]->links = (struct neuron_link **)malloc(network->neurons[i]->compartments[j]->num_links * sizeof(struct neuron_link *));
	    }
	}
    }

  while(node != NULL)
    {
      network->neurons[node->to]->compartments[node->to_compartment]->links[num_created[node->to * network->compartments + node->to_compartment]++] = create_link(node->from, node->from_compartment, node->to, node->to_compartment, node->weight, node->ctime);
      old_node = node;
      node = node->next;
      free(old_node);
      link_queue->head = node;
    }

  for(i = 0; i < network->size; i++)
    {
      for(j = 0; i < network->compartments; j++)
	{
	  if(network->neurons[i]->compartments[j]->num_links != num_created[network->compartments * i + j])
	    {
	      printf("Error creating links.\n");
	      exit(-1);
	    }
	}
    }

  free(link_queue);
  free(num_created);
}

// takes the current state of the network and outputs it into a file (with filename given) formatted
// as a config file so that we can analyze the data or use it to progress the simulation later
void output_state(struct network * network, struct init_compartment_states * states, struct neuron_params * params, char * filename)
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
      sprintf(line,"%s %lf\n", params->names[i], params->values[i]);
      write_to_file(fp, line);
    }

  for(i = 0; i < states->num_states; i++)
    {
      sprintf(line, "\n@DEFAULT_STATE %ld\n%d\n", i, states->states[i]->num_params);
      write_to_file(fp, line);

      for(j = 0; j < states->states[i]->num_params; j++)
	{
	  sprintf(line,"%s %lf\n", states->states[i]->names[j], states->states[i]->values[j]);
	  write_to_file(fp, line);
	}
    }

  sprintf(line, "\n@INIT_STATES\n");
  write_to_file(fp, line);

  for(i = 0; i < network->size; i++)
    {
      for(k = 0; k < network->compartments; k++)
	{
	  sprintf(line, "%ld %ld", i, k);
	  //for(j = 0; j < states->states[k]->num_params; j++)
	  for(j = 0; j < network->neurons[i]->compartments[k]->state->num_params; j++)
	    {
	      sprintf(tmp, " %lf", network->neurons[i]->compartments[k]->state->values[j]);
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
      for(j = 0; j < network->compartments; j++)
	{
	  for(k = 0; k < network->neurons[i]->compartments[j]->num_links; k++)
	    {
	      sprintf(line, "%ld %ld %ld %ld %lf %lf\n", network->neurons[i]->compartments[j]->links[k]->from, network->neurons[i]->compartments[j]->links[k]->from_compartment, network->neurons[i]->compartments[j]->links[k]->to, network->neurons[i]->compartments[j]->links[k]->to_compartment, network->neurons[i]->compartments[j]->links[k]->weight, network->neurons[i]->compartments[j]->links[k]->conduction_time);
	      write_to_file(fp, line);
	    }
	}
    }
  fclose(fp);

}

// self-explanatory
void write_to_file(FILE * fp, char * line)
{
  fwrite(line, sizeof(char), strlen(line), fp);
}

// self-explanatory
void print_network(struct network * network)
{
  long i,j,k;
  for(i = 0; i < network->size; i++)
    {
      printf("\n=================\n\n* * NEURON %ld\n", i);
      for(j = 0; j < network->compartments; j++)
	{
	  printf("\n* COMPARTMENT %ld\n", j);
	  for(k = 0; k < network->neurons[i]->compartments[j]->state->num_params; k++)
	    {
	      printf("* %s = %lf\n",network->neurons[i]->compartments[j]->state->names[k], network->neurons[i]->compartments[j]->state->values[k]);
	    }
	}
    }
}
