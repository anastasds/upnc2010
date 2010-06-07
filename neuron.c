struct network * create_network(char * filename)
{
  long i;
  struct network * network = malloc(sizeof(struct network));
  FILE * fp = open_file_to_section(filename, "@NETWORK_SIZE");

  fscanf(fp, "%ld", &(network->size));
  fclose(fp);

  if(DEBUG > 0)
    printf("* network->size set to %ld\n",network->size);

  network->neurons = malloc(network->size*sizeof(struct neuron *));
  for(i = 0; i < network->size; i++)
      network->neurons[i] = create_neuron();
  return network;
}

struct neuron * create_neuron()
{
  struct neuron * new_neuron = malloc(sizeof(struct neuron));
  new_neuron->num_links = 0;
  new_neuron->links = NULL;
  return new_neuron;
}

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
  printf("* init_neuron_params successful\n");
  return params;
}

struct neuron_state * init_init_neuron_state(char * filename)
{
  int i = 0;
  struct neuron_state * state = malloc(sizeof(struct neuron_state));
  FILE * fp = open_file_to_section(filename, "@DEFAULT_STATE");
  char line[MAX_LINE_LEN];

  fscanf(fp, "%d",&(state->num_params));
  if(DEBUG > 0)
    printf("* init_neuron_state->num_params set to %d\n", state->num_params);
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
  printf("* init_init_neuron_state successful\n");
  return state;
}

void init_network_states(struct network * network, struct neuron_state * init_neuron_state)
{
  long i;
  for(i = 0; i < network->size; i++)
    network->neurons[i]->state = copy_neuron_state(init_neuron_state);
}

struct neuron_state * copy_neuron_state(struct neuron_state * init_neuron_state)
{
  int i = 0;
  struct neuron_state * new_state = malloc(sizeof(struct neuron_state));

  new_state->num_params = init_neuron_state->num_params;
  new_state->names = malloc(new_state->num_params * sizeof(char *));
  new_state->values = malloc(new_state->num_params * sizeof(float));
  for(i = 0; i < new_state->num_params; i++)
    {
      new_state->names[i] = malloc((strlen(init_neuron_state->names[i]) + 1) * sizeof(char));
      strcpy(new_state->names[i],init_neuron_state->names[i]);
      new_state->values[i] = init_neuron_state->values[i];
    }
  return new_state;
}

void assoc_network_params(struct network * network, struct neuron_params * params)
{
  long i;
  for(i = 0; i < network->size; i++)
    network->neurons[i]->params = params;
}

void destroy_network_params(struct neuron_params * params)
{
  int i;
  for(i = 0; i < params->num_params; i++)
    free(params->names[i]);
  free(params->names);
  free(params->values);
  free(params);
}

void destroy_init_neuron_state(struct neuron_state * init_neuron_state)
{
  int i;
  for(i = 0; i < init_neuron_state->num_params; i++)
    free(init_neuron_state->names[i]);
  free(init_neuron_state->names);
  free(init_neuron_state->values);
  free(init_neuron_state);
}

void destroy_network(struct network * network)
{
  long i, j;
  for(i = 0; i < network->size; i++)
    {
      for(j = 0; j < network->neurons[i]->state->num_params; j++)
	  free(network->neurons[i]->state->names[j]);
      free(network->neurons[i]->state->names);
      free(network->neurons[i]->state->values);
      free(network->neurons[i]->state);
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

void cleanup(struct network * network, struct neuron_state * init_neuron_state, struct neuron_params * params)
{
  destroy_network_params(params);
  destroy_init_neuron_state(init_neuron_state);
  destroy_network(network);
}

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

void remove_newline(char * line)
{
  int i = 0;
  while(*(line + i++) != '\n');
  if(i-1 <= MAX_LINE_LEN)
    line[i-1] = '\0';
}

void link_neurons(struct network * network, char * filename)
{
  int random_weights, random_ctimes;
  long i, j, from, to;
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
	      network->neurons[i]->links[j] = create_link(j, weight, ctime);

	  for(j = i+1; j < network->size; j++)
	      network->neurons[i]->links[j-1] = create_link(j, weight, ctime);
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
		  queue_link(link_queue, i, j, weight, ctime);
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
	  sscanf(line, "%ld", &from);
	  if(DEBUG > 1)
	    printf("found link from neuron #%ld\n",from);
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
	  sscanf(line, "%ld %ld %f %f", &from, &to, &weight, &ctime);
	  network->neurons[from]->links[num_created[from]++] = create_link(to, weight, ctime);
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

struct neuron_link * create_link(long to, float weight, float ctime)
{
  struct neuron_link * link = malloc(sizeof(struct neuron_link));
  link->to = to;
  link->weight = weight;
  link->conduction_time = ctime;

  if(DEBUG > 1)
    printf("* link to #%ld, wieght %f, ctime %f\n", to, weight, ctime);

  return link;
}

void queue_link(struct link_queue * link_queue, long from, long to, float weight, float ctime)
{
  struct link_node * node = malloc(sizeof(struct link_node));
  node->from = from;
  node->to = to;
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
      network->neurons[node->from]->links[num_created[node->from]++] = create_link(node->to, node->weight, node->ctime);
      old_node = node;
      node = node->next;
      free(old_node);
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
