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
