#include "definitions.h"
#include "includes.h"
#include "math_includes.h"
#include "neuron.h"
#include "ode.h"
#include "currents.h"
#include "stimulate.h"
#include "plasticity.h"

#ifdef THREADED
void * ode_update_neurons_threaded(void * thread_params)
{
    struct thread_params * params = (struct thread_params *)thread_params;
    ode_update_neurons(params->network, params->start, params->num, params->y, params->f, params->t);
    pthread_exit(NULL);
}
#endif

void ode_update_neurons(struct network * network, long start, long num, const double * y, double * f, double t)
{
  double C_m, I_e, i_m;
  double i_L, i_Kdr, i_A, i_KCa, i_CaL,i_Na, i_NMDA, i_AMPA, i_in, i_coup, i_CaT;
  long i, j, k, l, limit = start + num;
  double coupling_factor;
  long neuron_offset, compartment_offset, link_offset, from_neuron, from_compartment;
  struct neuron_params * network_params;
  struct stimulus * stimulus;

  if(network->size < limit)
    limit = network->size;

  // get general network parameters that are the same for each neuron
  network_params = network->neurons[start]->params;
  C_m = network_params->values[0];
  coupling_factor = network_params->values[28];

  for(i = start; i < limit; i++)
    for(j = 0; j < network->compartments; j++)
      {
	compartment_offset = network->neurons[i]->compartments[j]->ode_system_offset;
	for(k = 0; k < network->neurons[i]->compartments[j]->state->num_params; k++)
	  network->neurons[i]->compartments[j]->state->values[k] = y[compartment_offset + k];   
	for(k = 0; k < network->neurons[i]->compartments[j]->num_links; k++)
	  {
	    link_offset = network->neurons[i]->compartments[j]->links[k]->ode_system_offset;
	    for(l = 0; l < network->neurons[i]->compartments[j]->links[k]->state->num_params; l++)
	      network->neurons[i]->compartments[j]->links[k]->state->values[l] = y[link_offset + l];
	  }
      }
  
  for(i = start; i < limit; i++)
    {
      neuron_offset = network->neurons[i]->ode_system_offset;
      for(j = 0; j < network->compartments; j++)
	{
	  compartment_offset = network->neurons[i]->compartments[j]->ode_system_offset;
	  stimulus = NULL;
	  I_e = network_params->values[26];
	  if((stimulus = apply_stimulus(network, i, j, t)) != NULL)
	    {
	      if(stimulus->direct == 1)
		I_e = stimulus->current;
	    }

	  // figure out individual currents
	  i_Na = Na_current(network,i,j,f,y,t);
	  i_Kdr = Kdr_current(network,i,j,f,y,t);
	  i_A = A_current(network,i,j,f,y,t);
	  i_KCa = KCa_current(network,i,j,f,y,t);
	  i_CaL = CaL_current(network,i,j,f,y,t);
	  i_CaT = CaT_current(network,i,j,f,y,t);
	  i_L = L_current(network,i,j,f,y,t);
	  i_NMDA = NMDA_current(network,i,j,f,y,t,i_CaL + i_CaT);
	  //i_NMDA = NMDA_current(network,i,j,f,y,t,i_CaL);
	  i_AMPA = AMPA_current(network,i,j,f,y,t);
	  diffuse_calcium(network,i,j,f,y,t);

	  // compartment 0 is the spine
	  if(j == 0)
	    {
	      i_in = 1.0*(i_NMDA + i_AMPA);
	      //i_coup = 1.0*coupling_factor*(y[offset + num_state_params] - y[offset]);
	      i_coup = coupling_factor * (y[network->neurons[i]->compartments[1]->ode_system_offset] - y[compartment_offset]);
	    }
	  else
	    {
	      i_in = I_e;
	      //i_coup = 1.0*coupling_factor*(y[offset - num_state_params] - y[offset]);
	      i_coup = coupling_factor * (y[network->neurons[i]->compartments[0]->ode_system_offset] - y[compartment_offset]);
	    }

	  
	  //debug
	  network->neurons[i]->compartments[j]->buffer->values[0] = i_Na;
	  network->neurons[i]->compartments[j]->buffer->values[1] = i_Kdr;
	  network->neurons[i]->compartments[j]->buffer->values[2] = i_A;
	  network->neurons[i]->compartments[j]->buffer->values[3] = i_KCa;
	  network->neurons[i]->compartments[j]->buffer->values[4] = i_CaL;
	  network->neurons[i]->compartments[j]->buffer->values[5] = i_L;
	  network->neurons[i]->compartments[j]->buffer->values[6] = i_NMDA;
	  network->neurons[i]->compartments[j]->buffer->values[7] = i_AMPA;
	  network->neurons[i]->compartments[j]->buffer->values[8] = i_in;
	  network->neurons[i]->compartments[j]->buffer->values[9] = i_coup;
	  

	  // membrane current
	  i_m = i_L + i_Kdr + i_A + i_KCa + i_CaL + i_Na + i_in + i_coup;
	  
	  // update derivatives (first state variable is voltage)
	  if(network->neurons[i]->compartments[j]->flag == FALSE)
	    {
	      if(y[compartment_offset] > 0)
		{
		  network->neurons[i]->compartments[j]->flag = TRUE;
		  network->neurons[i]->compartments[j]->spike_count++;
		}
	    }
	  else
	    {
	      if(y[compartment_offset] < -10)
		network->neurons[i]->compartments[j]->flag = FALSE;
	    }
	  f[compartment_offset] = i_m/C_m;
	  
	  // keep conductances constant
	  //for(k = 1; k < 10; k++)
	  //f[compartment_offset + k] = 0.0;

	  // detector system (compartment 0 is the spine)
	  if(j == 0)
	    {
	      for(k = 0; k < network->neurons[i]->compartments[j]->num_links; k++)
		{
		  from_neuron = network->neurons[i]->compartments[j]->links[k]->from;
		  from_compartment = network->neurons[i]->compartments[j]->links[k]->from_compartment;
		  if(network->neurons[from_neuron]->compartments[from_compartment]->flag == TRUE)
		    {
		      network->neurons[i]->compartments[j]->links[k]->recently_fired = TRUE;
		      network->neurons[i]->compartments[j]->links[k]->last_fired = t;
		    }

		  if(abs(network->neurons[i]->compartments[j]->links[k]->last_fired - t) > WINDOW)
		    network->neurons[i]->compartments[j]->links[k]->recently_fired = FALSE;

		  link_offset = network->neurons[i]->compartments[j]->links[k]->ode_system_offset;
		  f[link_offset + 0] = evolve_P(network,i,j,k,y);
		  f[link_offset + 1] = evolve_V(network,i,j,k,y);
		  f[link_offset + 2] = evolve_A(network,i,j,k,y);
		  f[link_offset + 3] = evolve_B(network,i,j,k,y);
		  f[link_offset + 4] = evolve_D(network,i,j,k,y);
		  f[link_offset + 5] = evolve_W(network,i,j,k,y);
		}
	    }
	}
    }
}

long get_ode_system_dimension(struct network * network)
{
  long i = 0, j = 0, k = 0, total = 0;
  for(i = 0; i < network->size; i++)
    {
      for(j = 0; j < network->compartments; j++)
	{
	  total += network->neurons[i]->compartments[j]->state->num_params;
	  for(k = 0; k < network->neurons[i]->compartments[j]->num_links; k++)
	      total += network->neurons[i]->compartments[j]->links[k]->state->num_params;
	}
    }

  return total;
}

int ode_run(struct network * network, double t, double t1, double step_size, double error)
{
  double end_runtime = t1;

  //const gsl_odeiv_step_type * T = gsl_odeiv_step_rk8pd;
  const gsl_odeiv_step_type * T = gsl_odeiv_step_rkf45;
  long i, j, k, l, dimension = get_ode_system_dimension(network), counter = 0;
  int status;
  
  gsl_odeiv_step * s = gsl_odeiv_step_alloc(T, dimension);
  gsl_odeiv_control * c = gsl_odeiv_control_y_new(error, error);
  gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc(dimension);
  
  // params: function, [jacobian], dimension, void * params
  #ifdef THREADED
    gsl_odeiv_system sys = {hh_ode_threaded, NULL, dimension, network};
  #else
    gsl_odeiv_system sys = {hh_ode, NULL, dimension, network};
  #endif

  double y[dimension];
  for(i = 0; i < network->size; i++)
    {
      network->neurons[i]->ode_system_offset = counter;
      for(j = 0; j < network->compartments; j++)
	{
	  network->neurons[i]->compartments[j]->ode_system_offset = counter;
	  for(k = 0; k < network->neurons[i]->compartments[j]->state->num_params; k++)
	    {
	      y[counter] = network->neurons[i]->compartments[j]->state->values[k];
	      counter++;
	    }
	  for(k = 0; k < network->neurons[i]->compartments[j]->num_links; k++)
	    {
	      network->neurons[i]->compartments[j]->links[k]->ode_system_offset = counter;
	      for(l = 0; l < network->neurons[i]->compartments[j]->links[k]->state->num_params; l++)
		{
		  y[counter] = network->neurons[i]->compartments[j]->links[k]->state->values[l];
		  counter++;
		}
	    }
	}
    }

  if(network->num_discontinuities != 0)
    {
      t1 = network->discontinuities[0];
      for(i = 0; i < network->num_discontinuities; i++)
	if(network->discontinuities[i] > end_runtime)
	  network->discontinuities[i] = end_runtime;
    }

  while (t < t1)
    {
      status = gsl_odeiv_evolve_apply(e, c, s, &sys, &t, t1, &step_size, y);
      if(status != GSL_SUCCESS)
	break;

      output_data(network,t,y);

      if(t >= t1 && network->num_discontinuities != 0)
	{
	  if(network->num_discontinuities - network->passed_discontinuities != 1)
	    {
	      network->passed_discontinuities++;
	      t1 = network->discontinuities[network->passed_discontinuities];
	    }
	  else
	    {
	      t1 = end_runtime;
	    }
	}
    }
  gsl_odeiv_evolve_free(e);
  gsl_odeiv_control_free(c);
  gsl_odeiv_step_free(s);
  
  return 0;
}

void output_data(struct network * network, double t, const double * y)
{
  int print_spine_v = 1;
  int print_spine_ca = 0;
  int print_spine_gating = 0;
  int print_spine_currents = 0;
  int print_soma_v = 1;
  int print_soma_ca = 0;
  int print_soma_gating = 0;
  int print_soma_currents = 0;
  int print_link_plasticity_all = 0;
  int print_link_plasticity_w = 0;
  int print_link_ca = 0;
  int print_link_gating = 0;
  long print_after_t = -1;
  long print_before_t = -1;
  
  long i, j, k, l, a, offset;
  
  if(t < print_after_t || (print_before_t != -1 && t > print_before_t))
    return;

  printf("%lf ", t);
  for(i = 0; i < network->size; i++)
    {
      for(j = 0; j < network->compartments; j++)
	{
	  if((print_spine_v == 1 && j == 0) || (print_soma_v == 1 && j == 1))
	    printf("%lf ",y[network->neurons[i]->compartments[j]->ode_system_offset]);
	  if((print_spine_ca == 1 && j == 0) || (print_soma_ca == 1 && j == 1))
	    printf("%lf ",y[network->neurons[i]->compartments[j]->ode_system_offset + 9]);
	  
	  if((print_spine_currents == 1 && j == 0) || (print_soma_currents == 1 && j == 1))
	    for(a = 0; a < 10; a++)
	      printf("%lf ",network->neurons[i]->compartments[j]->buffer->values[a]);
	  
	  if((print_spine_gating == 1 && j == 0) || (print_soma_gating == 1 && j == 1))
	    {
	      for(a = 1; a < 13; a++)
		if(a == 9 || (j == 1 && a == 12))
		  continue;
		else
		  printf("%lf ",y[network->neurons[i]->compartments[j]->ode_system_offset + a]);
	    }
	  
	  for(k = 0; k < network->neurons[i]->compartments[j]->num_links; k++)
	    {
	      if(print_link_ca == 1)
		printf("%lf ",y[network->neurons[i]->compartments[j]->links[k]->ode_system_offset + 12]);
	      if(print_link_plasticity_all == 1)
		for(l = 0; l < 6; l++)
		  {
		    offset = network->neurons[i]->compartments[j]->links[k]->ode_system_offset;
		    printf("%lf ",y[offset + l]);
		  }
	      else if(print_link_plasticity_w == 1)
		{
		  offset = network->neurons[i]->compartments[j]->links[k]->ode_system_offset;
		  printf("%lf",y[offset + 5]);
		}
	      if(print_link_gating == 1)
		for(l = 6; l < 12; l++)
		  {
		    offset = network->neurons[i]->compartments[j]->links[k]->ode_system_offset;
		    printf("%lf ",y[offset + l]);
		  }
	    }
	}
    }
  printf("\n");
}

int hh_ode(double t, const double y[], double f[], void *params)
{
  struct network * network = (struct network *)params;
  if(network->size < 1)
    {
      printf("[in ODE simulation] network size was less than one\n");
      exit(-1);
    }

  ode_update_neurons(network, 0, network->size, y, f, t);
  return GSL_SUCCESS;
}

#ifdef THREADED
int hh_ode_threaded(double t, const double y[], double f[], void *params)
{
  struct thread_params ** thread_params = (struct thread_params **)malloc(NUM_THREADS * sizeof(struct thread_params *));
  struct network * network = (struct network *)params;
  long i, neurons_per_thread, threads_spawned;
  pthread_t * threads;
  void * status;
  int rc;

  // make sure everything's there
  if(network->size < 1)
    {
      printf("[in ODE simulation] network size was less than one\n");
      exit(-1);
    }

  // make threads joinable
  //pthread_attr_init(&attr);
  //pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  neurons_per_thread = (long)network->size / NUM_THREADS;
  if(neurons_per_thread == 0)
    {
      neurons_per_thread = 1;
      threads_spawned = network->size;
    }
  else
    {
      threads_spawned = (long)network->size / neurons_per_thread + (network->size % neurons_per_thread == 0 ? 0 : 1);
    }
  threads = (pthread_t *)malloc(threads_spawned * sizeof(pthread_t));

  i = 0;
  while(i*neurons_per_thread < network->size)
    {
      thread_params[i] = (struct thread_params *)malloc(sizeof(struct thread_params));
      thread_params[i]->network = network;
      thread_params[i]->start = i * neurons_per_thread;
      thread_params[i]->num = neurons_per_thread;
      thread_params[i]->y = y;
      thread_params[i]->f = f;
      thread_params[i]->t = t;
      rc = pthread_create(&(threads[i]), NULL, ode_update_neurons_threaded, (void*)thread_params[i]);

      if(rc)
	{
	  printf("pthread_create() error: %d\n",rc);
	  exit(-1);
	}
      i++;
    }

  for(i = 0; i < threads_spawned; i++)
    {
      rc = pthread_join(threads[i], &status);
      if(rc)
	{
	  printf("pthread_join() error: %d\n",rc);
	  exit(-1);
	}
    }

  for(i = 0; i < threads_spawned; i++)
    free(thread_params[i]);
  free(thread_params);
  free(threads);
  return GSL_SUCCESS;
}
#endif
