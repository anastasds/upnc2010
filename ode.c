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
  double i_L, i_Kdr, i_A, i_KCa, i_CaT, i_Na, i_NMDA, i_AMPA, i_in, i_coup;
  long i, j, k, offset, num_state_params, limit = start + num;
  double coupling_factor;
  struct neuron_params * network_params;
  struct stimulus * stimulus;

  if(network->size < limit)
    limit = network->size;

  // get general network parameters that are the same for each neuron
  network_params = network->neurons[start]->params;
  C_m = network_params->values[0];
  coupling_factor = network_params->values[28];

  num_state_params = network->neurons[start]->compartments[0]->state->num_params;

  for(i = start; i < limit; i++)
    for(j = 0; j < network->compartments; j++)
      {
	offset = num_state_params * network->compartments * i + num_state_params * j;
	for(k = 0; k < network->neurons[i]->compartments[j]->state->num_params; k++)
	  network->neurons[i]->compartments[j]->state->values[k] = y[offset + k];   
      }
  
  for(i = start; i < limit; i++)
    {
      for(j = 0; j < network->compartments; j++)
	{
	  stimulus = NULL;
	  I_e = network_params->values[26];
	  if((stimulus = apply_stimulus(network, i, j, t)) != NULL)
	    {
	      if(stimulus->direct == 1)
		I_e = stimulus->current;
	    }

	  offset = num_state_params * network->compartments * i + num_state_params * j;

	  // figure out individual currents
	  i_Na = Na_current(network,i,j,f,y,t);
	  i_Kdr = Kdr_current(network,i,j,f,y,t);
	  i_A = A_current(network,i,j,f,y,t);
	  i_KCa = KCa_current(network,i,j,f,y,t);
	  i_CaT = CaT_current(network,i,j,f,y,t);
	  i_L = L_current(network,i,j,f,y,t);
	  i_NMDA = NMDA_current(network,i,j,f,y,t,i_CaT);
	  i_AMPA = AMPA_current(network,i,j,f,y,t);

	  // compartment 0 is the spine
	  if(j == 0)
	    {
	      i_in = 1.0*(i_NMDA + i_AMPA);
	      i_coup = 1.0*coupling_factor*(y[offset + num_state_params] - y[offset]);
	    }
	  else
	    {
	      i_in = I_e;
	      i_coup = 1.0*coupling_factor*(y[offset - num_state_params] - y[offset]);
	    }

	  // membrane current
	  i_m = i_L + i_Kdr + i_A + i_KCa + i_CaT + i_Na + i_in + i_coup;
	  
	  // update derivatives (first state variable is voltage)
	  if(network->neurons[i]->compartments[j]->flag == -1)
	    {
	      if(y[offset] > 10)
		{
		  network->neurons[i]->compartments[j]->flag = 1;
		  network->neurons[i]->compartments[j]->spike_count++;
		}
	    }
	  else
	    {
	      if(y[offset] < -10)
		network->neurons[i]->compartments[j]->flag = -1;
	    }
	  f[offset] = i_m/C_m;
	  
	  // keep conductances constant
	  for(k = 1; k < 10; k++)
	    f[k] = 0.0;

	  // detector system (compartment 0 is the spine)
	  if(j == 0)
	    {
	      f[offset + 26] = evolve_P(network,y,offset);
	      f[offset + 27] = evolve_V(network,y,offset);
	      f[offset + 28] = evolve_A(network,y,offset);
	      f[offset + 29] = evolve_B(network,y,offset);
	      f[offset + 30] = evolve_D(network,y,offset);
	      f[offset + 31] = evolve_W(network,y,offset);
	    }
	  else
	    {
	      for(k = 0; k < 6; k++)
		{
		  f[offset + 26 + k] = 0.0;
		}
	    }
	}
    }
}

int ode_run(struct network * network, double t, double t1, double step_size, double error)
{
  const gsl_odeiv_step_type * T = gsl_odeiv_step_rk8pd;
  long i, j, k, num_state_params = network->neurons[0]->compartments[0]->state->num_params;
  long dimension = network->size * network->compartments * num_state_params;
  int status;
  
  gsl_odeiv_step * s = gsl_odeiv_step_alloc(T, dimension);
  gsl_odeiv_control * c = gsl_odeiv_control_y_new(error, 0.0);
  gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc(dimension);
  
  // params: function, [jacobian], dimension, void * params
  #ifdef THREADED
    gsl_odeiv_system sys = {hh_ode_threaded, NULL, dimension, network};
  #else
    gsl_odeiv_system sys = {hh_ode, NULL, dimension, network};
  #endif

  double y[num_state_params * network->size * network->compartments];
  for(i = 0; i < network->size; i++)
    for(j = 0; j < network->compartments; j++)
      for(k = 0; k < num_state_params; k++)
	y[num_state_params * network->compartments * i + num_state_params * j + k] = network->neurons[i]->compartments[j]->state->values[k];

  while (t < t1)
    {
      status = gsl_odeiv_evolve_apply(e, c, s, &sys, &t, t1, &step_size, y);
      if(status != GSL_SUCCESS)
	break;
                
      printf("%lf ", t);
      for(i = 0; i < network->size * network->compartments; i++)
	{
	  // print voltage, [Ca]
	  printf("%lf %lf ",y[num_state_params * i], y[num_state_params * i + 18]);

	  // print synaptic plasticity values
	  for(j = 0; j < 6; j++)
	    printf("%lf ",y[num_state_params*i + 26 + j]);

	  // print change in synaptic weight
	  //printf("%lf ",y[num_state_params*i + 31]);
	}
      printf("\n");
      
    }
  
  gsl_odeiv_evolve_free(e);
  gsl_odeiv_control_free(c);
  gsl_odeiv_step_free(s);
  
  /*
  for(i = 0; i < network->size; i++)
    for(j = 0; j < network->compartments; j++)
      printf("%ld %ld %ld\n",i,j,network->neurons[i]->compartments[j]->spike_count);
  */

  return 0;
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
