#include "definitions.h"
#include "includes.h"
#include "math_includes.h"
#include "neuron.h"
#include "ode.h"

void * ode_update_neurons_threaded(void * thread_params)
{
    struct thread_params * params = (struct thread_params *)thread_params;
    ode_update_neurons(params->network, params->start, params->num, params->y, params->f);
}

void ode_update_neurons(struct network * network, long start, long num, const double * y, double * f)
{
  double C_m, g_bar_Na, g_bar_K, g_bar_L, E_Na, E_K, E_L, I_e;
  double alpha_n, beta_n, alpha_m, beta_m, alpha_h, beta_h, i_m;
  double tau_n, tau_m, tau_h, n_inf, m_inf, h_inf;
  struct neuron_params * network_params;
  long i, limit = start + num;

  if(network->size < limit)
    limit = network->size;

  // get general network parameters that are the same for each neuron
  network_params = network->neurons[0]->params;
  C_m = network_params->values[0];
  g_bar_Na = network_params->values[1];
  g_bar_K = network_params->values[2];
  g_bar_L = network_params->values[3];
  E_Na = network_params->values[4];
  E_K = network_params->values[5];
  E_L = network_params->values[6];
  I_e = network_params->values[7];

  for(i = start; i < limit; i++)
    {
      //update params in neuron_state
      network->neurons[i]->state->values[0] = y[4*i];
      network->neurons[i]->state->values[1] = y[4*i + 1];
      network->neurons[i]->state->values[2] = y[4*i + 2];
      network->neurons[i]->state->values[3] = y[4*i + 3];

      // find dn/dt, dm/dt, dh/dt usig some intermediate values
      alpha_n = 0.01*(y[4*i] + 55)/(1 - exp(-0.1*(y[4*i]+55.0)));
      beta_n = 0.125*exp(-0.0125*(y[4*i]+65.0));
      n_inf = alpha_n / (alpha_n + beta_n);
      tau_n = 1.0/(alpha_n + beta_n);
      
      alpha_m = 0.1*(y[4*i] + 40.0) / (1.0 - exp(-0.1*(y[4*i]+40.0)));
      beta_m = 4.0*exp(-0.0556*(y[4*i]+65.0));
      m_inf = alpha_m / (alpha_m + beta_m);
      tau_m = 1.0/(alpha_m + beta_m);
      
      alpha_h = 0.07*exp(-0.05*(y[4*i]+65.0));
      beta_h = 1.0 / (1.0 + exp(-0.1*(y[4*i]+35.0)));
      h_inf = alpha_h / (alpha_h + beta_h);
      tau_h = 1.0/(alpha_h + beta_h);
      
      // figure out membrane current, i_m
      i_m = g_bar_L*(y[4*i] - E_L) + g_bar_K*pow(y[4*i + 1],4.0)*(y[4*i] - E_K) + g_bar_Na*pow(y[4*i + 2],3)*y[4*i + 3]*(y[4*i] - E_Na);
      
      // update dx_i/dt, x_i \in {v,n,m,h}
      f[4*i] = I_e + -1.0*i_m/C_m;
      f[4*i + 1] = (n_inf - y[4*i + 1])/(tau_n);
      f[4*i + 2] = (m_inf - y[4*i + 2])/(tau_m);
      f[4*i + 3] = (h_inf - y[4*i + 3])/(tau_h);
    }
  #ifdef THREADED
    pthread_exit(NULL);
  #endif
}

int ode_run(struct network * network, double t, double t1, double step_size, double error)
{
  const gsl_odeiv_step_type * T = gsl_odeiv_step_rk8pd;
  long dimension = network->size * 4;
  long i = 0;
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

  // set up ode system, four odes per neuron
  // 0 = V, 1 = n, 2 = m, 3 = h
  double y[4*network->size];
  for(i = 0; i < network->size; i++)
    {
      y[4*i] = network->neurons[i]->state->values[0];
      y[4*i + 1] = network->neurons[i]->state->values[1];
      y[4*i + 2] = network->neurons[i]->state->values[2];
      y[4*i + 3] = network->neurons[i]->state->values[3];
    }

  // GSL IS EVOLVING!
  while (t < t1)
    {
      status = gsl_odeiv_evolve_apply(e, c, s, &sys, &t, t1, &step_size, y);
      if(status != GSL_SUCCESS)
	break;
    }
  // GSL HAS EVOLVED INTO... HODGKIN-HUXLEY!
  
  gsl_odeiv_evolve_free(e);
  gsl_odeiv_control_free(c);
  gsl_odeiv_step_free(s);

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

  ode_update_neurons(network, 0, network->size, y, f);
  return GSL_SUCCESS;
}

#ifdef THREADED
int hh_ode_threaded(double t, const double y[], double f[], void *params)
{
  struct thread_params ** thread_params = (struct thread_params **)malloc(NUM_THREADS * sizeof(struct thread_params *));
  struct network * network = (struct network *)params;
  long i, neurons_per_thread, limit, threads_spawned;
  pthread_t * threads;
  pthread_attr_t attr;
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
      thread_params[i]->thread_id = i;
      thread_params[i]->network = network;
      thread_params[i]->start = i * neurons_per_thread;
      thread_params[i]->num = neurons_per_thread;
      thread_params[i]->y = y;
      thread_params[i]->f = f;
      thread_params[i]->t = t;
      rc = pthread_create(&(threads[i]), NULL, ode_update_neurons_threaded, (void*)thread_params[i]);

      if(rc)
	{
	  printf("pthread_create() error: $d\n",rc);
	  exit(-1);
	}
      i++;
    }

  //pthread_attr_destroy(&attr);
  for(i = 0; i < threads_spawned; i++)
    {
      rc = pthread_join(threads[i], &status);
      if(rc)
	{
	  printf("pthread_join() error: %d\n",rc);
	  exit(-1);
	}
    }
  
  return GSL_SUCCESS;
}
#endif
