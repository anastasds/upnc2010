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
  double C_m, I_e, i_m, Mg_conc, f_pre;
  double g_bar_Na, g_bar_K, g_bar_L, g_bar_NMDA_Ca, g_bar_NMDA_syn, g_bar_AMPA;
  double E_Na, E_K, E_L, E_NMDA_Ca, E_NMDA_syn, E_AMPA_syn;
  double m_NMDA_Ca, m_NMDA_syn, s_NMDA, s_NMDA_rise, s_NMDA_fast, s_NMDA_slow;
  double s_AMPA, s_AMPA_rise, s_AMPA_fast, s_AMPA_slow;
  double Phi_NMDA, tau_NMDA_rise, tau_NMDA_fast, tau_NMDA_slow;
  double Phi_AMPA, tau_AMPA_rise, tau_AMPA_fast, tau_AMPA_slow;

  double alpha_n, beta_n, alpha_m, beta_m, alpha_h, beta_h;
  double tau_n, tau_m, tau_h, n_inf, m_inf, h_inf;
  struct neuron_params * network_params;
  long i, j, num_state_params, limit = start + num;

  if(network->size < limit)
    limit = network->size;

  // get general network parameters that are the same for each neuron
  network_params = network->neurons[0]->params;
  C_m = network_params->values[0];
  g_bar_Na = network_params->values[1];
  g_bar_K = network_params->values[2];
  g_bar_L = network_params->values[3];
  g_bar_NMDA_Ca = network_params->values[4];
  g_bar_NMDA_syn = network_params->values[5];
  g_bar_AMPA = network_params->values[6];
  E_Na = network_params->values[7];
  E_K = network_params->values[8];
  E_L = network_params->values[9];
  E_NMDA_Ca = network_params->values[10];
  E_NMDA_syn = network_params->values[11];
  E_AMPA_syn = network_params->values[12];
  Phi_NMDA = network_params->values[13];
  tau_NMDA_rise = network_params->values[14];
  tau_NMDA_fast = network_params->values[15];
  tau_NMDA_slow = network_params->values[16];
  Phi_AMPA = network_params->values[17];
  tau_AMPA_rise = network_params->values[18];
  tau_AMPA_fast = network_params->values[19];
  tau_AMPA_slow = network_params->values[20];
  Mg_conc = network_params->values[21];
  I_e = network_params->values[22];

  num_state_params = network->neurons[start]->state->num_params;

  for(i = start; i < limit; i++)
    {
      //update params in neuron_state
      for(j = 0; j < network->neurons[i]->state->num_params; j++)
	network->neurons[i]->state->values[j] = y[num_state_params*i + j];

      // find dn/dt, dm/dt, dh/dt usig some intermediate values
      // from theroetical neuroscience, by dayan & abbott

      // potassium 
      alpha_n = 0.01*(y[num_state_params*i] + 55)/(1 - exp(-0.1*(y[num_state_params*i]+55.0)));
      beta_n = 0.125*exp(-0.0125*(y[num_state_params*i]+65.0));
      tau_n = 1.0/(alpha_n + beta_n);
      n_inf = tau_n * alpha_n;
      
      //sodium
      alpha_m = 0.1*(y[num_state_params*i] + 40.0) / (1.0 - exp(-0.1*(y[num_state_params*i]+40.0)));
      beta_m = 4.0*exp(-0.0556*(y[num_state_params*i]+65.0));
      tau_m = 1.0/(alpha_m + beta_m);
      m_inf = tau_m * alpha_m;
      
      alpha_h = 0.07*exp(-0.05*(y[num_state_params*i]+65.0));
      beta_h = 1.0 / (1.0 + exp(-0.1*(y[num_state_params*i]+35.0)));
      tau_h = 1.0/(alpha_h + beta_h);
      h_inf = tau_h * alpha_h;

      // NMDAR
      m_NMDA_Ca =  1.0 / (1.0 + 0.3 * Mg_conc * exp(-0.124 * y[num_state_params*i]));
      m_NMDA_syn =  1.0 / (1.0 + 0.3 * Mg_conc * exp(-0.062 * y[num_state_params*i]));

      s_NMDA_rise = y[num_state_params * i + 4];
      s_NMDA_fast = y[num_state_params * i + 5];
      s_NMDA_slow = y[num_state_params * i + 6];
      s_NMDA = s_NMDA_rise + s_NMDA_fast + s_NMDA_slow;

      // AMPAR
      s_AMPA_rise = y[num_state_params * i + 7];
      s_AMPA_fast = y[num_state_params * i + 8];
      s_AMPA_slow = y[num_state_params * i + 9];
      s_AMPA = s_AMPA_rise + s_AMPA_fast + s_AMPA_slow;

      // figure out membrane current, i_m
      i_m = g_bar_L * (y[num_state_params*i] - E_L)
	  + g_bar_K * pow(y[num_state_params*i + 1],4.0) * (y[num_state_params*i] - E_K)
	  + g_bar_Na * pow(y[num_state_params*i + 2],3) * y[num_state_params*i + 3] * (y[num_state_params*i] - E_Na)
	  + g_bar_NMDA_Ca * s_NMDA * m_NMDA_Ca * (y[num_state_params*i] - E_NMDA_Ca)
	  + g_bar_NMDA_syn * s_NMDA * m_NMDA_syn * (y[num_state_params*i] - E_NMDA_syn)
	  + g_bar_AMPA * s_AMPA * (y[num_state_params*i] - E_AMPA_syn);
      
      f_pre = 0.0;

      // update derivatives
      f[num_state_params*i] = I_e + -1.0*i_m/C_m;
      f[num_state_params*i + 1] = (n_inf - y[num_state_params*i + 1])/(tau_n);
      f[num_state_params*i + 2] = (m_inf - y[num_state_params*i + 2])/(tau_m);
      f[num_state_params*i + 3] = (h_inf - y[num_state_params*i + 3])/(tau_h);
      f[num_state_params*i + 4] = -1.0 * Phi_NMDA * (1.0 - s_NMDA_fast - s_NMDA_slow) * f_pre - s_NMDA_rise / tau_NMDA_rise;
      f[num_state_params*i + 5] = Phi_NMDA * (0.527 - s_NMDA_fast) * f_pre - s_NMDA_fast / tau_NMDA_fast;
      f[num_state_params*i + 6] = Phi_NMDA * (0.472 - s_NMDA_slow) * f_pre - s_NMDA_slow / tau_NMDA_slow;
      f[num_state_params*i + 7] = -1.0 * Phi_AMPA * (1.0 - s_AMPA_fast - s_AMPA_slow) * f_pre - s_AMPA_rise / tau_AMPA_rise;
      f[num_state_params*i + 8] = Phi_AMPA * (0.903 - s_AMPA_fast) * f_pre - s_AMPA_fast / tau_AMPA_fast;
      f[num_state_params*i + 9] = Phi_AMPA * (0.097 - s_AMPA_slow) * f_pre - s_AMPA_slow / tau_AMPA_slow;
    }
  #ifdef THREADED
    pthread_exit(NULL);
  #endif
}

int ode_run(struct network * network, double t, double t1, double step_size, double error)
{
  const gsl_odeiv_step_type * T = gsl_odeiv_step_rk8pd;
  long i, j, num_state_params = network->neurons[0]->state->num_params;;
  long dimension = network->size * num_state_params;
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
  double y[num_state_params * network->size];
  for(i = 0; i < network->size; i++)
      for(j = 0; j < network->neurons[i]->state->num_params; j++)
	y[num_state_params*i + j] = network->neurons[i]->state->values[j];

  while (t < t1)
    {
      status = gsl_odeiv_evolve_apply(e, c, s, &sys, &t, t1, &step_size, y);
      if(status != GSL_SUCCESS)
	break;
      printf("%f %f\n",t,y[0]);
    }
  
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
