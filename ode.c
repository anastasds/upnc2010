#include "definitions.h"
#include "includes.h"
#include "math_includes.h"
#include "neuron.h"
#include "ode.h"

#ifdef THREADED
void * ode_update_neurons_threaded(void * thread_params)
{
    struct thread_params * params = (struct thread_params *)thread_params;
    ode_update_neurons(params->network, params->start, params->num, params->y, params->f);
    pthread_exit(NULL);
}
#endif

void ode_update_neurons(struct network * network, long start, long num, const double * y, double * f)
{
  double C_m, I_e, i_m, Mg_conc, Ca_resting_conc_soma, Ca_resting_conc_spine, f_pre, tau_Ca, i_Ca_soma, i_Ca_spine;
  double g_bar_Na, g_bar_K, g_bar_A, g_bar_KCa, g_bar_L, g_bar_NMDA_Ca, g_bar_NMDA_syn, g_bar_AMPA;
  double E_Na, E_K, E_A, E_L, E_NMDA_Ca, E_NMDA_syn, E_AMPA_syn;
  double m_NMDA_Ca, m_NMDA_syn, s_NMDA, s_NMDA_rise, s_NMDA_fast, s_NMDA_slow;
  double s_AMPA, s_AMPA_rise, s_AMPA_fast, s_AMPA_slow;
  double Phi_NMDA, tau_NMDA_rise, tau_NMDA_fast, tau_NMDA_slow;
  double Phi_AMPA, tau_AMPA_rise, tau_AMPA_fast, tau_AMPA_slow;
  double Phi_Ca, beta_soma, beta_spine, beta_buff, Ca_diffusion_rate, eta_buff;

  double alpha_n, beta_n, alpha_m, beta_m, alpha_h, beta_h;
  double tau_n, tau_m, tau_h, tau_a, tau_b, tau_c, n_inf, m_inf, h_inf, a_inf, b_inf, c_inf;
  long i, j, offset, num_state_params, limit = start + num;
  struct neuron_params * network_params;

  if(network->size < limit)
    limit = network->size;

  // get general network parameters that are the same for each neuron
  network_params = network->neurons[start]->params;
  C_m = network_params->values[0];
  g_bar_Na = network_params->values[1];
  g_bar_K = network_params->values[2];
  g_bar_A = network_params->values[3];
  g_bar_KCa = network_params->values[4];
  g_bar_L = network_params->values[5];
  g_bar_NMDA_Ca = network_params->values[6];
  g_bar_NMDA_syn = network_params->values[7];
  g_bar_AMPA = network_params->values[8];
  E_Na = network_params->values[9];
  E_K = network_params->values[10];
  E_A = network_params->values[11];
  E_L = network_params->values[12];
  E_NMDA_Ca = network_params->values[13];
  E_NMDA_syn = network_params->values[14];
  E_AMPA_syn = network_params->values[15];
  Phi_NMDA = network_params->values[16];
  tau_NMDA_rise = network_params->values[17];
  tau_NMDA_fast = network_params->values[18];
  tau_NMDA_slow = network_params->values[19];
  Phi_AMPA = network_params->values[20];
  tau_AMPA_rise = network_params->values[21];
  tau_AMPA_fast = network_params->values[22];
  tau_AMPA_slow = network_params->values[23];
  tau_Ca = network_params->values[24];
  Mg_conc = network_params->values[25];
  Ca_resting_conc_soma = network_params->values[26];
  Ca_resting_conc_spine = network_params->values[27];
  Phi_Ca = network_params->values[28];
  beta_soma = network_params->values[29];
  beta_spine = network_params->values[30];
  eta_buff = network_params->values[31];
  beta_buff = network_params->values[32];
  Ca_diffusion_rate = network_params->values[33];
  I_e = network_params->values[34];

  num_state_params = network->neurons[start]->state->num_params;

  for(i = start; i < limit; i++)
    {
      offset = num_state_params * i;

      // update params in neuron_state
      for(j = 0; j < network->neurons[i]->state->num_params; j++)
	network->neurons[i]->state->values[j] = y[num_state_params*i + j];

      // potassium, from Dayan & Abbott
      alpha_n = 0.02*(y[offset] + 45.7)/(1.0 - exp(-0.1*(y[offset]+45.7)));
      beta_n = 0.25*exp(-0.0125*(y[offset]+55.7));
      tau_n = 1.0/(alpha_n + beta_n);
      n_inf = tau_n * alpha_n;
      
      // A-current, from Dayan & Abbott
      a_inf = pow((0.0761 * exp(0.0314 * (y[offset] + 94.22))) / (1.0 + exp(0.0346 * (y[offset] + 1.17))), (1.0/3.0));
      tau_a = 0.3632 + 1.158 / (1.0 + exp(0.0497 * (y[offset] + 55.96)));
      
      b_inf = pow((1.0 / (1.0 + exp(0.0688 * (y[offset] + 53.3)))),4.0);
      tau_b = 1.24 + 2.678 / (1.0 + exp(0.0624 * (y[offset] + 50.0)));

      // sodium, from Dayan & Abbott
      alpha_m = 0.38*(y[offset] + 29.7) / (1.0 - exp(-0.1*(y[offset]+29.7)));
      beta_m = 15.2*exp(-0.0556*(y[offset]+54.7));
      tau_m = 1.0/(alpha_m + beta_m);
      m_inf = tau_m * alpha_m;
      
      alpha_h = 0.266*exp(-0.05*(y[offset]+48.0));
      beta_h = 3.8 / (1.0 + exp(-0.1*(y[offset]+18.0)));
      tau_h = 1.0/(alpha_h + beta_h);
      h_inf = tau_h * alpha_h;

      // calcium-dependent potassium current, from Dayan & Abbott
      c_inf = y[offset + 7] / (y[offset + 7] + 3.0) * (1.0 / (1.0 + exp((y[offset] + 28.3) / -12.6)));
      tau_c = 90.3 - 75.1 / (1.0 + exp((y[offset] + 46.0) / -22.7));

      // NMDAR, from Rubin et al. 2005
      m_NMDA_Ca =  1.0 / (1.0 + 0.3 * Mg_conc * exp(-0.124 * y[offset]));
      m_NMDA_syn =  1.0 / (1.0 + 0.3 * Mg_conc * exp(-0.062 * y[offset]));

      s_NMDA_rise = y[offset + 6];
      s_NMDA_fast = y[offset + 7];
      s_NMDA_slow = y[offset + 8];
      s_NMDA = s_NMDA_rise + s_NMDA_fast + s_NMDA_slow;

      // AMPAR, from Rubin et al. 2005
      s_AMPA_rise = y[offset + 9];
      s_AMPA_fast = y[offset + 10];
      s_AMPA_slow = y[offset + 11];
      s_AMPA = s_AMPA_rise + s_AMPA_fast + s_AMPA_slow;

      // membrane current
      i_m = g_bar_L * (y[offset] - E_L)
	  + g_bar_K * pow(y[offset + 1],4.0) * (y[offset] - E_K)
	  + g_bar_A * pow(y[offset + 4],3) * y[offset + 5] * (y[offset] - E_A)
	  + g_bar_Na * pow(y[offset + 2],3) * y[offset + 3] * (y[offset] - E_Na)
	  + g_bar_KCa * pow(y[offset + 6],4) * (y[offset] - E_K)
	  + g_bar_NMDA_Ca * s_NMDA * m_NMDA_Ca * (y[offset] - E_NMDA_Ca)
	  + g_bar_NMDA_syn * s_NMDA * m_NMDA_syn * (y[offset] - E_NMDA_syn)
	  + g_bar_AMPA * s_AMPA * (y[offset] - E_AMPA_syn);
      
      f_pre = 2.0;

      // i_Ca here is i_CaL from Poirazi et al. 2003
      i_Ca_soma = 0.0;
      i_Ca_spine = 0.0;

      // update derivatives
      f[offset] = I_e + -1.0*i_m/C_m;
      f[offset + 1] = (n_inf - y[offset + 1])/(tau_n);
      f[offset + 2] = (m_inf - y[offset + 2])/(tau_m);
      f[offset + 3] = (h_inf - y[offset + 3])/(tau_h);
      f[offset + 4] = (a_inf - y[offset + 4])/(tau_a);
      f[offset + 5] = (b_inf - y[offset + 5])/(tau_b);
      f[offset + 6] = (c_inf - y[offset + 6])/(tau_c);
      f[offset + 7] = Phi_Ca * i_Ca_soma - beta_soma * (y[offset + 7] - Ca_resting_conc_soma) + (y[offset + 8] - y[offset + 7]) / Ca_diffusion_rate - beta_soma / eta_buff * pow(y[offset + 7], 2.0);
      f[offset + 8] = Phi_Ca * (i_Ca_spine + i_Ca_soma) - beta_spine * (y[offset + 8] - Ca_resting_conc_spine) - beta_spine / eta_buff * pow(y[offset + 8],2.0) - beta_buff * y[offset + 8];
      f[offset + 9] = -1.0 * Phi_NMDA * (1.0 - s_NMDA_fast - s_NMDA_slow) * f_pre - s_NMDA_rise / tau_NMDA_rise;
      f[offset + 10] = Phi_NMDA * (0.527 - s_NMDA_fast) * f_pre - s_NMDA_fast / tau_NMDA_fast;
      f[offset + 11] = Phi_NMDA * (0.472 - s_NMDA_slow) * f_pre - s_NMDA_slow / tau_NMDA_slow;
      f[offset + 12] = -1.0 * Phi_AMPA * (1.0 - s_AMPA_fast - s_AMPA_slow) * f_pre - s_AMPA_rise / tau_AMPA_rise;
      f[offset + 13] = Phi_AMPA * (0.903 - s_AMPA_fast) * f_pre - s_AMPA_fast / tau_AMPA_fast;
      f[offset + 14] = Phi_AMPA * (0.097 - s_AMPA_slow) * f_pre - s_AMPA_slow / tau_AMPA_slow;
    }
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
