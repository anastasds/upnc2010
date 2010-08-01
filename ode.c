#include "definitions.h"
#include "includes.h"
#include "math_includes.h"
#include "neuron.h"
#include "ode.h"

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
  double C_m, I_e, i_m, Mg_conc, Ca_resting_conc_soma, Ca_resting_conc_spine, f_pre, tau_Ca, i_Ca_soma, i_Ca_spine, i_Ca_NMDA;
  double g_bar_Na, g_bar_K, g_bar_A, g_bar_KCa, g_bar_CaT, g_bar_L, g_bar_NMDA_Ca, g_bar_NMDA_syn, g_bar_AMPA;
  double E_Na, E_K, E_A, E_L, E_NMDA_Ca, E_NMDA_syn, E_AMPA_syn;
  double m_NMDA_Ca, m_NMDA_syn, s_NMDA, s_NMDA_rise, s_NMDA_fast, s_NMDA_slow;
  double s_AMPA, s_AMPA_rise, s_AMPA_fast, s_AMPA_slow;
  double Phi_NMDA, tau_NMDA_rise, tau_NMDA_fast, tau_NMDA_slow;
  double Phi_AMPA, tau_AMPA_rise, tau_AMPA_fast, tau_AMPA_slow;
  double Phi_Ca, beta_soma, beta_spine, beta_buff, Ca_diffusion_rate, eta_buff;
  double i_L, i_Kdr, i_A, i_KCa, i_CaT, i_Na, i_NMDA, i_AMPA, i_in;

  double alpha_n, beta_n, alpha_m, beta_m, alpha_h, beta_h;
  double tau_n, tau_m, tau_h, tau_a, tau_b, tau_c, tau_M, tau_H, n_inf, m_inf, h_inf, a_inf, b_inf, c_inf, M_inf, H_inf;
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
  g_bar_CaT = network_params->values[5];
  g_bar_L = network_params->values[6];
  g_bar_NMDA_Ca = network_params->values[7];
  g_bar_NMDA_syn = network_params->values[8];
  g_bar_AMPA = network_params->values[9];
  E_Na = network_params->values[10];
  E_K = network_params->values[11];
  E_A = network_params->values[12];
  E_L = network_params->values[13];
  E_NMDA_Ca = network_params->values[14];
  E_NMDA_syn = network_params->values[15];
  E_AMPA_syn = network_params->values[16];
  Phi_NMDA = network_params->values[17];
  tau_NMDA_rise = network_params->values[18];
  tau_NMDA_fast = network_params->values[19];
  tau_NMDA_slow = network_params->values[20];
  Phi_AMPA = network_params->values[21];
  tau_AMPA_rise = network_params->values[22];
  tau_AMPA_fast = network_params->values[23];
  tau_AMPA_slow = network_params->values[24];
  tau_Ca = network_params->values[25];
  Mg_conc = network_params->values[26];
  Ca_resting_conc_soma = network_params->values[27];
  Ca_resting_conc_spine = network_params->values[28];
  Phi_Ca = network_params->values[29];
  beta_soma = network_params->values[30];
  beta_spine = network_params->values[31];
  eta_buff = network_params->values[32];
  beta_buff = network_params->values[33];
  Ca_diffusion_rate = network_params->values[34];
  I_e = network_params->values[35];

  num_state_params = network->neurons[start]->compartments[0]->state->num_params;

  for(i = start; i < limit; i++)
    {
      offset = num_state_params * network->compartments * i;

      // update params in neuron_state
      for(j = 0; j < network->neurons[i]->compartments[0]->state->num_params; j++)
	network->neurons[i]->compartments[0]->state->values[j] = y[offset + j];

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
      c_inf = y[offset + 9] / (y[offset + 9] + 3.0) * (1.0 / (1.0 + exp((y[offset] + 28.3) / -12.6)));
      tau_c = 90.3 - 75.1 / (1.0 + exp((y[offset] + 46.0) / -22.7));

      // NMDAR, from Rubin et al. 2005
      m_NMDA_Ca =  1.0 / (1.0 + 0.3 * Mg_conc * exp(-0.124 * y[offset]));
      m_NMDA_syn =  1.0 / (1.0 + 0.3 * Mg_conc * exp(-0.062 * y[offset]));

      s_NMDA_rise = y[offset + 11];
      s_NMDA_fast = y[offset + 12];
      s_NMDA_slow = y[offset + 13];
      s_NMDA = s_NMDA_rise + s_NMDA_fast + s_NMDA_slow;

      // AMPAR, from Rubin et al. 2005
      s_AMPA_rise = y[offset + 14];
      s_AMPA_fast = y[offset + 15];
      s_AMPA_slow = y[offset + 16];
      s_AMPA = s_AMPA_rise + s_AMPA_fast + s_AMPA_slow;

      // transient calcium current as mentioned in Porazi et al. 2003 supplement but as defined in Dayan & Abbot
      M_inf = 1.0 / (1.0 + exp((y[offset] + 57.0) / -6.2));
      tau_M = 0.612 + 1.0 / (exp((y[offset] + 132.0) / -16.7) + exp((y[offset] + 16.8) / 18.2));
      
      H_inf = 1.0 / (1.0 + exp((y[offset] + 81.0) / 4.0));
      if(y[offset] < -80.0)
	tau_H = exp((y[offset] + 467.0) / 66.6);
      else
	tau_H = 28.0 + exp((y[offset] + 22.0) / -10.5);
      
      
      // various currents
      i_L = -1.0*(g_bar_L * (y[offset] - E_L));
      i_Kdr = -1.0*(g_bar_K * pow(y[offset + 1],4.0) * (y[offset] - E_K));
      i_A = -1.0*(g_bar_A * pow(y[offset + 4],3.0) * y[offset + 5] * (y[offset] - E_A));
      i_Na = -1.0*(g_bar_Na * pow(y[offset + 2],3.0) * y[offset + 3] * (y[offset] - E_Na));
      i_KCa = -1.0*(g_bar_KCa * pow(y[offset + 6],4.0) * (y[offset] - E_K));
      i_CaT = -1.0*(g_bar_CaT * pow(y[offset + 7], 2.0) * y[offset + 8] * (y[offset] - E_NMDA_Ca));
      i_NMDA = -1.0*(g_bar_NMDA_syn * s_NMDA * m_NMDA_syn * (y[offset] - E_NMDA_syn));
      i_AMPA = -1.0*(g_bar_AMPA * s_AMPA * (y[offset] - E_AMPA_syn));

      // presynaptic input
      i_in = i_NMDA + i_AMPA;

      // i_Ca
      i_Ca_soma = i_CaT;
      i_Ca_spine = i_CaT;
      i_Ca_NMDA = g_bar_NMDA_Ca * s_NMDA * m_NMDA_Ca * (y[offset] - E_NMDA_Ca);

      // membrane current
      i_m = i_L + i_Kdr + i_A + i_KCa + i_CaT + i_Na + i_in;

      i_m = i_L + i_Na + i_Kdr + i_A + i_CaT + i_KCa;
            
      /*
      i_m = g_bar_L * (y[offset] - E_L)
	  + g_bar_K * pow(y[offset + 1],4.0) * (y[offset] - E_K)
	  + g_bar_A * pow(y[offset + 4],3.0) * y[offset + 5] * (y[offset] - E_A)
	  + g_bar_Na * pow(y[offset + 2],3.0) * y[offset + 3] * (y[offset] - E_Na)
	  + g_bar_KCa * pow(y[offset + 6],4.0) * (y[offset] - E_K)
	  + g_bar_NMDA_syn * s_NMDA * m_NMDA_syn * (y[offset] - E_NMDA_syn)
	  + g_bar_AMPA * s_AMPA * (y[offset] - E_AMPA_syn);
      */

      f_pre = 0.0;

      // update derivatives
      f[offset] = I_e + i_m/C_m;
      f[offset + 1] = (n_inf - y[offset + 1])/(tau_n);
      f[offset + 2] = (m_inf - y[offset + 2])/(tau_m);
      f[offset + 3] = (h_inf - y[offset + 3])/(tau_h);
      f[offset + 4] = (a_inf - y[offset + 4])/(tau_a);
      f[offset + 5] = (b_inf - y[offset + 5])/(tau_b);
      f[offset + 6] = (c_inf - y[offset + 6])/(tau_c);
      f[offset + 7] = (M_inf - y[offset + 7])/(tau_M);
      f[offset + 8] = (H_inf - y[offset + 8])/(tau_H);
      f[offset + 9] = Phi_Ca * i_Ca_soma - beta_soma * (y[offset + 9] - Ca_resting_conc_soma) + (y[offset + 10] - y[offset + 9]) / Ca_diffusion_rate - beta_soma / eta_buff * pow(y[offset + 9], 2.0);
      f[offset + 10] = Phi_Ca * (i_Ca_spine + i_Ca_NMDA) - beta_spine * (y[offset + 10] - Ca_resting_conc_spine) - beta_spine / eta_buff * pow(y[offset + 10],2.0) - beta_buff * y[offset + 10];
      f[offset + 11] = -1.0 * Phi_NMDA * (1.0 - s_NMDA_fast - s_NMDA_slow) * f_pre - s_NMDA_rise / tau_NMDA_rise;
      f[offset + 12] = Phi_NMDA * (0.527 - s_NMDA_fast) * f_pre - s_NMDA_fast / tau_NMDA_fast;
      f[offset + 13] = Phi_NMDA * (0.472 - s_NMDA_slow) * f_pre - s_NMDA_slow / tau_NMDA_slow;
      f[offset + 14] = -1.0 * Phi_AMPA * (1.0 - s_AMPA_fast - s_AMPA_slow) * f_pre - s_AMPA_rise / tau_AMPA_rise;
      f[offset + 15] = Phi_AMPA * (0.903 - s_AMPA_fast) * f_pre - s_AMPA_fast / tau_AMPA_fast;
      f[offset + 16] = Phi_AMPA * (0.097 - s_AMPA_slow) * f_pre - s_AMPA_slow / tau_AMPA_slow;

      // now that everything for this neuron has been updated, 
      // let's update all neurons on which this neuron synapses
      /*
      for(j = 0; j < network->neurons[i]->num_links; j++)
	{
	  next = network->neurons[i]->links[j]->to;
	  network->neurons[next]->
	}
      */
    }
}

int ode_run(struct network * network, double t, double t1, double step_size, double error)
{
  const gsl_odeiv_step_type * T = gsl_odeiv_step_rk8pd;
  long i, j, k, num_state_params = network->neurons[0]->compartments[0]->state->num_params;;
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

  // set up ode system, four odes per neuron
  // 0 = V, 1 = n, 2 = m, 3 = h
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
      printf("%f ", t);
      /*
      for(j = 0; j < network->size; j++)
	printf("%f ", y[num_state_params*j]);
      */
      for(j = 0; j < num_state_params; j++)
	printf("%f ",y[j]);
      printf("\n");
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
