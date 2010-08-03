#include "definitions.h"
#include "includes.h"
#include "math_includes.h"
#include "neuron.h"
#include "currents.h"

double Na_current(struct network * network, long num_neuron, long num_compartment, double * f, const double * y, double t)
{
  double alpha_m, beta_m, tau_m, m_inf;
  double alpha_h, beta_h, tau_h, h_inf;
  double g_bar_Na, E_Na;
  long offset = num_neuron * network->compartments * network->neurons[num_neuron]->compartments[num_compartment]->state->num_params + num_compartment * network->neurons[num_neuron]->compartments[num_compartment]->state->num_params;
  
  g_bar_Na = y[offset + 2];
  E_Na = network->neurons[0]->params->values[1];

  // sodium, from Dayan & Abbott
  alpha_m = 0.38*(y[offset] + 29.7) / (1.0 - exp(-0.1*(y[offset]+29.7)));
  beta_m = 15.2*exp(-0.0556*(y[offset]+54.7));
  tau_m = 1.0/(alpha_m + beta_m);
  m_inf = tau_m * alpha_m;
  
  alpha_h = 0.266*exp(-0.05*(y[offset]+48.0));
  beta_h = 3.8 / (1.0 + exp(-0.1*(y[offset]+18.0)));
  tau_h = 1.0/(alpha_h + beta_h);
  h_inf = tau_h * alpha_h;

  f[offset + 11] = (m_inf - y[offset + 11])/(tau_m);
  f[offset + 12] = (h_inf - y[offset + 12])/(tau_h);
  
  return -1.0*(g_bar_Na * pow(y[offset + 11],3.0) * y[offset + 12] * (y[offset] - E_Na));
}

double Kdr_current(struct network * network, long num_neuron, long num_compartment, double * f, const double * y, double t)
{
  double alpha_n, beta_n, tau_n, n_inf;
  double g_bar_K, E_K;
  long offset = num_neuron * network->compartments * network->neurons[num_neuron]->compartments[num_compartment]->state->num_params + num_compartment * network->neurons[num_neuron]->compartments[num_compartment]->state->num_params;

  g_bar_K = y[offset + 4];
  E_K = network->neurons[0]->params->values[2];

  // potassium, from Dayan & Abbott
  alpha_n = 0.02*(y[offset] + 45.7)/(1.0 - exp(-0.1*(y[offset]+45.7)));
  beta_n = 0.25*exp(-0.0125*(y[offset]+55.7));
  tau_n = 1.0/(alpha_n + beta_n);
  n_inf = tau_n * alpha_n;
  
  f[offset + 10] = (n_inf - y[offset + 10])/(tau_n);

  return -1.0*(g_bar_K * pow(y[offset + 10],4.0) * (y[offset] - E_K));
}

double A_current(struct network * network, long num_neuron, long num_compartment, double * f, const double * y, double t)
{
  double a_inf, tau_a, b_inf, tau_b;
  double g_bar_A, E_A;
  long offset = num_neuron * network->compartments * network->neurons[num_neuron]->compartments[num_compartment]->state->num_params + num_compartment * network->neurons[num_neuron]->compartments[num_compartment]->state->num_params;

  g_bar_A = y[offset + 3];
  E_A = network->neurons[0]->params->values[3];

  // A-current, from Dayan & Abbott
  a_inf = pow((0.0761 * exp(0.0314 * (y[offset] + 94.22))) / (1.0 + exp(0.0346 * (y[offset] + 1.17))), (1.0/3.0));
  tau_a = 0.3632 + 1.158 / (1.0 + exp(0.0497 * (y[offset] + 55.96)));
  
  b_inf = pow((1.0 / (1.0 + exp(0.0688 * (y[offset] + 53.3)))),4.0);
  tau_b = 1.24 + 2.678 / (1.0 + exp(0.0624 * (y[offset] + 50.0)));

  f[offset + 13] = (a_inf - y[offset + 13])/(tau_a);
  f[offset + 14] = (b_inf - y[offset + 14])/(tau_b);
  
  return -1.0*(g_bar_A * pow(y[offset + 13],3.0) * y[offset + 14] * (y[offset] - E_A));
}

double KCa_current(struct network * network, long num_neuron, long num_compartment, double * f, const double * y, double t)
{
  double c_inf, tau_c;
  double g_bar_KCa, E_K;
  long offset = num_neuron * network->compartments * network->neurons[num_neuron]->compartments[num_compartment]->state->num_params + num_compartment * network->neurons[num_neuron]->compartments[num_compartment]->state->num_params;

  g_bar_KCa = y[offset + 5];
  E_K = network->neurons[0]->params->values[2];

  // calcium-dependent potassium current, from Dayan & Abbott
  c_inf = y[offset + 18] / (y[offset + 18] + 3.0) * (1.0 / (1.0 + exp((y[offset] + 28.3) / -12.6)));
  tau_c = 90.3 - 75.1 / (1.0 + exp((y[offset] + 46.0) / -22.7));
  
  f[offset + 15] = (c_inf - y[offset + 15])/(tau_c);

  return -1.0*(g_bar_KCa * pow(y[offset + 15],4.0) * (y[offset] - E_K));
}

double CaT_current(struct network * network, long num_neuron, long num_compartment, double * f, const double * y, double t)
{
  double M_inf, tau_M, H_inf, tau_H;
  double g_bar_CaT, E_Ca;
  long offset = num_neuron * network->compartments * network->neurons[num_neuron]->compartments[num_compartment]->state->num_params + num_compartment * network->neurons[num_neuron]->compartments[num_compartment]->state->num_params;

  g_bar_CaT = y[offset + 6];
  E_Ca = network->neurons[0]->params->values[5];

  // transient calcium current as mentioned in Porazi et al. 2003 supplement but as defined in Dayan & Abbot
  M_inf = 1.0 / (1.0 + exp((y[offset] + 57.0) / -6.2));
  tau_M = 0.612 + 1.0 / (exp((y[offset] + 132.0) / -16.7) + exp((y[offset] + 16.8) / 18.2));
  
  H_inf = 1.0 / (1.0 + exp((y[offset] + 81.0) / 4.0));
  if(y[offset] < -80.0)
    tau_H = exp((y[offset] + 467.0) / 66.6);
  else
    tau_H = 28.0 + exp((y[offset] + 22.0) / -10.5);
  
  f[offset + 16] = (M_inf - y[offset + 16])/(tau_M);
  f[offset + 17] = (H_inf - y[offset + 17])/(tau_H);

  return -1.0*(g_bar_CaT * pow(y[offset + 16], 2.0) * y[offset + 17] * (y[offset] - E_Ca));
}

double L_current(struct network * network, long num_neuron, long num_compartment, double * f, const double * y, double t)
{
  double g_bar_L, E_L;
  long offset = num_neuron * network->compartments * network->neurons[num_neuron]->compartments[num_compartment]->state->num_params + num_compartment * network->neurons[num_neuron]->compartments[num_compartment]->state->num_params;

  g_bar_L = y[offset + 1];
  E_L = network->neurons[0]->params->values[4];

  // leak current
  return -1.0*(g_bar_L * (y[offset] - E_L));
}

double NMDA_current(struct network * network, long num_neuron, long num_compartment, double * f, const double * y, double t, double i_CaT)
{
  double m_NMDA_Ca, m_NMDA_syn, s_NMDA_rise, s_NMDA_fast, s_NMDA_slow, s_NMDA;
  double Phi, beta_soma, Ca_resting_conc_soma, Ca_diffusion_rate, eta_buff;
  double i_Ca, i_Ca_NMDA, beta_spine, Ca_resting_conc_spine, beta_buff;
  double Mg_conc, g_bar_NMDA_Ca, g_bar_NMDA_syn, E_NMDA_syn, E_Ca;
  double Phi_NMDA, tau_NMDA_rise, tau_NMDA_fast, tau_NMDA_slow, tau_Ca;
  double f_pre;

  long num_state_params = network->neurons[num_neuron]->compartments[num_compartment]->state->num_params;

  long offset = num_neuron * network->compartments * num_state_params + num_compartment * num_state_params;

  g_bar_NMDA_Ca = y[offset + 7];
  g_bar_NMDA_syn = y[offset + 8];
  E_Ca = network->neurons[0]->params->values[5];
  E_NMDA_syn = network->neurons[0]->params->values[6];
  Mg_conc = network->neurons[0]->params->values[17];
  f_pre = 0.0;

  // NMDAR, from Rubin et al. 2005
  m_NMDA_Ca =  1.0 / (1.0 + 0.3 * Mg_conc * exp(-0.124 * y[offset]));
  m_NMDA_syn =  1.0 / (1.0 + 0.3 * Mg_conc * exp(-0.062 * y[offset]));

  s_NMDA_rise = y[offset + 19];
  s_NMDA_fast = y[offset + 20];
  s_NMDA_slow = y[offset + 21];
  s_NMDA = s_NMDA_rise + s_NMDA_fast + s_NMDA_slow;

  Phi_NMDA = network->neurons[0]->params->values[8];
  tau_NMDA_rise = network->neurons[0]->params->values[9];
  tau_NMDA_fast = network->neurons[0]->params->values[10];
  tau_NMDA_slow = network->neurons[0]->params->values[11];
  tau_Ca = network->neurons[0]->params->values[16];
  Ca_resting_conc_soma = network->neurons[0]->params->values[18];
  Ca_resting_conc_spine = network->neurons[0]->params->values[19];
  Phi = network->neurons[0]->params->values[20];
  beta_soma = network->neurons[0]->params->values[21];
  beta_spine = network->neurons[0]->params->values[22];
  eta_buff = network->neurons[0]->params->values[23];
  beta_buff = network->neurons[0]->params->values[24];
  Ca_diffusion_rate = network->neurons[0]->params->values[25];

  // i_Ca
  i_Ca = i_CaT; // i_CaT passed as argument
  i_Ca_NMDA = -1.0 * (g_bar_NMDA_Ca * s_NMDA * m_NMDA_Ca * (y[offset] - E_Ca));

  // compartment 0 is spine, 1 is soma
  if(num_compartment == 0)
    {
      f[offset + 18] = Phi * (i_Ca + i_Ca_NMDA) - beta_spine * (y[offset + 18] - Ca_resting_conc_spine) - beta_spine / eta_buff * pow(y[offset + 18],2.0) - beta_buff * y[offset + 18];
    }
  else
    {
      long num_state_params = network->neurons[num_neuron]->compartments[0]->state->num_params;
      f[offset + 18] = Phi * i_Ca - beta_soma * (y[offset + 18] - Ca_resting_conc_soma) + (y[offset + 18 - num_state_params] - y[offset + 18]) / Ca_diffusion_rate - beta_soma / eta_buff * pow(y[offset + 18], 2.0);
    }

  f[offset + 19] = -1.0 * Phi_NMDA * (1.0 - s_NMDA_fast - s_NMDA_slow) * f_pre - s_NMDA_rise / tau_NMDA_rise;

  f[offset + 20] = Phi_NMDA * (0.527 - s_NMDA_fast) * f_pre - s_NMDA_fast / tau_NMDA_fast;
  f[offset + 21] = Phi_NMDA * (0.473 - s_NMDA_slow) * f_pre - s_NMDA_slow / tau_NMDA_slow;

  return 1.0*(g_bar_NMDA_syn * s_NMDA * m_NMDA_syn * (y[offset] - E_NMDA_syn));
}

double AMPA_current(struct network * network, long num_neuron, long num_compartment, double * f, const double * y, double t)
{
  double s_AMPA_rise, s_AMPA_fast, s_AMPA_slow, s_AMPA;
  double Phi_AMPA, tau_AMPA_rise, tau_AMPA_fast, tau_AMPA_slow;
  double g_bar_AMPA, E_AMPA_syn, f_pre;
  long offset = num_neuron * network->compartments * network->neurons[num_neuron]->compartments[num_compartment]->state->num_params + num_compartment * network->neurons[num_neuron]->compartments[num_compartment]->state->num_params;

  g_bar_AMPA = y[offset + 9];
  E_AMPA_syn = network->neurons[0]->params->values[7];
  f_pre = 0.0;

  // AMPAR, from Rubin et al. 2005
  s_AMPA_rise = y[offset + 22];
  s_AMPA_fast = y[offset + 23];
  s_AMPA_slow = y[offset + 24];
  s_AMPA = s_AMPA_rise + s_AMPA_fast + s_AMPA_slow;
  
  Phi_AMPA = network->neurons[0]->params->values[12];
  tau_AMPA_rise = network->neurons[0]->params->values[13];
  tau_AMPA_fast = network->neurons[0]->params->values[14];
  tau_AMPA_slow = network->neurons[0]->params->values[15];

  f[offset + 22] = -1.0 * Phi_AMPA * (1.0 - s_AMPA_fast - s_AMPA_slow) * f_pre - s_AMPA_rise / tau_AMPA_rise;
  f[offset + 23] = Phi_AMPA * (0.903 - s_AMPA_fast) * f_pre - s_AMPA_fast / tau_AMPA_fast;
  f[offset + 24] = Phi_AMPA * (0.097 - s_AMPA_slow) * f_pre - s_AMPA_slow / tau_AMPA_slow;

  return 1.0*(g_bar_AMPA * s_AMPA * (y[offset] - E_AMPA_syn));
}
