#include "definitions.h"
#include "includes.h"
#include "math_includes.h"

int hh_ode (double t, const double y[], double f[], void *params)
{
  
  double V_T, C_m, g_bar_Na, g_bar_K, g_bar_L, E_Na, E_K, E_L, A_alpha, B_alpha, A_beta, B_beta;
  double alpha_n, beta_n, alpha_m, beta_m, alpha_h, beta_h;
  double i_m, A, I_e;
  double tau_n, tau_m, tau_h, n_inf, m_inf, h_inf;
  struct neuron_ode_params * state = (struct neuron_ode_params *)params;
  // f[0] - V', f[1] - n', f[2] - m', f[3] - h'
  
  V_T = state->params->values[0];
  C_m = state->params->values[1];
  g_bar_Na = state->params->values[2];
  g_bar_K = state->params->values[3];
  g_bar_L = state->params->values[4];
  E_Na = state->params->values[5];
  E_K = state->params->values[6];
  E_L = state->params->values[7];
  A_alpha = state->params->values[8];
  B_alpha = state->params->values[9];
  A_beta = state->params->values[10];
  B_beta = state->params->values[11];
  A = state->params->values[12];
  I_e = state->params->values[13];

  //alpha_n = A_alpha*exp(-1.0*B_alpha*y[0]/V_T);
  //beta_n = A_beta*exp(-1.0*B_beta*y[0]/V_T);
  alpha_n = 0.01*(y[0] + 55)/(1 - exp(-0.1*(y[0]+55.0)));
  beta_n = 0.125*exp(-0.0125*(y[0]+65.0));
  n_inf = alpha_n / (alpha_n + beta_n);
  tau_n = 1.0/(alpha_n + beta_n);

  alpha_m = 0.1*(y[0] + 40.0) / (1.0 - exp(-0.1*(y[0]+40.0)));
  beta_m = 4.0*exp(-0.0556*(y[0]+65.0));
  m_inf = alpha_m / (alpha_m + beta_m);
  tau_m = 1.0/(alpha_m + beta_m);
  
  alpha_h = 0.07*exp(-0.05*(y[0]+65.0));
  beta_h = 1.0 / (1.0 + exp(-0.1*(y[0]+35.0)));
  h_inf = alpha_h / (alpha_h + beta_h);
  tau_h = 1.0/(alpha_h + beta_h);

  i_m = g_bar_L*(y[0] - E_L) + g_bar_K*pow(y[1],4.0)*(y[0] - E_K) + g_bar_Na*pow(y[2],3)*y[3]*(y[0] - E_Na);

  //f[0] = (1.0/c_m) * (-1.0*i_m + I_e/A);
  f[0] = I_e + -1.0*i_m/C_m;
  f[1] = (n_inf - y[1])/(tau_n);
  f[2] = (m_inf - y[2])/(tau_m);
  f[3] = (h_inf - y[3])/(tau_h);
  /*
  f[1] = alpha_n*(1.0-y[1]) - beta_n*y[1];
  f[2] = alpha_m*(1-y[2]) - beta_m*y[2];
  f[3] = alpha_h*(1.0-y[3]) - beta_h*y[3];
  */
  //printf("%.5f %.5f %.5f %.5f\n",f[0],f[1],f[2],f[3]);

  return GSL_SUCCESS;
}

int hh_ode_jac (double t, const double y[], double *dfdy, double dfdt[], void *params)
{
  return GSL_SUCCESS;
}
    
int ode_test(struct neuron_state * neuron_state, struct neuron_params * neuron_params)
{
  int i = 0;
  const gsl_odeiv_step_type * T = gsl_odeiv_step_rk8pd;

  gsl_odeiv_step * s = gsl_odeiv_step_alloc(T, 4);
  gsl_odeiv_control * c = gsl_odeiv_control_y_new(1e-6, 0.0);
  gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc(4);

  struct neuron_ode_params * params = (struct neuron_ode_params *)malloc(sizeof(struct neuron_ode_params));
  params->state = neuron_state;
  params->params = neuron_params;

  // params: function, jacobian, dimension, void * params
  gsl_odeiv_system sys = {hh_ode, NULL, 4, params};
  
  double t = 0.0, t1 = 100.0;
  double h = 1e-6;
  double y[4] = { neuron_state->values[0], neuron_state->values[1], neuron_state->values[2], neuron_state->values[3]};
  
  while (t < t1)
    {
      int status = gsl_odeiv_evolve_apply (e, c, s,
					   &sys, 
					   &t, t1,
					   &h, y);
      
      if (status != GSL_SUCCESS)
	break;
      if(i++ % 1 == 0)
	printf("%.5e, %.5e, %.5e, %.5e, %.5e\n", t, y[0], y[1], y[2], y[3]);
    }
  
  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free (c);
  gsl_odeiv_step_free (s);
  free(params);
  return 0;
}
