#include "definitions.h"
#include "includes.h"
#include "math_includes.h"
#include "neuron.h"
#include "plasticity.h"

double evolve_P(struct network * network, const double * y, long offset)
{
  struct neuron_params * params = network->neurons[0]->params;
  double c_p = params->values[39];
  double tau_p = params->values[34];

  double chi_spine = y[offset + 18];
  double A = y[offset + 28];
  double P = y[offset + 26];

  return (p_sigma(network,chi_spine) - c_p*A*P)/tau_p;
}

double evolve_V(struct network * network, const double * y, long offset)
{
  struct neuron_params * params = network->neurons[0]->params;
  double tau_v = params->values[36];

  double chi_spine = y[offset + 18];
  double V = y[offset + 27];

  return (v_sigma(network,chi_spine) - V)/tau_v;
}

double evolve_A(struct network * network, const double * y, long offset)
{
  struct neuron_params * params = network->neurons[0]->params;
  double tau_a = params->values[35];

  double chi_spine = y[offset + 18];
  double A = y[offset + 28];

  return (a_sigma(network,chi_spine) - A)/tau_a;
}

double evolve_B(struct network * network, const double * y, long offset)
{
  struct neuron_params * params = network->neurons[0]->params;
  double c_d = params->values[53];
  double tau_b = params->values[38];

  double A = y[offset + 28];
  double B = y[offset + 29];
  double V = y[offset + 27];

  return (b_sigma(network,A) - B - c_d*B*V)/tau_b;
}

double evolve_D(struct network * network, const double * y, long offset)
{
  struct neuron_params * params = network->neurons[0]->params;
  double tau_d = params->values[37];

  double B = y[offset + 29];
  double D = y[offset + 30];

  return (d_sigma(network,B) - D)/tau_d;
}

double evolve_W(struct network * network, const double * y, long offset)
{
  struct neuron_params * params = network->neurons[0]->params;
  double alpha_w = params->values[41];
  double p = params->values[54];
  double k_p = params->values[48];
  double beta_w = params->values[42];
  double d = params->values[55];
  double k_d = params->values[49];
  double tau_w = params->values[40];

  double P = y[offset + 26];
  double D = y[offset + 30];
  double W = y[offset + 31];

  return (alpha_w/(1.0 + exp((P - p)/k_p)) - beta_w/(1.0 + exp((D - d)/k_d)) - W) / tau_w;
}

double p_sigma(struct network * network, double x)
{
  struct neuron_params * params = network->neurons[0]->params;
  double pHC = params->values[29];
  double pHN = params->values[43];

  return (10.0 * pow(x/pHC,pHN)) / (1.0 + pow(x/pHC,pHN));
}

double a_sigma(struct network * network, double x)
{
  struct neuron_params * params = network->neurons[0]->params;
  double aHC = params->values[30];
  double aHN = params->values[44];

  return pow(x/aHC,aHN) / (1.0 + pow(x/aHC,aHN));
}

double v_sigma(struct network * network, double x)
{
  struct neuron_params * params = network->neurons[0]->params;
  double alpha_v = params->values[50];
  double theta_v = params->values[31];
  double sigma_v = params->values[45];

  return alpha_v / (1.0 + exp((x - theta_v)/sigma_v));
}

double d_sigma(struct network * network, double x)
{
  struct neuron_params * params = network->neurons[0]->params;
  double alpha_d = params->values[51];
  double theta_d = params->values[32];
  double sigma_d = params->values[46];

  return alpha_d / (1.0 + exp((x - theta_d)/sigma_d));
}

double b_sigma(struct network * network, double x)
{
  struct neuron_params * params = network->neurons[0]->params;
  double alpha_b = params->values[52];
  double theta_b = params->values[33];
  double sigma_b = params->values[47];

  return alpha_b / (1.0 + exp((x - theta_b)/sigma_b));
}
