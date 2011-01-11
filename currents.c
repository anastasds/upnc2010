#include "definitions.h"
#include "includes.h"
#include "math_includes.h"
#include "neuron.h"
#include "currents.h"
#include "stimulate.h"

double max(double a, double b)
{
  return a >= b ? a : b;
}

double min(double a, double b)
{
  return a <= b ? a : b;
}

double Na_current(struct network * network, long num_neuron, long num_compartment, double * f, const double * y, double t)
{
  double alpha_m, beta_m, tau_m, m_inf;
  double alpha_h, beta_h, tau_h, h_inf;
  double g_bar_Na, E_Na;
  long offset = calc_y_array_offset(network, num_neuron, num_compartment);
  
  g_bar_Na = y[offset + 2];
  E_Na = network->neurons[num_neuron]->params->values[1];

  /*
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
  */

  double alphams, betams, alphahs, betahs, M_Spininf, M_Spintau, H_Somainf, H_Spininf, H_Somatau, H_Spintau, I_Spininf, I_Spintau, Q, M_Soma;
  double tempC = network->neurons[num_neuron]->params->values[27];

  // Q=96480/(8.315*(273.16+tempC))
  Q = 96480.0 / (8.315 * (273.16 + tempC));

  // compartment 0 is the spine
  if(num_compartment == 0)
    {
      // M_Spininf=1/(1+exp((-vSpin-40)/3))
      M_Spininf = 1.0 / (1.0 + exp((y[offset] + 40.0)/ -3.0));

      // M_Spintau=max(0.1,0.05)
      M_Spintau = 0.1;

      // H_Spininf=1/(1+exp((vSpin+45)/3))
      H_Spininf = 1.0 / (1.0 + exp((y[offset] + 45.0) / 3.0));

      // H_Spintau=0.5
      H_Spintau = 0.5;
      
      // Na attenuation. Set to 1 since it's not clear if the spine should show attenuation.
      // par natt=0
      // I_Spininf=(1+natt*exp((vSpin+60)/2))/(1+exp((vSpin+60)/2))
      I_Spininf = 1.0 / (1.0 + exp((y[offset] + 60.0) / 2.0));

      // I_Spintau=max(0.1,(0.00333*exp(0.0024*(vSpin+60)*Q))/(1+exp(0.0012*(vSpin+60)*Q)))
      I_Spintau = max(0.1,(0.00333 * exp(0.0024 * (y[offset] + 60.0) * Q))/(1.0 + exp(0.0012 * (y[offset] + 60.0) * Q)));
            
      // M_Spin'=(M_Spininf-M_Spin)/(M_Spintau)
      f[offset + 11] = (M_Spininf - y[offset + 11]) / (M_Spintau);

      // H_Spin'=(H_Spininf-H_Spin)/(H_Spintau)
      f[offset + 12] = (H_Spininf - y[offset + 12]) / (H_Spintau);

      // I_Spin'=(I_Spininf-I_Spin)/(I_Spintau)
      f[offset + 25] = (I_Spininf - y[offset + 25]) / (I_Spintau);

      // iNaspin=-gNaSpin*M_Spin^2*H_Spin*I_Spin*(vspin-vNa)
      return -1.0 * g_bar_Na * y[offset + 11] * y[offset + 11] * y[offset + 12] * y[offset + 25] * (y[offset] - E_Na);
    }
  else
    {
      // alphams(v)=  0.32*(-46.9-v)/(exp((-46.9-v)/4.0)-1.0)
      alphams = 0.32 * (-46.9 - y[offset])/(exp( (-46.9 - y[offset]) / 4.0) - 1.0);

      // betams(v)=   0.28*(v+19.9)/(exp((v+19.9)/5.0)-1.0)
      betams = 0.28 * (y[offset] + 19.9)/(exp((y[offset] + 19.9) / 5.0) - 1.0);

      // alphahs(v)=  0.128*exp((-43.0-v)/18.0)
      alphahs = 0.128 * exp((-43.0 - y[offset]) / 18.0);

      // betahs(v)=   4.0/(1.0+exp((-20.0-v)/5.0))
      betahs = 4.0 / (1.0 + exp((-20.0 - y[offset]) / 5.0));
      
      // M_Soma(v)=alphams(v)/(alphams(v)+betams(v))
      M_Soma = alphams / (alphams + betams);

      // H_Somainf=1/(1+exp((vSoma+49)/3.5))
      H_Somainf = 1.0 / (1.0 + exp((y[offset] + 49.0) / 3.5));

      // H_Somatau=1
      H_Somatau = 1.0;

      // no M_Soma' in code
      f[offset + 11] = 0.0;

      // H_Soma'=alphahs(vSoma)-(alphahs(vSoma)+betahs(vSoma))*H_Soma
      f[offset + 12] = alphahs - (alphahs + betahs) * y[offset + 12];

      // no I for soma
      f[offset + 25] = 0.0;

      // iNaSoma=-gNaSoma*M_Soma(vSoma)^2*H_Soma*(vSoma-vNa)
      return -1.0 * g_bar_Na * M_Soma * M_Soma * y[offset + 12] * (y[offset] - E_Na);
    }
}

double Kdr_current(struct network * network, long num_neuron, long num_compartment, double * f, const double * y, double t)
{
  double alpha_n, beta_n, tau_n, n_inf;
  double g_bar_K, E_K;
  long offset = calc_y_array_offset(network, num_neuron, num_compartment);

  g_bar_K = y[offset + 4];
  E_K = network->neurons[num_neuron]->params->values[2];

  /*
  // potassium, from Dayan & Abbott
  alpha_n = 0.02*(y[offset] + 45.7)/(1.0 - exp(-0.1*(y[offset]+45.7)));
  beta_n = 0.25*exp(-0.0125*(y[offset]+55.7));
  tau_n = 1.0/(alpha_n + beta_n);
  n_inf = tau_n * alpha_n;
  
  f[offset + 10] = (n_inf - y[offset + 10])/(tau_n);

  return -1.0*(g_bar_K * pow(y[offset + 10],4.0) * (y[offset] - E_K));
  */

  // from Rubin et al. 2005
  double N_Somainf, N_Spininf, N_Somatau, N_Spintau, alphans, betans;

  // compartment 0 is the spine
  if(num_compartment == 0)
    {
      // N_Spininf=1/(1+exp((-vSpin-42)/2))
      N_Spininf = 1.0/(1.0 + exp((y[offset] + 42.0) / -2.0));

      // N_Spintau=2.2
      N_Spintau = 2.2;

      // N_Spin'=(N_Spininf-N_Spin)/(N_Spintau)
      f[offset + 10] = (N_Spininf - y[offset + 10]) / (N_Spintau);

      // iKdrSpin=-gKdrSpin*N_Spin^2*(vspin-vK)
      return -1.0 * g_bar_K * y[offset + 10] * y[offset + 10] * (y[offset] - E_K);
    }
  else
    {
      // alphans(v)=  0.016*(-24.9-v)/(exp((-24.9-v)/5.0)-1.0)
      alphans = 0.016*(-24.9-y[offset])/(exp((-24.9-y[offset])/5.0)-1.0);

      // betans(v)=   0.25*exp(-1.0-0.025*v)
      betans =  0.25*exp(-1.0-0.025*y[offset]);

      // N_Somainf=1/(1+exp((-vSoma-46.3)/3))
      N_Somainf = 1.0/(1.0 + exp((y[offset] + 46.3) / -3.0));

      // N_Somatau=3.5
      N_Somatau = 3.5;

      // N_Soma'=alphans(vSoma)-(alphans(vSoma)+betans(vSoma))*N_Soma
      f[offset + 10] = alphans - (alphans + betans) * y[offset + 10];

      // iKdrSoma=-gKdrSoma*N_Soma*(vSoma-vK)
      return -1.0 * g_bar_K * y[offset + 10]  * (y[offset] - E_K);
    }

}

double zeta(double v)
{
  // par Btaumod=7 inact=72 inact2=0.11 inact3=2 inact4=64 inact5=1 zetap=30
  double zetap = 30.0;

  // zeta(v)=-1.5-(1/(1+exp((v+zetap)/5)))
  return -1.5 - (1.0 / (1.0 + exp((v + zetap) / 5.0)));
}

double zeta2(double v)
{
  // zeta2(v)=-1.8-(1/(1+exp((v+40)/5)))
  return -1.8 - (1.0 / (1.0 + exp((v + 40.0) / 5.0)));
}


double A_current(struct network * network, long num_neuron, long num_compartment, double * f, const double * y, double t)
{
  double a_inf, tau_a, b_inf, tau_b;
  double g_bar_A, E_A;
  long offset = calc_y_array_offset(network, num_neuron, num_compartment);

  g_bar_A = y[offset + 3];
  E_A = network->neurons[num_neuron]->params->values[3];

  /*
  // A-current, from Dayan & Abbott
  a_inf = pow((0.0761 * exp(0.0314 * (y[offset] + 94.22))) / (1.0 + exp(0.0346 * (y[offset] + 1.17))), (1.0/3.0));
  tau_a = 0.3632 + 1.158 / (1.0 + exp(0.0497 * (y[offset] + 55.96)));
  
  b_inf = pow((1.0 / (1.0 + exp(0.0688 * (y[offset] + 53.3)))),4.0);
  tau_b = 1.24 + 2.678 / (1.0 + exp(0.0624 * (y[offset] + 50.0)));

  f[offset + 13] = (a_inf - y[offset + 13])/(tau_a);
  f[offset + 14] = (b_inf - y[offset + 14])/(tau_b);
  
  return -1.0*(g_bar_A * pow(y[offset + 13],3.0) * y[offset + 14] * (y[offset] - E_A));
  */

  double A_SpnAlf, A_SpnBet, A_Spininf, A_Spintau, B_Spininf, B_Spintau, A_SmaAlf, A_SmaBet, A_Somainf, A_Somatau, B_Somainf, B_Somatau;
  double tempC = network->neurons[num_neuron]->params->values[27];

  // Q=96480/(8.315*(273.16+tempC))
  double Q = 96480.0 / (8.315 * (273.16 + tempC));

  // QT=5^((tempC-24)/10)
  double QT = pow(5.0,(tempC-24.0)/10.0);

  // par asap=0.001
  double asap = 0.001;

  // par Btaumod=7 inact=72 inact2=0.11 inact3=2 inact4=64 inact5=1 zetap=30
  double btaumod=7, inact=72, inact2=0.11, inact3=2, inact4=64, inact5=1, zetap=30;


  // compartment 0 is the spine
  if(num_compartment == 0)
    {
      // A_SpnAlf=exp(asap*zeta(vSpin)*(vSpin+1)*Q)
      A_SpnAlf = exp(asap*zeta(y[offset])*(y[offset]+1.0)*Q);

      // A_SpnBet=exp(0.00039*Q*(vSpin+1)*zeta2(vSpin))
      A_SpnBet = exp(0.00039*Q*(y[offset]+1.0)*zeta2(y[offset]));

      // A_Spininf=1/(1+A_SpnAlf)
      A_Spininf = 1.0/(1.0 + A_SpnAlf);

      // A_Spintau=max(A_SpnBet/((1+A_SpnAlf)*QT*0.1),0.1)
      A_Spintau = max(A_SpnBet/((1.0+A_SpnAlf)*QT*0.1),0.1);
  
      // B_Spininf=0.3+0.7/(1+exp(inact2*(vSpin+inact)*Q))
      B_Spininf = 0.3 + 0.7/(1.0 + exp(inact2 * (y[offset] + inact)*Q));

      // B_Spintau=btaumod*max(inact3*(vSpin+inact4),inact5)
      B_Spintau = btaumod * max(inact3 * (y[offset] + inact4),inact5);

      // A_Spin'=(A_Spininf-A_Spin)/(A_Spintau)
      f[offset + 13] = (A_Spininf - y[offset + 13]) / (A_Spintau);

      // B_Spin'=(B_Spininf-B_Spin)/(B_Spintau)
      f[offset + 14] = (B_Spininf - y[offset + 14]) / (B_Spintau);
    }
  else
    {
      // A_SmaAlf=exp(0.001*zeta(vSoma)*(vSoma-11)*Q)
      A_SmaAlf = exp(0.001*zeta(y[offset])*(y[offset]-11.0)*Q);

      // A_SmaBet=exp(0.00055*Q*(vSoma-11)*zeta(vSoma))
      A_SmaBet = exp(0.00055*Q*(y[offset]-11.0)*zeta(y[offset]));

      // A_Somainf=1/(1+A_SmaAlf)
      A_Somainf = 1.0/(1.0+A_SmaAlf);

      // A_Somatau=max(A_SmaBet/((1+A_SmaAlf)*QT*0.05),0.1)
      A_Somatau = max(A_SmaBet/((1.0+A_SmaAlf)*QT*0.05),0.1);

      // B_Somainf=0.3+0.7/(1+exp(0.02*(vSoma+63.5)*Q))
      B_Somainf = 0.3 + 0.7/(1.0 + exp(0.02*(y[offset] + 63.5)*Q));

      // B_Somatau=btaumod*max(0.11*(vSoma+62),2)
      B_Somatau = btaumod * max(0.11*(y[offset]+62.0),2.0);

      // A_Soma'=(A_Somainf-A_Soma)/(A_Somatau)
      f[offset + 13] = (A_Somainf - y[offset + 13])/(A_Somatau);

      // B_Soma'=(B_Somainf-B_Soma)/(B_Somatau)
      f[offset + 14] = (B_Somainf - y[offset + 14])/(B_Somatau);
    }

  // iKaSoma=-gKaSoma*A_Soma*B_Soma*(vSoma-vK)
  // iKaSpin=-gKaSpin*A_Spin*B_Spin*(vspin-vK)
  return -1.0 * g_bar_A * y[offset + 13] * y[offset + 14] * (y[offset] - E_A);

}

double KCa_current(struct network * network, long num_neuron, long num_compartment, double * f, const double * y, double t)
{
  double c_inf, tau_c;
  double g_bar_KCa, E_K;
  double Q, Qm_Smainf, Qm_SmaAlf, Qm_SmaBet, Qm_Smatau;
  // par qma=.00048,qmb=.28,qhat=1
  double qma = .00048 , qmb = 0.28, qhat = 1.0;
  double temperature;
  long offset = calc_y_array_offset(network, num_neuron, num_compartment);

  g_bar_KCa = y[offset + 5];
  E_K = network->neurons[num_neuron]->params->values[2];
  temperature = network->neurons[num_neuron]->params->values[27];

  /*
  // calcium-dependent potassium current, from Dayan & Abbott
  c_inf = y[offset + 18] / (y[offset + 18] + 3.0) * (1.0 / (1.0 + exp((y[offset] + 28.3) / -12.6)));
  tau_c = 90.3 - 75.1 / (1.0 + exp((y[offset] + 46.0) / -22.7));
  
  f[offset + 15] = (c_inf - y[offset + 15])/(tau_c);

  return -1.0*(g_bar_KCa * pow(y[offset + 15],4.0) * (y[offset] - E_K));
  */
  
  // calcium-dependent potassium current, from Rubin et al. 2005
  if(num_compartment == 0)
    {
      f[offset + 15] = 0.0;
      return 0.0;
    }
  else
    {
      // Q=96480/(8.315*(273.16+tempC))
      Q = 96480.0 / (8.315 * (273.16 + temperature));

      // Qm_SmaAlf=qma*Chisoma/(0.001*Chisoma+0.18*exp(-1.68*vSoma*Q))
      Qm_SmaAlf = qma * y[offset + 18] / (0.001 * y[offset + 18] + 0.18 * exp(-1.68 * y[offset] * Q));

      // Qm_SmaBet=(qmb*exp(-0.022*vSoma*Q))/(exp(-0.022*vSoma*Q)+0.001*Chisoma)
      Qm_SmaBet = (qmb * exp(-0.022 * y[offset] * Q)) / (exp(-0.022 * y[offset] * Q) + 0.001 * y[offset + 18]);
      
      // Qm_Smatau=1/(Qm_SmaAlf+Qm_SmaBet)
      Qm_Smatau = 1.0 / (Qm_SmaAlf + Qm_SmaBet);

      // Qm_Smainf=qhat*Qm_SmaAlf*Qm_Smatau
      Qm_Smainf = qhat * Qm_SmaAlf * Qm_Smatau;
      
      // Qm_Soma'=(Qm_Smainf-Qm_Soma)/(Qm_Smatau)
      f[offset + 15] = (Qm_Smainf - y[offset + 15]) / (Qm_Smatau);
      
      // iKahpmSma=-gKmahpSma*Qm_Soma*(vSoma-vK)
      return -1.0 * g_bar_KCa * y[offset + 15] * (y[offset] - E_K);
    }
  
}

double CaT_current(struct network * network, long num_neuron, long num_compartment, double * f, const double * y, double t)
{
  double S, T;
  double alpha_S, beta_S, S_inf, tau_S;
  double T_inf, tau_T;
  double g_bar_CaT, E_Ca, temperature;
  double xx, ghk;
  long offset = calc_y_array_offset(network, num_neuron, num_compartment);

  g_bar_CaT = y[offset + 6];
  E_Ca = network->neurons[num_neuron]->params->values[5];
  temperature = network->neurons[num_neuron]->params->values[27];

  // transient calcium current as mentioned in Poirazi et al. 2003 supplement but as defined in Dayan & Abbot
  /*
  double M_inf, tau_M, H_inf, tau_H;

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
  */

  // CaL from Rubin et al. 2005
  if(num_compartment == 0)
    {
      // S_Spininf = 1 / (1 + exp(-vSpin - 37))
      S_inf = 1.0 / (1.0 + exp(-1.0 * y[offset] - 37.0));

      // par s1 = 0, s2 = 40, s3 = 3.6
      // S_Spintau = s3 + s1/(1 + exp(vSpin + s2))
      // S_Spintau = 3.6 + 0/(1 + exp(vSpin + 40))
      tau_S = 3.6 + 0.0 / (1.0 + exp(y[offset] + 40.0));

      // T_Spininf = 1 / (1 + exp((vSpin + 41)/.5))
      T_inf = 1.0 / (1.0 + exp((y[offset] + 41.0) / 0.5));

      // T_Spintau = 29
      tau_T = 29.0;

      // S_Spin' = (S_Spininf - S_Spin) / (S_Spintau)
      // T_Spin' = (T_Spininf - T_Spin) / (T_Spintau)
      f[offset + 16] = (S_inf - y[offset + 16]) / tau_S;
      f[offset + 17] = (T_inf - y[offset + 17]) / tau_T;

      // iCaLSpin = -gCaLSpin * S_Spin^3 * T_Spin * (vSpin - vCa)
      return -1.0 * g_bar_CaT * pow(y[offset + 16], 3.0) * y[offset + 17] * (y[offset] - E_Ca);
    }
  else
    {
      //Salfa(v)=-0.055*(v+27.01)/(exp((-v-27.01)/3.8)-1)
      alpha_S = -0.055 * (y[offset] + 27.01) / (exp((y[offset] + 27.01)/-3.8) - 1.0);

      //Sbeta(v)=0.94*exp((-v-63.01)/17)
      beta_S = 0.94 * exp((y[offset] + 63.01) / -17.0);

      //S_Somainf=Salfa(vSoma)/(Salfa(vSoma)+Sbeta(vSoma))
      S_inf = alpha_S / (alpha_S + beta_S);

      //S_Somatau=1/(5*(Salfa(vSoma)+Sbeta(vSoma)))
      tau_S = 1.0 / (5.0 * (alpha_S + beta_S));

      // T doesn't seem to be used for the soma
      //T_inf = 1.0 / (1.0 + exp((y[offset] + 41.0) / 0.5));
      //tau_T = 29.0;
      T_inf = 0;
      tau_T = 0;

      // S_Soma'=(S_Somainf-S_Soma)/(S_Somatau)
      f[offset + 16] = (S_inf - y[offset + 16]) / tau_S;

      //f[offset + 17] = (T_inf - y[offset + 17]) / tau_T;
      f[offset + 17] = 0;

      //xx=0.0853*(273.16+tempC)/2
      xx = 0.0853 * (273.16 + temperature) / 2.0;

      //ghk(v,chi)=-xx*(1-((chi/Ca)*exp(v/xx)))*eff(v/xx)
      ghk = -1.0 * xx * (1.0 - ((y[offset + 18] / 2.0) * exp(y[offset] / xx))) * eff(y[offset] / xx);

      //iCaLSoma=-gCaLSoma*S_Soma*ghk(vSoma,chiSoma)*(1/(1+chisoma))
      return -1.0 * g_bar_CaT * y[offset + 16] * ghk / (1.0 + y[offset + 18]);
    }
}

double eff(double z)
{
  //eff(z)=(1-z/2)*eff2(z)+(z/(exp(z)-1))*eff3(z)
  return (1.0 - z / 2.0) * eff2(z) + (z / (exp(z) - 1.0)) * eff3(z);
}

double eff2(double z)
{
  //eff2(z)=heav(0.0001-abs(z))
  return heav(0.0001 - abs(z));
}

double eff3(double z)
{
  //eff3(z)=heav(abs(z)-0.0001)
  return heav(abs(z) - 0.0001);
}

double heav(double x)
{
  if(x < 0)
    return 0.0;
  return 1.0;
}

double L_current(struct network * network, long num_neuron, long num_compartment, double * f, const double * y, double t)
{
  double g_bar_L, E_L;
  long offset = calc_y_array_offset(network, num_neuron, num_compartment);

  g_bar_L = y[offset + 1];
  E_L = network->neurons[num_neuron]->params->values[4];

  // iLeakSoma=-gLeakSoma*(vSoma-vLeak)
  // iLeakSpin=-gLeakSpin*(vspin-vLeak)
  return -1.0 * g_bar_L * (y[offset] - E_L);
}

double NMDA_current(struct network * network, long num_neuron, long num_compartment, double * f, const double * y, double t, double i_CaT)
{
  double m_NMDA_Ca, m_NMDA_syn, s_NMDA_rise, s_NMDA_fast, s_NMDA_slow, s_NMDA;
  double Phi, beta_soma, Ca_resting_conc_soma, Ca_diffusion_rate, eta_buff;
  double i_Ca, i_Ca_NMDA, beta_spine, Ca_resting_conc_spine, beta_buff;
  double Mg_conc, g_bar_NMDA_Ca, g_bar_NMDA_syn, E_NMDA_syn, E_Ca;
  double Phi_NMDA, tau_NMDA_rise, tau_NMDA_fast, tau_NMDA_slow, tau_Ca;
  double f_pre;

  long offset = calc_y_array_offset(network, num_neuron, num_compartment);

  g_bar_NMDA_Ca = y[offset + 7];
  g_bar_NMDA_syn = y[offset + 8];
  E_Ca = network->neurons[num_neuron]->params->values[5];
  E_NMDA_syn = network->neurons[num_neuron]->params->values[6];
  Mg_conc = network->neurons[num_neuron]->params->values[17];

  Phi_NMDA = network->neurons[num_neuron]->params->values[8];
  tau_NMDA_rise = network->neurons[num_neuron]->params->values[9];
  tau_NMDA_fast = network->neurons[num_neuron]->params->values[10];
  tau_NMDA_slow = network->neurons[num_neuron]->params->values[11];
  tau_Ca = network->neurons[num_neuron]->params->values[16];
  Ca_resting_conc_soma = network->neurons[num_neuron]->params->values[18];
  Ca_resting_conc_spine = network->neurons[num_neuron]->params->values[19];
  Phi = network->neurons[num_neuron]->params->values[20];
  beta_soma = network->neurons[num_neuron]->params->values[21];
  beta_spine = network->neurons[num_neuron]->params->values[22];
  eta_buff = network->neurons[num_neuron]->params->values[23];
  beta_buff = network->neurons[num_neuron]->params->values[24];
  Ca_diffusion_rate = network->neurons[num_neuron]->params->values[25];

  // i_Ca
  i_Ca = i_CaT; // i_CaT passed as argument

  // compartment 0 is spine, 1 is soma
  if(num_compartment == 0)
    {
      // NMDAR, from Rubin et al. 2005

      // par block=0.062 blockCa=0.124
      // mCaNMDA=  1/(1.0+0.3*Mg*exp(-blockCa*vSpin))
      m_NMDA_Ca =  1.0 / (1.0 + 0.3 * Mg_conc * exp(-0.124 * y[offset]));

      // par block=0.062 blockCa=0.124
      // mNMDA=    1/(1.0+0.3*Mg*exp(-block*vSpin))
      m_NMDA_syn =  1.0 / (1.0 + 0.3 * Mg_conc * exp(-0.062 * y[offset]));
      // note: this actually first apperas in destexhe et al 1994

      s_NMDA_rise = y[offset + 19];
      s_NMDA_fast = y[offset + 20];
      s_NMDA_slow = y[offset + 21];

      // sNMDA=sNMDAfast+sNMDAslow+sNMDArise
      s_NMDA = s_NMDA_rise + s_NMDA_fast + s_NMDA_slow;

      // iCaNMDA=-gCaNMDA*sNMDA*mCaNMDA*(vSpin-vCaNMDA)
      i_Ca_NMDA = -1.0 * g_bar_NMDA_Ca * s_NMDA * m_NMDA_Ca * (y[offset] - E_Ca);

      f_pre = presynaptic_activity(network, num_neuron, num_compartment, f, y, t);

      network->neurons[num_neuron]->compartments[num_compartment]->buffer->values[10] = f_pre;
      network->neurons[num_neuron]->compartments[num_compartment]->buffer->values[11] = i_Ca_NMDA;
      network->neurons[num_neuron]->compartments[num_compartment]->buffer->values[12] = m_NMDA_Ca;
      network->neurons[num_neuron]->compartments[num_compartment]->buffer->values[13] = m_NMDA_syn;
      network->neurons[num_neuron]->compartments[num_compartment]->buffer->values[14] = s_NMDA_rise;
      network->neurons[num_neuron]->compartments[num_compartment]->buffer->values[15] = s_NMDA_fast;
      network->neurons[num_neuron]->compartments[num_compartment]->buffer->values[16] = s_NMDA_slow;

      // iCaSpin=iCaLSpin
      // par nonc=6
      // par chiRstSo=0.05 chiRstSp=0.07 PhiSoma=0.1 PhiSpin=0.1 betaSoma=0.083 betaSpin=0.083 CaTauSoSp=1000
      // par st=-100, buff=0
      // ChiSpin' = PhiSpin*(iCaSpin+iCaNMDA)-betaspin*(ChiSpin-ChiRstSp)-(betaspin/nonc)*ChiSpin^2-buff*ChiSpin
      f[offset + 18] = Phi * (i_Ca + i_Ca_NMDA) - beta_spine * (y[offset + 18] - Ca_resting_conc_spine) - beta_spine / eta_buff * pow(y[offset + 18],2.0) - beta_buff * y[offset + 18];

      // par speedup=20 nmdarate=2 ndf=10 nds=45
      // sNMDArise'=-speedup*(1-sNMDAfast-sNMDAslow)*fpre(t)-(1/nmdarate)*sNMDArise
      f[offset + 19] = -1.0 * Phi_NMDA * (1.0 - s_NMDA_fast - s_NMDA_slow) * f_pre - s_NMDA_rise / tau_NMDA_rise;

      // sNMDAfast'=speedup*(0.527-sNMDAfast)*fpre(t)-(1/ndf)*sNMDAfast
      f[offset + 20] = Phi_NMDA * (0.527 - s_NMDA_fast) * f_pre - s_NMDA_fast / tau_NMDA_fast;

      // sNMDAslow'=speedup*(0.473-sNMDAslow)*fpre(t)-(1/nds)*sNMDAslow
      f[offset + 21] = Phi_NMDA * (0.473 - s_NMDA_slow) * f_pre - s_NMDA_slow / tau_NMDA_slow;

      // iNMDA=-gNMDA*sNMDA*mNMDA*(vSpin-vSynapRev)
      return -1.0 * g_bar_NMDA_syn * s_NMDA * m_NMDA_syn * (y[offset] - E_NMDA_syn);
    }
  else
    {

      network->neurons[num_neuron]->compartments[num_compartment]->buffer->values[10] = 0;
      network->neurons[num_neuron]->compartments[num_compartment]->buffer->values[11] = 0;
      network->neurons[num_neuron]->compartments[num_compartment]->buffer->values[12] = 0;
      network->neurons[num_neuron]->compartments[num_compartment]->buffer->values[13] = 0;
      network->neurons[num_neuron]->compartments[num_compartment]->buffer->values[14] = 0;
      network->neurons[num_neuron]->compartments[num_compartment]->buffer->values[15] = 0;
      network->neurons[num_neuron]->compartments[num_compartment]->buffer->values[16] = 0;

      // iCaSoma=iCaLSoma
      // par chiRstSo=0.05 chiRstSp=0.07 PhiSoma=0.1 PhiSpin=0.1 betaSoma=0.083 betaSpin=0.083 CaTauSoSp=1000
      // ChiSoma' = PhiSoma*(iCaSoma)-(betaSoma*(ChiSoma-ChiRstSo))+(ChiSpin-ChiSoma)/CaTauSoSp-(betaSoma/nonc)*ChiSoma^2
      f[offset + 18] = Phi * i_Ca - beta_soma * (y[offset + 18] - Ca_resting_conc_soma) + (y[offset + 18 - network->neurons[num_neuron]->compartments[0]->state->num_params] - y[offset + 18]) / Ca_diffusion_rate - beta_soma / eta_buff * pow(y[offset + 18], 2.0);
      f[offset + 19] = 0.0;
      f[offset + 20] = 0.0;
      f[offset + 21] = 0.0;
      return 0;
    }
}

double AMPA_current(struct network * network, long num_neuron, long num_compartment, double * f, const double * y, double t)
{
  double s_AMPA_rise, s_AMPA_fast, s_AMPA_slow, s_AMPA;
  double Phi_AMPA, tau_AMPA_rise, tau_AMPA_fast, tau_AMPA_slow;
  double g_bar_AMPA, E_AMPA_syn, f_pre;
  long offset = calc_y_array_offset(network, num_neuron, num_compartment);

  if(num_compartment == 0)
    {
      g_bar_AMPA = y[offset + 9];
      E_AMPA_syn = network->neurons[num_neuron]->params->values[7];
      f_pre = presynaptic_activity(network, num_neuron, num_compartment, f, y, t);
      
      // AMPAR, from Rubin et al. 2005
      s_AMPA_rise = y[offset + 22];
      s_AMPA_fast = y[offset + 23];
      s_AMPA_slow = y[offset + 24];

      // sAMPA=sAMPAfast+sAMPAslow+sAMPArise
      s_AMPA = s_AMPA_rise + s_AMPA_fast + s_AMPA_slow;
      
      Phi_AMPA = network->neurons[num_neuron]->params->values[12];
      tau_AMPA_rise = network->neurons[num_neuron]->params->values[13];
      tau_AMPA_fast = network->neurons[num_neuron]->params->values[14];
      tau_AMPA_slow = network->neurons[num_neuron]->params->values[15];
      
      network->neurons[num_neuron]->compartments[num_compartment]->buffer->values[17] = s_AMPA_rise;
      network->neurons[num_neuron]->compartments[num_compartment]->buffer->values[18] = s_AMPA_fast;
      network->neurons[num_neuron]->compartments[num_compartment]->buffer->values[19] = s_AMPA_slow;

      // par speedup=20 nmdarate=2 ndf=10 nds=45
      // sAMPArise'=-speedup*(1-sAMPAfast-sAMPAslow)*fpre(t)-(1/0.58)*sAMPArise
      f[offset + 22] = -1.0 * Phi_AMPA * (1.0 - s_AMPA_fast - s_AMPA_slow) * f_pre - s_AMPA_rise / tau_AMPA_rise;

      // sAMPAfast'=speedup*(0.903-sAMPAfast)*fpre(t)-(1/7.6)*sAMPAfast
      f[offset + 23] = Phi_AMPA * (0.903 - s_AMPA_fast) * f_pre - s_AMPA_fast / tau_AMPA_fast;

      // sAMPAslow'=speedup*(0.097-sAMPAslow)*fpre(t)-(1/25.69)*sAMPAslow
      f[offset + 24] = Phi_AMPA * (0.097 - s_AMPA_slow) * f_pre - s_AMPA_slow / tau_AMPA_slow;
      
      // iAMPA=-allowAMPA*gAMPA*sAMPA*(vSpin-vSynapRev)
      return -1.0*g_bar_AMPA * s_AMPA  * (y[offset] - E_AMPA_syn);
    }
  else
    {

      network->neurons[num_neuron]->compartments[num_compartment]->buffer->values[17] = 0;
      network->neurons[num_neuron]->compartments[num_compartment]->buffer->values[18] = 0;
      network->neurons[num_neuron]->compartments[num_compartment]->buffer->values[19] = 0;

      f[offset + 22] = 0.0;
      f[offset + 23] = 0.0;
      f[offset + 24] = 0.0;
      return 0.0;
    }
}

double square_wave_f_pre(struct network * network, long num_neuron, long num_compartment, double * f, const double * y, double t)
{
  long i, from_neuron, from_compartment;
  double presyn_V;

  for(i = 0; i < network->neurons[num_neuron]->compartments[num_compartment]->num_links; i++)
    {
      from_neuron = network->neurons[num_neuron]->compartments[num_compartment]->links[i]->from;
      from_compartment = network->neurons[num_neuron]->compartments[num_compartment]->links[i]->from_compartment;

      //presyn_V = network->neurons[from_neuron]->compartments[from_compartment]->state->values[0];
      long from_offset = calc_y_array_offset(network, from_neuron, from_compartment);
      presyn_V = y[from_offset];

      if(presyn_V >= 38.0)
	{
	  return 1.0;
	}
    }
  return 0.0;
}

double presynaptic_activity(struct network * network, long num_neuron, long num_compartment, double * f, const double * y, double t)
{
  struct stimulus * stimulus = NULL;
  if((stimulus = apply_stimulus(network, num_neuron, num_compartment, t)) != NULL && stimulus->direct == 0)
    return 1.0;

  //return square_wave_f_pre(network, num_neuron, num_compartment, f, y, t);
  return destexhe_transmitter_release(network, num_neuron, num_compartment, f, y, t);
}

double destexhe_transmitter_release(struct network * network, long num_neuron, long num_compartment, double * f, const double * y, double t)
{
  // based on Destexhe et al. 1994, "Simplification of the Release Process" section, p. 208
  long from_offset, offset = calc_y_array_offset(network, num_neuron, num_compartment);
  double L = 0.0;
  double L_max = network->neurons[num_neuron]->params->values[61];
  double V_p = network->neurons[num_neuron]->params->values[60];
  double K_p = network->neurons[num_neuron]->params->values[59];
  double V_pre = 0.0;
  long i, from_neuron, from_compartment;

  for(i = 0; i < network->neurons[num_neuron]->compartments[num_compartment]->num_links; i++)
    {
      from_neuron = network->neurons[num_neuron]->compartments[num_compartment]->links[i]->from;
      from_compartment = network->neurons[num_neuron]->compartments[num_compartment]->links[i]->from_compartment;
      from_offset = calc_y_array_offset(network, from_neuron, from_compartment);
      V_pre = y[from_offset];
      L += L_max / (1.0 + exp(-(V_pre - V_p)/K_p));
    }
  
  L = min(L, L_max);
  network->neurons[num_neuron]->compartments[num_compartment]->buffer->values[10] = L;

  return L;
}

double dayan_abbott_f_pre(struct network * network, long num_neuron, long num_compartment, double * f, const double * y, double t)
{
  //return square_wave(network, num_neuron, num_compartment, t);
  // compartment 0 is the spine
  if(num_compartment != 0)
    return 0.0;

  long offset = calc_y_array_offset(network, num_neuron, num_compartment);
  double input = square_wave_f_pre(network, num_neuron, num_compartment, f, y, t);
  double alpha_s = network->neurons[num_neuron]->params->values[57];
  double beta_s = network->neurons[num_neuron]->params->values[58];

  if(input > 0.1)
    alpha_s = network->neurons[num_neuron]->params->values[56];

  f[offset + 26] = alpha_s*(1 - y[offset + 26]) - beta_s*y[offset + 26];

  network->neurons[num_neuron]->compartments[num_compartment]->buffer->values[10] = y[offset + 26];

  return y[offset + 26];
}

long calc_y_array_offset(struct network * network, long num_neuron, long num_compartment)
{
  return network->neurons[num_neuron]->compartments[num_compartment]->ode_system_offset;
}
