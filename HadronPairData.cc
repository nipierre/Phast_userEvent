
#include "HadronPairData.h"
#include "math.h"

HadronPairData::HadronPairData()
{
  Reset();
}


void HadronPairData::Reset()
{
  inPart=false;
  for(int i=0; i<3; ++i){
    primVertPos[i] = vertPos[i] =  p1[i] = p2[i] = 0.0;
  }

  tracks = tracks_up = tracks_down = 0;
  distToPV = sigmaPV = sigmaSV = 0.0;
  P1 = P2 = thetaR1 = thetaR2 = 0.0;
  P1_RICH = P2_RICH = 0.0;
  Q1 = Q2 = 0;
  Y = Q_2 = Xbj = W = 0.0;
  Mmiss = 0.0;
  extrapolation_pos_x_1 = extrapolation_pos_y_1 =0.0;
  extrapolation_pos_x_2 = extrapolation_pos_y_2 =0.0;
  //extrapolation_pos_x_center_1 = extrapolation_pos_y_center_1 =0.0;
  //extrapolation_pos_x_center_2 = extrapolation_pos_y_center_2 =0.0;
  photons1 = photons2=0.0;
  theta_CH1 = theta_CH2 =0.0;
  beta1_pi  = beta1_k  = beta1_p=0.0;
  beta2_pi  = beta2_k  = beta2_p=0.0;
  theta_CM1 = theta_CM2 =0.0;
  Phi1_vtx = Phi2_vtx =0.0; 
  Phi1 = Phi2 =0.0; 
  index_UV = index_vis =0.0;
  richInfo1 = richInfo2 = excl_phi= false;
  for(int i=0; i<6; ++i){
    LHs1[i] = LHs2[i] = -10.0;
  }
}


float HadronPairData::InvMass(const float& m) const
{

  float E1 = sqrt(m*m + P1*P1);
  float E2 = sqrt(m*m + P2*P2);

  float Etot = E1 + E2;
  float s = Etot*Etot;
  float ptot[3]; 
  for(int i=0; i<3; ++i){
    ptot[i] = p1[i] + p2[i];
    s -= ptot[i]*ptot[i];
  }
  
  return sqrt(s);
}

float HadronPairData::CosThetaVK() const
{
  if( inPart || distToPV==0 ) return 1.;

  float r_i,p_i,P=0.0,cosTh=0.0;
  for(int i=0; i<3; ++i){
    r_i = vertPos[i] - primVertPos[i];
    p_i = p1[i] + p2[i];
    P += p_i*p_i;
    cosTh += r_i*p_i;
  }
  cosTh /= distToPV*sqrt(P);
  return cosTh;
}


int HadronPairData::MaxLHMassHypo1() const
{
  return MaxLHMassHypo(richInfo1,LHs1);
}
int HadronPairData::MaxLHMassHypo2() const
{
  return MaxLHMassHypo(richInfo2,LHs2);
}
int HadronPairData::MaxLHMassHypo(const bool& richInfo,const float* LH) const
{
  if( !richInfo ) return -1;
  int max=0;
  for(int i=1; i<6; ++i) if( LH[i] > LH[max] ) max = i;
  return max;
}


bool HadronPairData::MaxLHPion(const int& part) const
{
  return MaxLHHypo(part, 0);
}
bool HadronPairData::MaxLHKaon(const int& part) const
{
  return MaxLHHypo(part, 1);
}
bool HadronPairData::MaxLHElectron(const int& part) const
{
  return MaxLHHypo(part, 3);
}
bool HadronPairData::MaxLHHypo(const int& part,const int& hypo) const
{
  if( part < 1 || part > 2 ) return false;
  int max = -1;
  switch (part){
  case 1: max = MaxLHMassHypo1(); break;
  case 2: max = MaxLHMassHypo2(); break;
  }
  return ( max == hypo ); 
}
