#ifndef HadronPairData_h
#define HadronPairData_h 1

#include "Rtypes.h"

struct HadronPairData{
  
  HadronPairData();
  virtual ~HadronPairData(){};

  bool inPart;

  float primVertPos[3];
  float vertPos[3];
  float extrapolation_pos_x_1;
  float extrapolation_pos_y_1;
  float extrapolation_pos_x_2;
  float extrapolation_pos_y_2;
  //float extrapolation_pos_x_center_1;
  //float extrapolation_pos_y_center_1;
  //float extrapolation_pos_x_center_2;
  //float extrapolation_pos_y_center_2;
  float distToPV;
  float sigmaPV;
  float sigmaSV;
  float Y;
  float Q_2;
  float Xbj;
  float W;
  float Mmiss;
  double index_UV;
  double index_vis;
    
  float p1[3];
  float p2[3];
  float P1,P2;
  float P1_RICH,P2_RICH;
  short Q1,Q2;
  float photons1, photons2;
  float theta_CH1,theta_CH2;
  float beta1_pi, beta1_k, beta1_p;
  float beta2_pi, beta2_k, beta2_p;
  float thetaR1,thetaR2;
  float theta_CM1,theta_CM2;
  float Phi1,Phi2;
  float Phi1_vtx,Phi2_vtx;
  bool excl_phi;
  bool richInfo1,richInfo2;  
  float LHs1[6];
  float LHs2[6];
  int tracks,tracks_up,tracks_down;
  void Reset();

  float InvMass(const float& mass1) const;
  float CosThetaVK() const;
  int MaxLHMassHypo1() const;
  int MaxLHMassHypo2() const;
  int MaxLHMassHypo(const bool& richInfo,const float* LH) const;

  bool MaxLHPion(const int& part) const;
  bool MaxLHKaon(const int& part) const;
  bool MaxLHElectron(const int& part) const;
  bool MaxLHHypo(const int& part,const int& hypo) const;


  ClassDef(HadronPairData,2);
};



#endif
