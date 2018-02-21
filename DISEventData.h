
#ifndef DISEventData_h
#define DISEventData_h 1

#include "Rtypes.h"
#include <vector>

struct DISEventData{

  DISEventData();
  virtual ~DISEventData(){};

  Int_t runNo,spillNo,evtInSpill;
  Int_t trigMask;
  Double_t  evNo,timeInSpill;
  Float_t x,y,z; // vertex position
  Double_t p0x,p0y,p0z; // beam momentum
  Double_t p1x,p1y,p1z; // mu1 momentum
  Double_t E_beam;    // beam energie
  Double_t E_mu_prim; // mu1 energie
  Double_t XX0; // no. of rad. len. crossed by mu1
  Double_t HM04x,HM04y,HM05x,HM05y;
  Double_t HO03x,HO03y,HO04x,HO04y;

  Bool_t saved;
  Bool_t cellsCrossed;
  Bool_t backPropFlag;

  void Reset();
  void CalcKin();

  Double_t     E0; //!
  Double_t     E1; //!
  Double_t    th0; //!
  Double_t    th1; //!
  Double_t  theta; //!
  Double_t     Q2; //!
  Double_t     nu; //!
  Double_t      Y; //!
  Double_t    xBj; //!
  Double_t     W2; //!
  Double_t      W; //!

  ClassDef(DISEventData,5);
};



struct HadronData{

  HadronData();
  virtual ~HadronData();

  Double_t P,th,ph,ph_pl;
  Double_t XX0;
  Bool_t  inHCALacc;
  Float_t HCAL;
  Short_t charge;
  Double_t thRICH;
  Float_t LH[6];

  Short_t MCpid;

  Double_t MM01x,MM01y,MM02x,MM02y,MM03x,MM03y;
  Double_t Z2Ax,Z2Ay,Z2Bx,Z2By;
  Double_t RICHx,RICHy;                       //Extrapolation at RICH by Quiela

  void Reset();
  void CalcVariables(const double& M=-1);

  DISEventData* dis; //!
  Double_t E; //!
  Double_t z; //!
  Double_t z_pi; //!
  Double_t z_p; //!
  Int_t richPID; //!

  ClassDef(HadronData,3);
};



struct DISEventMCData{

  DISEventMCData();
  virtual ~DISEventMCData(){};
  Float_t MC_vx,MC_vy,MC_vz; // vertex position
  Double_t MC_p0x,MC_p0y,MC_p0z; // beam momentum
  Double_t MC_p1x,MC_p1y,MC_p1z; // mu1 momentum
  Short_t irad;       // RC flag:
  Double_t w1; // RADGEN weight
  Double_t w2; // RADGEN weight
  Double_t mcWeight; // HEPGEN weight
  Short_t difType; // values 1 for dd and 0 for elastic (see documentation of HEPGEN) diffractive evnt
  float   MC_nuTr, MC_Q2Tr; // True photon kinematics
  float   MC_yTr, MC_xTr;
  float   MC_w;          // RC weight
  bool recons; // if true the event has been reconstructed

  void Reset();
  void CalcKin();

  Double_t     MC_E0; //!
  Double_t     MC_E1; //!
  Double_t    MC_th0; //!
  Double_t    MC_th1; //!
  Double_t  MC_theta; //!
  Double_t     MC_Q2; //!
  Double_t     MC_nu; //!
  Double_t      MC_Y; //!
  Double_t    MC_xBj; //!
  Double_t     MC_W2; //!
  Double_t      MC_W; //!
  ClassDef(DISEventMCData,2);
};


struct HadronMCData{

  HadronMCData();
  virtual ~HadronMCData(){};

  ///////////added by erin/////////////////////////
  Double_t P,th,ph;
  // Double_t XX0;

  Short_t charge;
  //Double_t thRICH;

  //Double_t MM01x,MM01y,MM02x,MM02y,MM03x,MM03y;
  //Double_t Z2Ax,Z2Ay,Z2Bx,Z2By;
  //Double_t RICHx,RICHy;                       //Extrapolation at RICH by Quiela
  ///////////////////////////////////////////////////
  Short_t pid;

  Bool_t recons; // if true the hadron has been detected
  Short_t recHadIdx; // index of the corresponding reconstructed particle;

  void Reset();
  void CalcVariables(const double& M=-1);

  int itrack; //! used for associating to reconstructed track

  DISEventMCData* dis; //!
  Double_t E; //!
  Double_t z; //!
  Double_t z_pi; //!
  Double_t z_p; //!

  ClassDef(HadronMCData,2)
};

/*
  ClassDef(std::vector\<HadronData\>,1)
*/

/*
  ClassDef(std::vector\<HadronMCData\>,1)
*/


// utility functions
const double& GetMassJetset(const short& id);
const double& GetMassGeant3(const short& id);


#endif
