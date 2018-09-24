#include "DISEventData.h"
#include "TMath.h"
#include "G3part.h"

DISEventData::DISEventData():
  runNo(0),spillNo(0),evtInSpill(0),
  trigMask(0),
  evNo(0),timeInSpill(0),
  x(0),y(0),z(0),       // vertex position
  p0x(0),p0y(0),p0z(0), // beam momentum
  p1x(0),p1y(0),p1z(0), // mu1 momentum
  E_beam(0), // beam energie
  E_mu_prim(0), // mu1 energie
  Charge(0),
  XX0(0),  // no. of rad. len. crossed by mu1
  saved(false),
  BPV(0),
  isMuPrim(false),
  MZfirst(0),
  beam_chi2(0),
  mu_prim_chi2(0),
  cellsCrossed(false),
  backPropFlag(false)
{}

void DISEventData::Reset()
{
  runNo=0; spillNo=0; evtInSpill=0;
  trigMask=0;
  evNo=0; timeInSpill=0;
  x=0; y=0; z=0;        // vertex position
  p0x=0; p0y=0; p0z=0;  // beam momentum
  p1x=0; p1y=0; p1z=0;  // mu1 momentum
  E_beam=0; // beam energie
  E_mu_prim=0; // mu1 energie
  Charge=0;
  XX0=0;   // no. of rad. len. crossed by mu1
  saved=false;
  BPV=0;
  isMuPrim=false;
  MZfirst=0;
  beam_chi2=0;
  mu_prim_chi2=0;
  cellsCrossed=false;
  backPropFlag=false;
}

static const double muMass  = G3partMass[5];
static const double muMass2 = muMass*muMass;
static const double pMass   = G3partMass[14]; // proton mas
static const double pMass2  = pMass*pMass;

void DISEventData::CalcKin()
{

  double P0 = sqrt( p0x*p0x + p0y*p0y + p0z*p0z );
  double P1 = sqrt( p1x*p1x + p1y*p1y + p1z*p1z );
  double p0p1 = p0x*p1x + p0y*p1y + p0z*p1z;
  double costh = p0p1 / P0 / P1;

  if (costh>-1 && costh<1) theta = acos(costh);
  else theta = (costh>0)? 0 : TMath::Pi();

  E0  = sqrt( muMass2 + P0*P0 );
  E1  = sqrt( muMass2 + P1*P1 );
  th0 = acos(p0z/P0);
  th1 = acos(p1z/P1);
  nu  = E0 - E1;
  Y   = nu / E0;
  Q2  = 2.*( E0*E1 - p0p1 - muMass2);
  xBj = Q2 / (2. * pMass * nu);
  W2  = pMass2 + 2.*pMass*nu - Q2;
  W   = sqrt( W2 );
}



//--------------- HadronData

HadronData::HadronData():
  P(0),pt(0),th(0),ph(0),ph_pl(0),
  XX0(0),
  inHCALacc(false),
  HCAL(0),
  charge(0),
  thRICH(0),
  thC(0),
  MCpid(-1),
  MM01x(0),MM01y(0),
  MM02x(0),MM02y(0),
  MM03x(0),MM03y(0),
  Z2Ax(0),Z2Ay(0),
  Z2Bx(0),Z2By(0),
  RICHx(0),RICHy(0), //Extrapolation at RICH by Quiela
  chi2_hadron(0),
  HZfirst(0), HZlast(0),
  dis(0),
  E(0),z(0),
  z_pi(0),z_p(0)
{
  for(int i=0; i<6; ++i) LH[i]=0;
}
HadronData::~HadronData()
{}


void HadronData::Reset()
{
  P=0; pt=0; th=0; ph=0; ph_pl=0;
  XX0=0;
  inHCALacc=false;
  HCAL=0;
  charge=0;
  thRICH=0;
  thC=0;
  MCpid=-1;
  dis=0;
  E=0; z=0;
  z_pi=0; z_p=0;
  for(int i=0; i<6; ++i) LH[i]=0;

  MM01x=MM01y=MM02x=MM02y=MM03x=MM03y=0.0;
  Z2Ax=Z2Ay=Z2Bx=Z2By=0.0;
  RICHx=RICHy=0.0;      //Extrapolation at RICH by Quiela
  chi2_hadron=0;
  HZfirst=HZlast=0;

}


static const double M_pi = G3partMass[8];  // pi+ mass

void HadronData::CalcVariables(const double& Mh)
{
  double M = Mh<0 ? M_pi:Mh;

  E = sqrt(P*P + M*M);
  z = E/dis->nu;
  double E_pi = sqrt(P*P + M_pi*M_pi);
  z_pi = E_pi/dis->nu;
  double E_p = sqrt(P*P + pMass*pMass);
  z_p = E_p/dis->nu;
}



//------------------------------------------- MC data
DISEventMCData::DISEventMCData():
  MC_vx(0),MC_vy(0),MC_vz(0),       // Vertex position
  MC_p0x(0),MC_p0y(0),MC_p0z(0), // beam momentum
  MC_p1x(0),MC_p1y(0),MC_p1z(0), // mu1 momentum
  irad(-1),              // RC flag
  mcWeight(0),              // HEPGEN weight
  w1(0),w2(0),              // RADGEN weights
  difType(0),              // values 1 for dd and 0 for elastic (see documentation of HEPGEN) diffractive evnt
  MC_nuTr(0),MC_Q2Tr(0),      // True photon kinematics
  MC_xTr(0),MC_yTr(0),      // True photon kinematics
  MC_w(1),                 // Weight
  recons(false)
{}

void DISEventMCData::Reset()
{
  MC_vx=0; MC_vy=0; MC_vz=0;        // Vertex position
  MC_p0x=0; MC_p0y=0; MC_p0z=0;  // beam momentum
  MC_p1x=0; MC_p1y=0; MC_p1z=0;  // mu1 momentum
  MC_nuTr=0; MC_Q2Tr=0;  MC_xTr=0; MC_yTr=0;
  MC_w = 1;                // Weight
  recons=false;
}
void DISEventMCData::CalcKin()
{
  float MC_P0 = sqrt( MC_p0x*MC_p0x + MC_p0y*MC_p0y + MC_p0z*MC_p0z );
  float MC_P1 = sqrt( MC_p1x*MC_p1x + MC_p1y*MC_p1y + MC_p1z*MC_p1z );
  float MC_p0p1 = MC_p0x*MC_p1x + MC_p0y*MC_p1y + MC_p0z*MC_p1z;
  float MC_costh = MC_p0p1 / MC_P0 / MC_P1;
  if (MC_costh>-1 && MC_costh<1) MC_theta = acos(MC_costh);
  else MC_theta = (MC_costh>0)? 0 : TMath::Pi();

  MC_E0  = sqrt( muMass2 + MC_P0*MC_P0 );
  MC_E1  = sqrt( muMass2 + MC_P1*MC_P1 );
  MC_th0 = acos(MC_p0z/MC_P0);
  MC_th1 = acos(MC_p1z/MC_P1);
  MC_nu  = MC_E0 - MC_E1;
  MC_Y   = MC_nu / MC_E0;
  MC_Q2  = 2.*( MC_E0*MC_E1 - MC_p0p1 - muMass2);
  MC_xBj = MC_Q2 / (2. * pMass * MC_nu);
  MC_W2  = pMass2 + 2.*pMass*MC_nu - MC_Q2;
  MC_W   = sqrt( MC_W2 );
  /* if (irad) {
    MC_xTr = MC_Q2Tr/2/pMass/MC_nuTr;
    MC_yTr = MC_nuTr/MC_E0;
    }*/
}


HadronMCData::HadronMCData():
  P(0),th(0),ph(0),
  /* XX0(0),*/
  charge(0),
  /*thRICH(0),
  MM01x(0),MM01y(0),
  MM02x(0),MM02y(0),
  MM03x(0),MM03y(0),
  Z2Ax(0),Z2Ay(0),
  Z2Bx(0),Z2By(0),
  RICHx(0),RICHy(0), //Extrapolation at RICH by Quiela*/
  pid(0),
  recons(false),recHadIdx(-1),
  itrack(-1)
{}

void HadronMCData::Reset()
{
  P=0; th=0; ph=0; pid=0; recons=false; recHadIdx=-1;
  itrack=-1;
  E=0.0; z=0.0;
  dis=0;
  /* XX0=0;*/
  charge=0;
  /*thRICH=0;
  MM01x=MM01y=MM02x=MM02y=MM03x=MM03y=0.0;
  Z2Ax=Z2Ay=Z2Bx=Z2By=0.0;
  RICHx=RICHy=0.0;      //Extrapolation at RICH by Quiela*/


}


static const double m_e  = 0.000511;
static const double m_mu = 0.1056583;
static const double m_K  = 0.493677;
static const double m_p  = 0.9382723;
static const double m_gamma = 0.0;

const double& GetMassJetset(const short& id)
{
  if(id == 22 ) return m_gamma; // gammas
  if(id == 11 || id == -11 ) return m_e;
  if(id == 13 || id == -13 ) return m_mu;
  if(id ==  211 || id ==  -211 ) return M_pi;
  if(id ==  321 || id ==  -321 ) return m_K;
  if(id == 2212 || id == -2212 ) return m_p;

  return M_pi;
}

const double& GetMassGeant3(const short& id)
{
  if(id == 1 ) return m_gamma; // gammas
  if(id == 2 || id == 3 ) return m_e;
  if(id == 5 || id == 6 ) return m_mu;
  if(id == 8 || id == 9 ) return M_pi;
  if(id == 11 || id == 12 ) return m_K;
  if(id == 14 || id == 15 ) return m_p;

  return M_pi;
}



void HadronMCData::CalcVariables(const double& Mh)
{
  double M = Mh<0 ? M_pi:Mh;

  E = sqrt(P*P + M*M);
  z = E/dis->MC_nu;
}
