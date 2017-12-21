#include "LCAnalysis.h"
#include "DISEventData.h"
#include "HadronPairData.h"
#include "TargetCell.h"

#include <stdio.h>
#include "Phast.h"
#include "PaEvent.h"
#include "PaAlgo.h"
#include "PaPid.h"

//#include "RescaleMom.cc"

#include "G3part.h"

#include "TFile.h"
#include "TTree.h"
#include "TEventList.h"
#include "TNtuple.h"
#include "TH1D.h"
#include "TH2D.h"
#include <TVector3.h>
#include <fstream>
#include <sstream>

using namespace std;
////////temp by erin//////
int count_bestpv =0;
int count_mup =0;
int count_mupp=0;
int count_munp=0;
int count_mup0=0;
int count_mun0=0;
///////////////////////
static const double M_mu = G3partMass[5];  // muon mass
static const double M_pi = G3partMass[8];  // pi+ mass
static const double M_K  = G3partMass[11]; // K+ mass
static const double M_p  = G3partMass[14]; // proton mass

static const char CONFFILENAME[]=
  "/afs/cern.ch/user/n/nipierre/workspace/PHAST/user/LC_configuration/lc_conf_file.txt";

static const char CUTFILENAME[]=
  "/afs/cern.ch/user/n/nipierre/workspace/PHAST/user/LC_configuration/lc_cut_file.txt";

static const char CONFFILENAMEVAR[]=
  "/afs/cern.ch/user/n/nipierre/workspace/PHAST/user/LC_configuration/lc_conf_file_%d.txt";

static const char CUTFILENAMEVAR[]=
  "/afs/cern.ch/user/n/nipierre/workspace/PHAST/user/LC_configuration/lc_cut_file_%d.txt";

static const char OUTFILENAME[]=
  "/afs/cern.ch/user/n/nipierre/workspace/HadronSelection/lc_output.root";

static const char LOGFILENAME[]=
  "/afs/cern.ch/user/n/nipierre/workspace/HadronSelection/lc_analysis.log";


static const float  Zrich= 615.6; //entrance of the RICH 580.0

static const float ZHCAL1 = 1195.0;
static const float ZHCAL2 = 3520.0;

void RescaleMom(PaEvent & e, bool faster_mode);




LCAnalysis::LCAnalysis():
  //--------------- utility flags
  fAnalysisType(0),
  fEventAccepted(false),
  fSaveEvents(true),
  //--------------- event selection counters
  fTotNEvts(0),
  fTrigger(0),
  //--------------- default values for kinematical cuts
  fQ2min(1.0),
  fEbeamMin(140.0),fEbeamMax(180.0),
  fyMin(0.1),fyMax(0.9),
  fxBjMin(0.),fxBjMax(1),
  //--------------- hadron selection counters
  fTotNHadr(0),
  fNsec(0),
  fNzlt1(0),
  fNZfirst(0),
  fNZlast(0),
  fNXX0(0),
  fNTrigHad(0),
  fNMChadr(0), fNMChadrMultAssoc(0),
  //--------------- hadron selection variables
  fZref(350.),
  fXX0max(15.),
  //--------------- DIS event info
  fDISEvt(0),
  fDISMCEvt(0),
  fHadronsPtr(&fHadrons),
  fMCHadronsPtr(&fMCHadrons),
  //--------------- Pid utility class object
  fPid(new PaPid),
  //--------------- 2012 target utility class object
  fTcell(new TargetCell),
  //--------------- Pointer initialisations
  fOutFile(0),
  fEvtTree(0),
  fEvtTreeMC(0),
  fHadrTree(0),
  fHadrTreeMC(0),
  fHadrPair(0),
  fHistInfo(0),
  fDISEvtTree(0),
  fDISEvtList(0),
  fUtilityTree(0),
  fdatatargetType(50000),
  fMCtargetType(-1) // -2 two cells (2004), -3 three cells (from 2006) -5 three cells(Yann) -1 one cell (2016)
{
  fEvent=0;

  //--- read config file
  char conffilename[1024];
  if(Phast::Ref().NUserFlag()>1){
    sprintf(conffilename,CONFFILENAMEVAR,Phast::Ref().UserFlag(1));
  }else{
    sprintf(conffilename,CONFFILENAME);
  }
  printf("\n\nLCAnalysis: reading conf. file %s ...\n",conffilename);
  FILE *conffile=fopen(conffilename,"r");
  char line[500];
  while(fgets(line,500,conffile)){
    if(line[0] == '#') continue;
    sscanf(line,"%d",&fAnalysisType);

  }
  fclose(conffile);
  const char anTypeMsg[][256]={"Select events",
			    "Look for hadrons",
			    "Select hadron pairs"};
  printf("LCAnalysis: analysis type is %d - %s\n",fAnalysisType,anTypeMsg[fAnalysisType]);

  if(Phast::Ref().NUserFlag()>0){
    //sprintf(fLogFileName,"/sps/compass/capozza/mc/2006/mcprod-2012-09-14/hadrons/ev_sel_logs/ev-sel-run-%d-%03d.log", PaSetup::Ref().RunNum(),Phast::Ref().UserFlag(0));
    sprintf(fLogFileName,
	    "/sps/compass/seder/hadron_selection/logs/ev-sel-run-%d-%03d.log",
	    PaSetup::Ref().RunNum(),Phast::Ref().UserFlag(0));
  }else{
    sprintf(fLogFileName,LOGFILENAME);
  }
  printf("LCAnalysis: logfile is %s\n\n",fLogFileName);

  // file name for saving the DISEvtTree
  string disEvtFileName(Phast::Ref().h_file->GetName());
  size_t pos=disEvtFileName.find_last_of('/');
  if(pos==string::npos) pos = -1;
  disEvtFileName.insert(pos+1,"disevt_");
  cout<<"LCAnalysis: DISEvtTree file is: "<<disEvtFileName<<endl;



  //--- setting handler pointer
  switch(fAnalysisType){
  case 0:
    fEvHandler = &LCAnalysis::SelectEvent;
    ReadCutFile(); // read cut configuration file lc_cut_file.txt

    //--- book ntuple
    Phast::Ref().h_file->cd();
    fHistInfo = new TNtuple("LCHistInfo","Histogramming info","Q2:radlen:zVTX:TMask");
    fDISEvt = new DISEventData;
    fDISEvtTree = new TTree("DISEvtTree","DIS event info");
    fDISEvtTree->Branch("DISEvt","DISEventData",&fDISEvt);
    fDISEvtTree->SetMaxTreeSize(1000000000);
    fDISEvtList = new TEventList("DISEvtList","List of selected DIS events");
    break;
  case 1:
    fEvHandler = &LCAnalysis::FindHadrons;
    fDISEvt = new DISEventData;

    fOutFile = new TFile(disEvtFileName.data(),"RECREATE");
    //cout<<"file open \t"<<fOutFile<<endl;
    fDISEvtTree = new TTree("DISEvtTree","DIS event and hadron info");
    fDISEvtTree->Branch("DISEvt","DISEventData",&fDISEvt);
    fDISEvtTree->Branch("Hadrons","std::vector<HadronData>",&fHadronsPtr);
    fDISEvtTree->SetMaxTreeSize(1000000000);
    break;
  case 2:
    fSaveEvents = false;
    fEvHandler = &LCAnalysis::SelectHadronPairs;
    ReadCutFile(); // read cut configuration file lc_cut_file.txt
    break;
  default: fEvHandler = &LCAnalysis::SelectEvent;
  }

  cout << "test" << endl;

  //--- get pointers to hadronic calorimeters
  const PaSetup& setup = PaSetup::Ref();
  fHadrCal1 = (PaCalorimeter*)&setup.Calorimeter(setup.iCalorim("HC01P1__"));
  fHadrCal2 = (PaCalorimeter*)&setup.Calorimeter(setup.iCalorim("HC02P1__"));

  cout << "test" << endl;

  //--- find MT/OT hodoscope z positions
  int idet = PaSetup::Ref().iDetector("HM04Y1_d");
  const PaDetect& HM04 =  PaSetup::Ref().Detector(idet);
  HM04z = HM04.Z();

  cout << "test" << endl;

  idet = PaSetup::Ref().iDetector("HM05Y1_d");
  const PaDetect& HM05 =  PaSetup::Ref().Detector(idet);
  HM05z = HM05.Z();


  idet = PaSetup::Ref().iDetector("HO03Y1_m");
  const PaDetect& HO03 =  PaSetup::Ref().Detector(idet);
  HO03z = HO03.Z();


  idet = PaSetup::Ref().iDetector("HO04Y1_m");
  const PaDetect& HO04 =  PaSetup::Ref().Detector(idet);
  HO04z = HO04.Z();


  //--- find MMs z positions
  // values taken from detectors.52959.plus.dat (2006) because not all detectors are present in 2004
  /*idet = PaSetup::Ref().iDetector("MM01X1__");
    MM01z =  PaSetup::Ref().Detector(idet).Z();*/
  //idet = PaSetup::Ref().iDetector("MM02X1__");
  //MM02z =  PaSetup::Ref().Detector(idet).Z();
  //idet = PaSetup::Ref().iDetector("MM03X1__");
  //MM03z =  PaSetup::Ref().Detector(idet).Z();

  idet = PaSetup::Ref().iDetector("GM01U1__"); // first detector in zone 2
  Z2Az = PaSetup::Ref().Detector(idet).Z();

  //idet = PaSetup::Ref().iDetector("DR01X1__"); // last detector in zone 2
  // value taken from detectors.52959.plus.dat (2006) because RW is not present in 2004
  Z2Bz = 994.4913; //PaSetup::Ref().Detector(idet).Z();
  //check if MC analysis
  //use flag -T mc , or -T MC to indicate input files are mc
  //if input is data do not use "-T" option
  if(Phast::Ref().TextUserFlag(0)=="mc3"){ fIsMC=true;fMCtargetType=-3;}
  else if(Phast::Ref().TextUserFlag(0)=="MC3"){fIsMC=true;fMCtargetType=-3;}
  else if(Phast::Ref().TextUserFlag(0)=="mc5"){ fIsMC=true;fMCtargetType=-5;}
  else if(Phast::Ref().TextUserFlag(0)=="MC5"){fIsMC=true;fMCtargetType=-5;}
  else if(Phast::Ref().TextUserFlag(0)=="MC2012"){fIsMC=true;fMCtargetType=-12;}
  else if(Phast::Ref().TextUserFlag(0)=="MC2016"){fIsMC=true;fMCtargetType=-1;}
  else fIsMC=false;

  SetMC(fIsMC);


  cout << "LCAnalysis: done!" << endl;
}
LCAnalysis::~LCAnalysis()
{}

void LCAnalysis::SetMC(bool isMC)
{

  fIsMC=isMC;
  if(isMC){
    //cout<< "is mc"<<endl;
    switch(fAnalysisType){
    case 1:
      fDISEvtTree->Branch("DISMCEvt","DISEventMCData",&fDISMCEvt);
      fDISEvtTree->Branch("MCHadrons","std::vector<HadronMCData>",&fMCHadronsPtr);
    }

    //cout<<"LCAnalysis::SetMC: MC job! Tree branches created. "<<endl
    cout<<"     target type: "<<fMCtargetType<<endl;//" (good for "<<(fMCtargetType>-3 ? "2004":"2006")<<")"<<endl;
  }
  //cout<<"check SetMC"<<endl;
}



void LCAnalysis::ReadCutFile()
{
  char cutfilename[1024];
  if(Phast::Ref().NUserFlag()>1){
    sprintf(cutfilename,CUTFILENAMEVAR,Phast::Ref().UserFlag(1));
  }else{
    sprintf(cutfilename,CUTFILENAME);
  }

  cout<<"LCAnalysis::ReadCutFile: reading cut file "<<cutfilename<<" ..."<<endl;

  ifstream fin(cutfilename);
  string line,name;
  double minVal,maxVal;
  fCutVector.clear();

  while( getline(fin,line).rdstate() == ios::goodbit ){
    cout<<line<<endl;
    istringstream str(line);
    str >> name;
    if( name[0] == '#' ) continue;  // skip comments
    if( (str >> minVal).fail() ){   // "check cut"
      fCutVector.push_back(new LCEventSelectionCut(name,this));
    }else{ // "range cut"
      if( (str >> maxVal).fail() ) maxVal=-1.0; // no upper limit for range
      fCutVector.push_back(new LCEventSelectionRangeCut(name,this,minVal,maxVal));
    }
  }
  cout<<"check ReadCutFile"<<endl;
}



//------------------------------------------------ general interface
void LCAnalysis::DoEvent(PaEvent& ev)
{
  ++fTotNEvts; // increment total number of events
  if(ev.iBestPrimaryVertex()>=0)count_bestpv++;
   if( !fIsMC &&(ev.RunNum() >52564 && ev.RunNum() <54639))RescaleMom(ev,true);

  // call event handler
  (this->*fEvHandler)(ev);
}

void LCAnalysis::JobEnd()
{
  switch(fAnalysisType){
  case 0: PrintEventStats(); break;
  case 1: SaveHadronTree(); break;
  case 2: PrintEventStats(); break;
  default: PrintEventStats();
  }
}



//------------------------------------------------ event selection
void LCAnalysis::SelectEvent(PaEvent& ev)
{
  ResetEventInfo();       // put flags to false and pointers to 0
  fEvent=&ev;
  if(fDISEvtTree) fDISEvt->Reset();
  ReadMinimumEventInfo(); // set some event info member variables

  fEventAccepted = false;

  //--- loop over cuts. If one check fails, it exits without saving
  vector<LCEventSelectionCut*>::iterator it;
  for(it=fCutVector.begin(); it!=fCutVector.end(); ++it){
    if( !(*it)->Check()
	//&& (*it)->GetName() != "intInTarget"
	) return;
  }

  //--- Check Trigger
  // CheckTriggerMask(ev);
  // if( fTrigOk  ) ++fTrigger;
  // retaining events with all trigger masks for now. This cut should
  // be incorporated into the dynamic cuts.
  //return;

  //--- Save Event
  fEventAccepted = true;
  fDISEvtList->Enter(Phast::Ref().iev_en);
  if(fSaveEvents) ev.TagToSave();
  //--- Discard some raw info (taken from Yann, hopefully he is right)
  //ev.Discard(0x4);
  //cout<<"check Select Event"<<endl;
}




//---------------------------------------------------- utility methods
void LCAnalysis::ResetEventInfo()
{
  fEvent=0;
  fBPV=0;
  fiBPV=-1;
  fimu0=fimu1=-1;
  fReconsEvent=false;
  fValidMu=false;
  //cout<<"check ResetEventInfo"<<endl;
}

void LCAnalysis::ReadMinimumEventInfo()
{
  fiBPV = fEvent->iBestPrimaryVertex();
  if(fiBPV==-1) return;

  fBPV = &(fEvent->vVertex(fiBPV));

  //--- Check if there is a mu'
  if(!fBPV->IsPrimary()){
    printf("Warning! Best primary vertex is not primary!\n");
    return;
  }
  fimu0 = fBPV->InParticle();
  //fimu1 = fBPV->iMuPrim(false,false,true,false);
  fimu1 = fBPV->iMuPrim();
  if(fimu0 == -1) return;
  if(fimu1 == -1) return;

  fReconsEvent = true;

  //--- Get trigger mask
  fTrigMask = fEvent->TrigMask();

  //--- Get Muon four vectors
  SetMuKinematics(*fEvent,fiBPV,fimu0,fimu1);


  //--- store info of DISevent
  if(fDISEvtTree &&( fIsMC || fReconsEvent)){
    CopyDISEvtData(fReconsEvent);
    fDISEvtTree->Fill();
  }

  //--- Fill hist. info
  if(fHistInfo) fHistInfo->Fill(fQ2,fXX0mu1,fBPV->Z(),fTrigMask);
  //cout<<"check ReadMinimumEventInfo"<<endl;
}


void LCAnalysis::CopyDISEvtData(int pReconsEvent)
{
  fDISEvt->runNo      = fEvent->RunNum();
  fDISEvt->spillNo    = fEvent->SpillNum();
  fDISEvt->evtInSpill = fEvent->EvInSpill();
  fDISEvt->trigMask   = fEvent->TrigMask();
  fDISEvt->evNo       = fEvent->UniqueEvNum();
  fDISEvt->timeInSpill= fEvent->TimeInSpill();
  if(!pReconsEvent) return;

  fDISEvt->x = fBPV->X();
  fDISEvt->y = fBPV->Y();
  fDISEvt->z = fBPV->Z();
  fDISEvt->p0x = fkMu0.X();
  fDISEvt->p0y = fkMu0.Y();
  fDISEvt->p0z = fkMu0.Z();
  fDISEvt->p1x = fkMu1.X();
  fDISEvt->p1y = fkMu1.Y();
  fDISEvt->p1z = fkMu1.Z();
  fDISEvt->E_beam = fkMu0.E();
  fDISEvt->E_mu_prim = fkMu1.E();
  fDISEvt->XX0 = fXX0mu1;
  fDISEvt->cellsCrossed = fCellsCrossed;
  fDISEvt->HM04x = HM04x;
  fDISEvt->HM04y = HM04y;
  fDISEvt->HM05x = HM05x;
  fDISEvt->HM05y = HM05y;
  fDISEvt->HO03x = HO03x;
  fDISEvt->HO03y = HO03y;
  fDISEvt->HO04x = HO04x;
  fDISEvt->HO04y = HO04y;
  fDISEvt->backPropFlag = fChi2CutFlag;
  // cout<<"Reconstructed Event : "<<pReconsEvent<<endl;
  // cout<<"check CopyDISEvtData"<<endl;
}



bool LCAnalysis::IsThereABestPV()
{
  //--- Check if there is a best primary vertex
  return fiBPV != -1;
  //cout<<"check IsThereABestPV"<<endl;
}
bool LCAnalysis::IsMu1Reconstructed()
{
  //--- Check if there is a reconstructed mu'
  return fimu1 != -1;
  //cout<<"check IsMu1Reconstructed"<<endl;
}
bool LCAnalysis::IsMu0Valid(PaTPar par)
{
  //--- Check if mu is valid from cov matrix
  // cout << par(5,5) << endl;
  // cout << ((par(5,5)<20e-9) ? par(5,5) : 0) << endl;
  return ((par(5,5)<2e-9) ? 1 : 0);
  //cout<<"check IsMu0Valid"<<endl;
}

bool LCAnalysis::InteractionInTarget()
{
  const PaTPar& ParamMu0 = fEvent->vParticle(fimu0).ParInVtx(fiBPV);
  int runno = fIsMC ? fMCtargetType : fEvent->RunNum();
  return
    PaAlgo::InTarget(ParamMu0,'U',runno) ||
    PaAlgo::InTarget(ParamMu0,'C',runno) ||
    PaAlgo::InTarget(ParamMu0,'D',runno);
  //cout<<"check InteractionInTarget"<<endl;
}

bool LCAnalysis::InteractionInTarget2009()
{
  // Check if the interaction is in the lH2 volume:
  // 1. Check longitudinal pos.  -68.5 cm < z < -28.5 cm (Nicole d'Hose)
  // 2. Transverse position within a circle of the radius of 1.6 cm
  // (see talk of M. Boer, analysis meeting of 8.7.2010)
  double x = fBPV->X();
  double y = fBPV->Y();
  double z = fBPV->Z();

  // longitudinal cut
  if( z < -68.5 || z > -28.5 ) return false;

  // radial cut (1.6^2 = 2.56)
  if( x*x + y*y > 2.56 ) return false;

  return true;
}

bool LCAnalysis::InteractionInTarget2016()
{
  const PaTPar& ParamMu0 = fEvent->vParticle(fimu0).ParInVtx(fiBPV);
  int runno = fIsMC ? fMCtargetType : fEvent->RunNum();

  return fTcell->TargetCell::InTarget(ParamMu0,runno);
}


bool LCAnalysis::CellsCrossed()
{
  const PaTPar& ParamMu0 = fEvent->vParticle(fimu0).ParInVtx(fiBPV);
  return PaAlgo::CrossCells(ParamMu0,fEvent->RunNum());
  //cout<<"check CellsCrossed"<<endl;
}

bool LCAnalysis::CellsCrossed2016()
{
  const PaTPar& ParamMu0 = fEvent->vParticle(fimu0).ParInVtx(fiBPV);
  return fTcell->TargetCell::CrossCells(ParamMu0,fEvent->RunNum());
  //cout<<"check CellsCrossed"<<endl;
}


void LCAnalysis::CheckTriggerMask(const PaEvent& ev)
{
  int trigMask = ev.TrigMask();
  //cout << trigMask << endl;
  fTrigMask = trigMask;
  //fTrigOk = trigMask&8 || trigMask&256 ;  // true if OT || incl. MT
  fTrigOk = trigMask&1 || trigMask&2 || trigMask&3 || trigMask&9; // 2016 triggers
  // cout<<"check CheckTriggerMask"<<endl;
}

void LCAnalysis::SetMuKinematics(const PaEvent& ev,const int& iVtx,
				 const int& imu0,const int& imu1)
{
  //cout<<"First check SetMuKinematics"<<endl;
  const PaParticle& Mu0   = ev.vParticle(imu0); // the beam muon
  const PaParticle& Mu1   = ev.vParticle(imu1); // the scattered muon
  const PaTPar& ParamMu0  = Mu0.ParInVtx(iVtx); // fitted mu  parameters in the primary vertex (TODO : element (5,5) check <20*10^-9)
  const PaTPar& ParamMu1  = Mu1.ParInVtx(iVtx); // fitted mu' parameters in the primary vertex
  fkMu0 = ParamMu0.LzVec(M_mu); // beam      mu Lorentz vector
  fkMu1 = ParamMu1.LzVec(M_mu); // scattered mu Lorentz vector

  fEbeam = fkMu0.E();
  fQ2    = PaAlgo::Q2(fkMu0, fkMu1);
  fnu    = (fEbeam - fkMu1.E());
  fy     = fnu/fEbeam;
  fxBj   = PaAlgo::xbj(fkMu0,fkMu1);
  fW2    = PaAlgo::W2(fkMu0,fkMu1);

  int itr = Mu1.iTrack();
  const PaTrack& track = ev.vTrack(itr);
  fXX0mu1 = track.XX0();

  fzVTX = fBPV->Z();

  // save mu1 pos at MT hodos
  int Npars = track.NTPar();
  PaTPar partr = track.vTPar(Npars-1); //
  PaTPar parHM;
  partr.Extrapolate(HM04z, parHM, false);
  HM04x = parHM(1);
  HM04y = parHM(2);
  partr.Extrapolate(HM05z, parHM, false);
  HM05x = parHM(1);
  HM05y = parHM(2);
  partr.Extrapolate(HO03z, parHM, false);
  HO03x = parHM(1);
  HO03y = parHM(2);
  partr.Extrapolate(HO04z, parHM, false);
  HO04x = parHM(1);
  HO04y = parHM(2);

  // check if all cells crossed
  //fCellsCrossed = PaAlgo::CrossCells(ParamMu0,fIsMC? fMCtargetType : ev.RunNum());
  if((ev.RunNum() >52564 && ev.RunNum() <54639)||fMCtargetType==-5){ //2006 data and MC
    fCellsCrossed = PaAlgo::CrossCells(ParamMu0,fMCtargetType);
    if(fCellsCrossed)fCellsCrossed = PaAlgo::CrossCells(ParamMu0,fdatatargetType);
  }
  else if( (ev.RunNum() > 269918)||fMCtargetType==-1){ //2016 data and MC
    fCellsCrossed = fTcell->TargetCell::CrossCells(ParamMu0,fMCtargetType);
    if(fCellsCrossed)fCellsCrossed = fTcell->TargetCell::CrossCells(ParamMu0,ev.RunNum());
  }
  else {
      cout << "Year not found.. " << endl;
  }

  // save "back propagation flag"
  if((ev.RunNum() >52564 && ev.RunNum() <54639)){  fChi2CutFlag = Mu0.Chi2CutFlag();}
  //cout<<"check SetMuKinematics "<<endl;

  //save 2012 "back prop flag"
  if( ev.RunNum() >107923 && ev.RunNum() <109082){
  const PaTrack& Mu0track   = ev.vTrack(imu0); // the beam muon track reference
  fChi2CutFlag = (Mu0track.NHitsFoundInDetect("BM")>3)?(true):(false);}

  //save 2016 "back prop flag"
  if((269918<ev.RunNum())){
  const PaTrack& Mu0track   = ev.vTrack(imu0); // the beam muon track reference
  fChi2CutFlag = (Mu0track.NHitsFoundInDetect("BM")>3)?(true):(false);}
}

double LCAnalysis::GetMassPid(int pid) const
{
  switch(pid){
  case 0: return M_pi;
  case 1: return M_K;
  case 2: return M_p;
  case 4: return M_mu;
  default: return M_pi;
  }
  //cout<<"check GetMassPid"<<endl;
}

void LCAnalysis::Calcz(const PaTrack& tr) // calculate a preliminary z
{
  fLikePid = fPid->LikePid(tr);
  double mass = GetMassPid(fLikePid);
  const PaTPar& param = tr.vTPar(0);
  TLorentzVector Ph = param.LzVec(mass);
  fz = Ph.E()/fnu;
  fMom = param.Mom();
  fthP = param.Theta();
  fphP = param.Phi();
  //cout<<"check Calcz"<<endl;
}


//------------------------------------------------ hadron selection
void LCAnalysis::FindHadrons(PaEvent& ev)
{
  //cout<<"First check FindHadrons"<<endl;
  fEvent = &ev;

  fDISEvt->Reset();
  fHadrons.clear();


  //--- get primary vertex and muon kinematics
  fiBPV = ev.iBestPrimaryVertex();
  int imu0=-1, imu1=-1;
  // cout << "fiBPV : " << fiBPV;

  if( fiBPV >= 0 ){
    fBPV = &(ev.vVertex(fiBPV));
    const PaVertex& v = ev.vVertex(fiBPV);

    imu0 = fimu0 = v.InParticle();
    //imu1 = v.iMuPrim(false,false,true,false);
    imu1 = fimu1 = v.iMuPrim();

    if(imu0)
    {
      const PaParticle& Mu0   = ev.vParticle(imu0);
      const PaTPar& ParamMu0  = Mu0.ParInVtx(fiBPV);

      fValidMu=IsMu0Valid(ParamMu0);
    }

    // cout << ", Vx : " << v.X() << ", Vy : " << v.Y() << ", Vz : " << v.Z() << ", mu' rec. :" << v.iMuPrim();
  }
  fReconsEvent = IsThereABestPV() && IsMu1Reconstructed();
  if(fReconsEvent)count_mup++;

  // cout << ", Recons. : " << fReconsEvent << endl;

  if( fIsMC ){ // read MC event info

    fDISMCEvt->Reset();
    ///LUND block data
    NLUDATA ld;                       // create NLUDATA structure
    double mcWeight = 1.00;
    vector<LUJET> lujets;         // create LUJET structure
    int nparticles;               // number of particles
    ev.MCgen(ld);                // get LUDATA
    ev.MCgen(nparticles,lujets);
     ///get mc weight
    if (ld.uservar[2] != 0)mcWeight = ld.uservar[2]; //HEPPGEN Weights
    fDISMCEvt->mcWeight = mcWeight;
    ///find dd events
    int difType=0;
    for(int i=0;i<nparticles;i++){ //newwww
      if(lujets[i].k[1]==2210 || lujets[i].k[1]==2110){  difType= 1; break;} } fDISMCEvt->difType = difType;
    //get type of event from RADGEN
    // irad=0 - no radiation
    // irad=1 - inelastic (1<x<0)
    // irad=2 - quasi-elastic on p or on nucleon (x=1)
    // irad=3 - elastic (x=M_A/M)
    // irad=4 is irad=1 but W^2<4 (i.e. not useful to call jetset for hadronization, SET BY Yann)
    int irad = int(ld.uservar[19]+.5); fDISMCEvt->irad = irad;
      fDISMCEvt->MC_nuTr = ld.u; fDISMCEvt->MC_Q2Tr = ld.q2;
      fDISMCEvt->MC_xTr = ld.x; fDISMCEvt->MC_yTr = ld.y;
    float w1 = ld.uservar[18], w2 =  ld.uservar[17];
    if (w1!=0 || w2!=0) {
      fDISMCEvt->w1 = w1; fDISMCEvt->w2 = w2;
    }
    //...else we not dealing w/ a RADGEN event: weight left =1
    //end LUND block data

    const PaMCvertex& mcVtx = ev.vMCvertex(0);
    if( !mcVtx.IsPrimary() ){
      cout<<"LCAnalysis::FindHadrons: MC vertex is not primary"<<endl;
    }

    if( mcVtx.iBeam()<0 || mcVtx.iBeam()>=mcVtx.NMCtrack() ){
      cout<<"mcVtx.iBeam() = "<<mcVtx.iBeam()<<"; skipping event"<<endl;
      mcVtx.Print();
      return;
    }

    const PaMCtrack& mu0 = fEvent->vMCtrack(mcVtx.iBeam());
    const PaMCtrack& mu1 = fEvent->vMCtrack(1);
    //if( mu1.Pid() != mu0.Pid() )
    // cout<<"LCAnalysis::FindHadronsMC: mu1 is not a muon!"<< mu1.Pid()<<" and "<< mu0.Pid()<<endl;
    if( mu1.Pid() == 5 )count_mupp++;
    if( mu1.Pid() == 6 )count_munp++;
    if( mu0.Pid() == 5 )count_mup0++;
    if( mu0.Pid() == 6 )count_mun0++;
    fDISMCEvt->MC_vx = mcVtx.Pos(0); fDISMCEvt->MC_vy = mcVtx.Pos(1); fDISMCEvt->MC_vz = mcVtx.Pos(2);
    TLorentzVector kmu0 = mu0.LzVec();
    TLorentzVector kmu1 = mu1.LzVec();
    ////////addition for different radgen/lepto
    //TLorentzVector gamm = kmu0; gamm -= kmu1;
    //const vector<PaMCgen> &vMCgen = ev.vMCgen();
    //static TLorentzVector gammTr; if (!vMCgen.empty()) {
    //  if (!parseMCgen(vMCgen,gammTr)) {
//	printf("** LCA:\a Evt #%f Missing LUDATA or LUJET\n",
//	     fDISEvt->evNo); exit(1);
//      }
//    }
    //if (fDISMCEvt->irad && mcVtx.NMCtrack()==4) {// Elastic
    //  fDISMCEvt->irad = 2;
    // }
    ////////////end addition

    fDISMCEvt->MC_p0x = -kmu0.X();
    fDISMCEvt->MC_p0y = -kmu0.Y();
    fDISMCEvt->MC_p0z = -kmu0.Z();

    fDISMCEvt->MC_p1x =  kmu1.X();
    fDISMCEvt->MC_p1y =  kmu1.Y();
    fDISMCEvt->MC_p1z =  kmu1.Z();

    fDISMCEvt->recons = fReconsEvent;

    int Ntr = mcVtx.NMCtrack();
    int Nmuons=0;
    fMCHadrons.clear();
    for(int ih=0; ih<Ntr; ih++){// loop on secondaries

      const PaMCtrack& hadr = fEvent->vMCtrack(mcVtx.iMCtrack(ih));

      if( hadr.Pid() == 5 || hadr.Pid() == 6 ){ // skip muons
	++Nmuons;
	continue;
	}

      HadronMCData mcHadron;

      TLorentzVector Ph = hadr.LzVec();
      mcHadron.P   = Ph.Rho();
      mcHadron.th  = Ph.Theta();
      mcHadron.ph  = Ph.Phi();
      mcHadron.charge = hadr.Q();
      mcHadron.pid = hadr.Pid();

      mcHadron.recons = false;
      std::set<int>::iterator it=hadr.sTrkRef().begin();
      for(; it!=hadr.sTrkRef().end(); ++it){
	if( ev.vTrack(*it).HasMom() ){
	  if(mcHadron.recons){
	    //cout<<"LCAnalysis::FindHadrons: More than one reconstructed track associated with MC track!"<<endl;
	    ++fNMChadrMultAssoc;
	  }
	  mcHadron.recons = true;
	  mcHadron.itrack = *it;
	}
      }

      ++fNMChadr;
      fMCHadrons.push_back(mcHadron);

    }// end loop on secondaries

    if( Nmuons != 2 ){
      cout<<"LCAnalysis::FindHadronsMC: not two muons in vertex ("<<Nmuons<<") !"<<endl;
      //fEvent->vMCgen()[0].Print(3);
    }

  } // end of MC part

  //cout<<"3rd check FindHadron"<<endl;
  if( fReconsEvent ){ // continue only if the event is reconstructed

    SetMuKinematics(ev,fiBPV,imu0,imu1);
          const PaVertex& v = ev.vVertex(fiBPV);


	  /*
    //--- loop on events photons added by erin
	  int NoutPhot = ev.NCaloClus();
	for(int i=0; i<NoutPhot; ++i){
	const PaCaloClus& outPhot = ev.vCaloClus(i);
	const string& nam = outPhot.CalorimName();
	if(nam.find("EC02P1") != 0) continue; //reject ecal2 photons
	TVector3 p_vec;
	p_vec.SetXYZ(0.,0.,0.);
	p_vec.SetXYZ(outPhot.X()-v.X(),outPhot.Y()-v.Y(),outPhot.Z()-v.Z());
      HadronData hadron;
      hadron.P  = outPhot.E();
      hadron.th = p_vec.Theta();
      hadron.ph = p_vec.Phi();
      hadron.charge = 0;

      //--- fill hadron vector
      if(hadron.P>0){
	fHadrons.push_back(hadron);
      }
	}//end photon loop


	  */


    //--- loop on primary vertex outgoing particles
    int NoutPart = v.NOutParticles();
    fTotNHadr += NoutPart;
    int ip=0,itr=0,Nh=0;
    for(int i=0; i<NoutPart; ++i){
      ip = v.iOutParticle(i);

      //--- check if it is mu'
      if(ip == imu1) continue;
      ++fNsec;

      const PaParticle& outPart = ev.vParticle(ip);
      if(outPart.IsBeam()) // this check should be always negative at this stage!
	printf("LCAnalysis::FindHadrons: outgoing beam particle!\n");

      itr = outPart.iTrack();
      if(itr == -1)
	printf("LCAnalysis::FindHadrons: outgoing particle without track!\n");

      HadronData hadron;
      const PaTPar& param = outPart.ParInVtx(fiBPV);
      hadron.P  = param.Mom();
      hadron.th = param.Theta();
      hadron.ph = param.Phi();


      //--- check z > 1
      const PaTrack& track = ev.vTrack(itr);
      Calcz(track);  // also fLikePid is set!
      if(fz < 1) ++fNzlt1;

      //--- check Zfirst
      fZfirst = track.ZFirst();
      if(fZfirst > fZref) continue;
      ++fNZfirst;

      //--- check Zlast
      fZlast = track.ZLast();
      if(fZlast < fZref) continue;
      ++fNZlast;

      //--- check X/X0
      hadron.XX0 = track.XX0();
      ++fNXX0;

      //--- check trigger
      // CheckTriggerMask(ev);
      // if(fTrigOk) {++fNTrigHad;++fTrigger;}

      // cout << "check end cut" << endl;

      //Inelasticity variable
      //if
      //TLorentzVector p(0,0,0,M_p)                // proton's 4-vector (in rest)
      //TLorentzVector q = fkMu0-fkMu1;            // virtual photon's 4-vector
      //TLorentzVector LzVc_h =param.LzVec(M_pi);  //rho vector (CM vector)

      //--- polar angle at RICH entrance
      PaTPar tParRich;
      //track.vTPar(0).Extrapolate(Zrich, tParRich);
      track.Extrapolate(Zrich, tParRich);
      fthRICH = hadron.thRICH = tParRich.Theta(false);
      hadron.RICHx=tParRich.Pos(0);                   //extrapolation RICH added by Quiela
      hadron.RICHy=tParRich.Pos(1);                   //extrapolation RICH added by Quiel

      //--- get likelihoods
      for(int i=0; i<6; ++i) hadron.LH[i] = fPid->GetLike(i,track);

      //--- check if track falls into hadronic calorimeter acceptance
      PaTPar tParHCAL;
      double xc,yc;
      track.vTPar(0).Extrapolate(ZHCAL1, tParHCAL, false);  // HCAL1
      hadron.inHCALacc = fHadrCal1->iCell(tParHCAL(1),tParHCAL(2),xc,yc) != -1;

      if( !hadron.inHCALacc ){
	track.vTPar(0).Extrapolate(ZHCAL2, tParHCAL, false);  // HCAL2
	hadron.inHCALacc = fHadrCal2->iCell(tParHCAL(1),tParHCAL(2),xc,yc) != -1;
      }

      //--- get hadronic calorimeter signals
      int iC;
      for(int i=0; i<outPart.NCalorim(); ++i){
	iC = outPart.iCalorim(i);
	const PaCaloClus & clus = ev.vCaloClus(iC);
	if( clus.CalorimName()[0] == 'H' ) // only HCALs
	  hadron.HCAL += clus.E();
      }


      //--- calculate track extrapolations
      PaTPar tParH;
      /*track.vTPar(0).Extrapolate(MM01z,tParH,false);
	hadron.MM01x = tParH(1); hadron.MM01y = tParH(2);*/
      track.vTPar(0).Extrapolate(MM02z,tParH,false);
      hadron.MM02x = tParH(1); hadron.MM02y = tParH(2);
      track.vTPar(0).Extrapolate(MM03z,tParH,false);
      hadron.MM03x = tParH(1); hadron.MM03y = tParH(2);

      track.vTPar(0).Extrapolate(Z2Az,tParH,false);
      hadron.Z2Ax = tParH(1); hadron.Z2Ay = tParH(2);
      track.vTPar(0).Extrapolate(Z2Bz,tParH,false);
      hadron.Z2Bx = tParH(1); hadron.Z2By = tParH(2);



      if( fIsMC ){// special treatment for MC

	//--- set vector position of this hadron to MC branch
	vector<HadronMCData>::iterator mch = fMCHadrons.begin();
	for( ; mch != fMCHadrons.end(); ++mch){
	  if( mch->itrack == itr ) mch->recHadIdx = Nh;
	}

	//--- find the true pid of this track if any
	int imctr = track.iMCtrack();
	if( imctr >= 0 ){
	  const PaMCtrack& mcTrack = ev.vMCtrack(imctr);
	  hadron.MCpid = mcTrack.Pid();
	}
      } //end special treatment for MC

      //--- fill hadron vector
      hadron.charge = outPart.Q();
      fHadrons.push_back(hadron);
      Nh = fHadrons.size();


    }//--- end loop on primary vertex secondaries

  // cout<< fReconsEvent <<endl;
  }// end if event reconstructed
  if(fIsMC || (fReconsEvent && fValidMu)){
  CopyDISEvtData(fReconsEvent);
  fDISEvtTree->Fill();
  // cout<<"check Event saved"<<endl;
  }
  // cout<<"check FindHadrons"<<endl;
} //end FindHadrons
// ***************************************************************************
// ******************************   parseMCgen  ******************************
// ***************************************************************************
bool LCAnalysis::parseMCgen(const vector<PaMCgen> &vMCgen,
			    TLorentzVector &gammTr)
// Method
// - Gets input from PaEvent::vMCgen
// - Sets several DISEventMCData members (which are otherwise init'd by ctor)
//   - irad = flag set !=0 if hard RC photon.
//   - True nu and Q2 (expected to differ (beyond rounding approx.) from
//    standard generated nu and Q2 only if hard RC photon emitted).
//   - Radgen weight (init'd = 1)
// - Fills arg. gammTr w/ hard photon 4-vector if hard RC photon.
// - Requires both LUDATA/LUJET (type code 110/102). Returns false if not.
{ //tt
  int igen; unsigned int match;
  for (igen = 0, match = 0; igen<(int)vMCgen.size(); igen++) {
    const PaMCgen &gen = vMCgen[igen];

    if (gen.vInt().back()==110) {
      //const NLUDATA &ld = gen.NLundData();
      //double mcWeight = ld.uservar[2]; fDISMCEvt->mcWeight = mcWeight; //added by erin to get event weight

      //int difType = ld.lst[23]; fDISMCEvt->difType = difType; //added by erin values 1 for dd and 0 for elastic (see documentation of HEPGEN) diffractive evnt

      //int irad = int(ld.uservar[19]+.5); fDISMCEvt->irad = irad; if (irad) {
      //	fDISMCEvt->MC_nuTr = ld.u; fDISMCEvt->MC_Q2Tr = ld.q2;
      //	fDISMCEvt->MC_w = ld.uservar[18];
      // }
      match |= 0x1; if (match&0x2) break;
    }
    else if (gen.vInt().back()==102) { // LUJET
      const LUJET &lj = vMCgen[igen].LundJet();
      if (abs(lj.k[1]==22) && !(match&0x2)) {
	gammTr.SetPxPyPzE(lj.p[0],lj.p[1],lj.p[2],lj.p[3]);
	match |= 0x2; if (match&0x1) break;
      }
    }
  }


  //#define LCA_DEBUG_RAD 2
#ifdef LCA_DEBUG_RAD
  fDISMCEvt->CalcKin();
  if (fDISMCEvt->irad) {
# if LCA_DEBUG_RAD > 1
    printf("Evt %d x,y,q2,nu %.4f,%.3f,%.4f,%.2f -> %.4f,%.3f,%.4f,%.2f\n",
	   (int)ev.UniqueEvNum(),
	   fDISMCEvt->MC_xBj, fDISMCEvt->MC_Y,  fDISMCEvt->MC_Q2,  fDISMCEvt->MC_nu,
	   fDISMCEvt->MC_xTr, fDISMCEvt->MC_yTr,fDISMCEvt->MC_Q2Tr,fDISMCEvt->MC_nuTr);
# endif
  }
  float Q2XC = -gammTr.M2(), Q2Orig = fDISMCEvt->irad?
    fDISMCEvt->MC_Q2Tr:fDISMCEvt->MC_Q2, diff = (Q2XC-Q2Orig)/Q2Orig;
  if (fabs(diff)>1e-2 || LCA_DEBUG_RAD>1)
    printf("Evt %d: XC %.4f -> %4f  %f\n",
	   (int)ev.UniqueEvNum(),Q2Orig,Q2XC,diff);
#endif
  return match==0x3;
}
// ***************************************************************************
// ****************************** FindHadronsMC ******************************
// ***************************************************************************
int LCAnalysis::FindHadronsMC(const PaMCvertex& vtx)
{

  //----------------- fill MC inclusive variables
  if( vtx.iBeam()<0 || vtx.iBeam()>=vtx.NMCtrack() ){
    cout<<"vtx.iBeam() = "<<vtx.iBeam()<<"; skipping event"<<endl;
    vtx.Print();
    return -1;
  }

  const PaMCtrack& mu0 = fEvent->vMCtrack(vtx.iBeam());
  const PaMCtrack& mu1 = fEvent->vMCtrack(1);
  if( mu1.Pid() != mu0.Pid() )
    cout<<"LCAnalysis::FindHadronsMC: mu1 is not a muon!"<< mu1.Pid()<<" and "<< mu0.Pid()<<endl;

  TLorentzVector kmu0 = mu0.LzVec();
  TLorentzVector kmu1 = mu1.LzVec();

  fkMu0 = TLorentzVector(-kmu0.X(),-kmu0.Y(),-kmu0.Z(),kmu0.T());
  fkMu1 = TLorentzVector(kmu1.X(),kmu1.Y(),kmu1.Z(),kmu1.T());

  fEbeam = kmu0.E();
  fEmu1  = kmu1.E();
  fThmu1 = kmu1.Theta();
  fPhmu1 = kmu1.Phi();
  fQ2    = PaAlgo::Q2(fkMu0, fkMu1);
  fnu    = (fEbeam - fEmu1);
  fy     = fnu/fEbeam;
  fxBj   = PaAlgo::xbj(fkMu0,fkMu1);
  fW2    = PaAlgo::W2(fkMu0,fkMu1);
  fzVTX  = vtx.Pos(2);


  fEvtTreeMC->Fill();

  //---------- find MC hadrons
  int Ntr = vtx.NMCtrack();
  int Nmuons=0;
  for(int ih=0; ih<Ntr; ih++){// loop on secondaries

    const PaMCtrack& hadr = fEvent->vMCtrack(vtx.iMCtrack(ih));
    if( hadr.Pid() == mu0.Pid() ){ // skip muons
      ++Nmuons;
      //cout<<"    "<<ih<<"  "<<vtx.iMCtrack(ih)<<endl;
      continue;
    }

    TLorentzVector Ph = hadr.LzVec();
    fMCPid = hadr.Pid();
    fCharge = hadr.Q();
    fz = Ph.E()/fnu;
    fMom = Ph.Rho();
    fthP = Ph.Theta();
    fphP = Ph.Phi();

    fHadrTreeMC->Fill();
  }// end loop on secondaries

  if( Nmuons != 2 ){
    cout<<"LCAnalysis::FindHadronsMC: not two muons in vertex ("<<Nmuons<<") !"<<endl;
    //fEvent->vMCgen()[0].Print(3);
  }

  return 1;
  //cout<<"check FindHadronsMC"<<endl;
}




void LCAnalysis::InitHadronTree()
{
  fOutFile  = new TFile(Phast::Ref().out_file_name.data(),
			"RECREATE");
  //cout<<"fOutFile open? \t"<<fOutFile<<endl;
  fEvtTree = new TTree("EvtTree","EvtTree");
  fEvtTree->Branch("xBj",&fxBj);
  fEvtTree->Branch("Q2",&fQ2);
  fEvtTree->Branch("W2",&fW2);
  fEvtTree->Branch("y",&fy);
  fEvtTree->Branch("zVTX",&fzVTX);
  fEvtTree->Branch("TrigOk",&fTrigOk);
  fEvtTree->Branch("TrigMask",&fTrigMask);
  //-- mu1 pos at MT hodos
  fEvtTree->Branch("HM04x",&HM04x);
  fEvtTree->Branch("HM04y",&HM04y);
  fEvtTree->Branch("HM05x",&HM05x);
  fEvtTree->Branch("HM05y",&HM05y);
  int idet = PaSetup::Ref().iDetector("HM04X1_d");
  const PaDetect& HM04 =  PaSetup::Ref().Detector(idet);
  HM04z = HM04.Z();
  idet = PaSetup::Ref().iDetector("HM05X1_d");
  const PaDetect& HM05 =  PaSetup::Ref().Detector(idet);
  HM05z = HM05.Z();


  fHadrTree = new TTree("HadrTree","analysisTree");
  fHadrTree->Branch("TrigMask",&fTrigMask);
  fHadrTree->Branch("TrigOk",&fTrigOk);
  fHadrTree->Branch("xBj",&fxBj);
  fHadrTree->Branch("Q2",&fQ2);
  fHadrTree->Branch("y",&fy);
  fHadrTree->Branch("W2",&fW2);
  //-- mu1 pos at MT hodos
  fHadrTree->Branch("HM04x",&HM04x);
  fHadrTree->Branch("HM04y",&HM04y);
  fHadrTree->Branch("HM05x",&HM05x);
  fHadrTree->Branch("HM05y",&HM05y);
  //-----------------------------
  fHadrTree->Branch("z",&fz);
  fHadrTree->Branch("Q",&fCharge);
  fHadrTree->Branch("P",&fMom);
  fHadrTree->Branch("thP",&fthP);
  fHadrTree->Branch("phP",&fphP);
  fHadrTree->Branch("thRICH",&fthRICH);
  fHadrTree->Branch("RICHx",&RICHx);
  fHadrTree->Branch("RICHy",&RICHy);
  fHadrTree->Branch("PID",&fLikePid);
  if(fIsMC) fHadrTree->Branch("truePID",&fTruePidMC);
  fHadrTree->Branch("LHpi",&fLpi);
  fHadrTree->Branch("LHK", &fLK );
  fHadrTree->Branch("LHp", &fLp );
  fHadrTree->Branch("LHe", &fLe );
  fHadrTree->Branch("LHmu",&fLmu );
  fHadrTree->Branch("LHBg",&fLBg);

  if(fIsMC){
    fEvtTreeMC = new TTree("EvtTreeMC","EvtTreeMC");
    fEvtTreeMC->Branch("xBj",&fxBj);
    fEvtTreeMC->Branch("Q2",&fQ2);
    fEvtTreeMC->Branch("W2",&fW2);
    fEvtTreeMC->Branch("y",&fy);
    fEvtTreeMC->Branch("E0",&fEbeam);
    fEvtTreeMC->Branch("E1",&fEmu1);
    fEvtTreeMC->Branch("th1",&fThmu1);
    fEvtTreeMC->Branch("ph1",&fPhmu1);
    fEvtTreeMC->Branch("zVTX",&fzVTX);
    fEvtTreeMC->Branch("TrigMask",&fTrigMask);
    fEvtTreeMC->Branch("recon",&fReconsEvent);

    fHadrTreeMC = new TTree("HadrTreeMC","analysisTreeMC");
    fHadrTreeMC->Branch("TrigMask",&fTrigMask);
    fHadrTreeMC->Branch("xBj",&fxBj);
    fHadrTreeMC->Branch("Q2",&fQ2);
    fHadrTreeMC->Branch("y",&fy);
    fHadrTreeMC->Branch("W2",&fW2);
    fHadrTreeMC->Branch("recon",&fReconsEvent);
    fHadrTreeMC->Branch("z",&fz);
    fHadrTreeMC->Branch("Q",&fCharge);
    fHadrTreeMC->Branch("P",&fMom);
    fHadrTreeMC->Branch("thP",&fthP);
    fHadrTreeMC->Branch("phP",&fphP);
    fHadrTreeMC->Branch("MCPid",&fMCPid);
  }
//cout<<"check InitHadronTree "<<endl;
}


short gNvtx;
float gNSecVtx;
TH1D *ghNvtx;
TH2D *ghSecVtx;

void LCAnalysis::SelectHadronPairs(PaEvent& ev)
{
  static PaPid *PID;

  if(fTotNEvts == 1){
    InitHadronPairTree();
    PID = new PaPid;
  }


  SelectEvent(ev);
  if(!fEventAccepted) return;

  //  cout<<"Event "<<fTotNEvts<<endl;

  int Nvtx = ev.NVertex();
  int NoutPart;
  ghNvtx->Fill(Nvtx);


  gNSecVtx=0.0;
  for(int iv=0; iv<Nvtx; ++iv){//--- loop on vertices
    const PaVertex& vtx = ev.vVertex(iv);
    bool primVtx = iv == fiBPV;
    //cout << iv <<"   "<< vtx.NOutParticles() <<endl;
    NoutPart = vtx.NOutParticles();
    gNSecVtx =+ NoutPart;
    ghSecVtx->Fill(Nvtx,NoutPart);


    if( primVtx && NoutPart < 3 ) continue; // case of primary vertex

    if( NoutPart > 1 ){ // only verteces with more than one secondary

      int ih1=-1,ih2=-1;
      int itr1=-1,itr2=-1;
      for(int ip=0; ip<NoutPart; ++ip){ // loop on paricles, finding pair candidates
	int ipart=vtx.iOutParticle(ip);
	int itrack=ev.vParticle(ipart).iTrack();
	const PaTrack& track = ev.vTrack(itrack);
	const PaTPar&  tpar  = ev.vParticle(ipart).ParInVtx(iv);
	if( itrack < 0                            ||  // skip if there is no track
	    !ev.vTrack(itrack).HasMom()           ||  // skip if no reconstr. mom.
	    track.CanBeMuon()                     ||  // skip muons
	    tpar.Mom() < 2. || tpar.Mom() > 70.       // range in momentum
	    ) continue;

	if( NoutPart > 3 ) continue; // particular case to be treated differently

	if(ih1<0) { ih1 = ipart; itr1 = itrack; }
	else      { ih2 = ipart; itr2 = itrack; }
      }// end loop on particles

      if( ih1<0 || ih2<0 ) continue; // pair candidate not found

      const PaTPar& tpar1 = ev.vParticle(ih1).ParInVtx(iv);
      const PaTPar& tpar2 = ev.vParticle(ih2).ParInVtx(iv);

      if( tpar1.Q() * tpar2.Q() > 0 ) continue; // skip pairs of the same sign


      //------- filling tree
      fHadrPair->Reset();
      fHadrPair->inPart = vtx.InParticle() >= 0;

      fHadrPair->distToPV = 0.0;
      float dvect[3],dist2=0.0;
      for(int i=0; i<3; ++i){// 3 component data

	//--------- positions
	fHadrPair->primVertPos[i] = fBPV->Pos(i);
	fHadrPair->vertPos[i] = vtx.Pos(i);
	dvect[i] = vtx.Pos(i) - fBPV->Pos(i);
	dist2 += dvect[i]*dvect[i];

	//--------- momenta
	fHadrPair->p1[i] = tpar1.Mom3()[i];
	fHadrPair->p2[i] = tpar2.Mom3()[i];

      } // end loop on 3 components


      if( dist2 > 0 ){
	fHadrPair->distToPV = sqrt(dist2); // distance to prim. vtx.

	// resolution of vertices
	int icov=0;
	float sigma2=0.0;
	for(int i=0; i<3; ++i){
	  for(int j=0; j<=i; ++j){
	    sigma2 +=
	      ( i==j ? 1.0:2.0 )*dvect[i]*( fBPV->Cov(icov) + vtx.Cov(icov) )*dvect[j];
	    ++icov;
	  }
	}
	sigma2 /= dist2;
	fHadrPair->sigmaPV = sigma2;
	fHadrPair->sigmaSV = sqrt(sigma2);
      }


      // hadron info
      fHadrPair->P1 = tpar1.Mom();
      fHadrPair->P2 = tpar2.Mom();
      fHadrPair->Q1 = tpar1.Q();
      fHadrPair->Q2 = tpar2.Q();

      // RICH info
      PaTPar tparRich;
      tpar1.Extrapolate(Zrich, tparRich);
      fHadrPair->thetaR1 = tparRich.Theta();
      tpar2.Extrapolate(Zrich, tparRich);
      fHadrPair->thetaR2 = tparRich.Theta();

      const PaTrack& tr1 =  ev.vTrack(itr1);
      fHadrPair->richInfo1 = PID->CheckRichInfo(tr1);
      const PaTrack& tr2 =  ev.vTrack(itr2);
      fHadrPair->richInfo2 = PID->CheckRichInfo(tr2);
      for(int i=0; i<6; ++i){
	if(fHadrPair->richInfo1) fHadrPair->LHs1[i] = PID->GetLike(i,tr1);
	if(fHadrPair->richInfo2) fHadrPair->LHs2[i] = PID->GetLike(i,tr2);
      }

      fHadrTree->Fill();

    }// end if NoutPart > 1

  }//--- end loop on vertices

  gNvtx=Nvtx;
  gNSecVtx /= Nvtx;
  fUtilityTree->Fill();
  //cout<<"check HadronPairs"<<endl;
  //cout<<"-------------------------------"<<endl;
} //end Find Hadron Paris

void LCAnalysis::InitHadronPairTree()
{
  fHadrPair = new HadronPairData;
  Phast::Ref().h_file->cd();
  fHadrTree = new TTree("hadrPairTree","Hadron pairs");
  fHadrTree->Branch("hadrPair","HadronPairData",&fHadrPair);

  fUtilityTree = new TTree("statTree","statTree");
  fUtilityTree->Branch("Nvtx",&gNvtx,"Nvtx/S");
  fUtilityTree->Branch("NSecVtx",&gNSecVtx,"NSecVtx/F");

  ghNvtx = new TH1D("hNvtx","hNvtx",40,-0.5,39.5);
  ghSecVtx = new TH2D("hSecVtx","hSecVtx",40,-0.5,39.5,40,-0.5,39.5);
//cout<<"check InitHadronPairTree"<<endl;
}



//------------------------------------------------ output methods
void LCAnalysis::PrintEventStats()
{

  int t=fTotNEvts,o=fTotNEvts,n;
  FILE *fout=fopen(fLogFileName,"w");
  printf("\n");
  for(int i=0; i<60; ++i)printf("-"); printf("\n");
  printf(" writing log output to file %s\n",fLogFileName);

  vector<LCEventSelectionCut*>::iterator it;
  char line[512];
  n=fTotNEvts;
  printf("test");
  printf("%-25s%10d%9.3f%%%9.3f%%  %p\n","Total",      n, 100.0*n/t,100.*(o-n)/o,fout);
  fprintf(fout,"%-25s%10d%9.3f%%%9.3f%%\n","Total",      n, 100.0*n/t,100.*(o-n)/o);
  for(it=fCutVector.begin(); it!=fCutVector.end(); ++it){
    (*it)->Output(line,fTotNEvts);
    cout<<line;
    fprintf(fout,"%s",line);
    o = (*it)->GetNpassed();
  }
  n=fTrigger;
  printf("%-25s%10d%9.3f%%%9.3f%%\n","Trigger",    n, 100.0*n/t,100.*(o-n)/o);
  fprintf(fout,"%-25s%10d%9.3f%%%9.3f%%\n","Trigger",    n, 100.0*n/t,100.*(o-n)/o);
  fclose(fout);

  for(int i=0; i<60; ++i)printf("-"); printf("\n");
  //cout<<"check PrintEventStats"<<endl;
}

void LCAnalysis::SaveHadronTree()
{
//   printf("\nLCAnalysis: saving hadron tree ... \n\n");
//   fEvtTree->AutoSave();
//   fHadrTree->AutoSave();
//   if(fIsMC){
//     fEvtTreeMC->AutoSave();
//     fHadrTreeMC->AutoSave();
//   }
//   fOutFile->Close();

  fDISEvtTree->AutoSave();
  //cout<<"First check SaveHadronTree"<<endl;
  //fDISEvtTree->GetCurrentFile()->Close(); //remove because segmentation falut
  //cout<<"Second check SaveHadronTree"<<endl;

  //----- Hadron selection stats
  int t=fTotNHadr,o=fTotNHadr,n;
  printf("\n");
  for(int i=0; i<50; ++i)printf("-"); printf("\n");
  n=fTotNHadr; printf("%-15s%10d%9.3f%%%9.3f%%\n","\"Hadrons\"",  n, 100.0*n/t,100.*(o-n)/o); o=n;
  n=fNsec    ; printf("%-15s%10d%9.3f%%%9.3f%%\n","Sec (no mu')", n, 100.0*n/t,100.*(o-n)/o); o=n;
  n=fNzlt1   ; printf("%-15s%10d%9.3f%%%9.3f%%\n","z < 1",        n, 100.0*n/t,100.*(o-n)/o); o=n;
  n=fNZfirst ; printf("%-15s%10d%9.3f%%%9.3f%%\n","Zfirst",       n, 100.0*n/t,100.*(o-n)/o); o=n;
  n=fNZlast  ; printf("%-15s%10d%9.3f%%%9.3f%%\n","Zlast",        n, 100.0*n/t,100.*(o-n)/o); o=n;
  n=fNXX0    ; printf("%-15s%10d%9.3f%%%9.3f%%\n","X/X0",         n, 100.0*n/t,100.*(o-n)/o); o=n;
  n=fNTrigHad;    printf("%-15s%10d%9.3f%%%9.3f%%\n","OT || inclMT", n, 100.0*n/t,100.*(o-n)/o);
  n=count_bestpv; printf("%-15s%10d%9.3f%%%9.3f%%\n","BPV",           n, 100.0*n/t,100.*(o-n)/o);;
  n=count_mup;    printf("%-15s%10d%9.3f%%%9.3f%%\n","mu'",       n, 100.0*n/t,100.*(o-n)/o);
  n=count_mupp;  printf("%-15s%10d%9.3f%%%9.3f%%\n","MC mu'is m+",   n, 100.0*n/t,100.*(o-n)/o);
  n=count_munp;  printf("%-15s%10d%9.3f%%%9.3f%%\n","MC mu'is m-",   n, 100.0*n/t,100.*(o-n)/o);
  n=count_mup0;  printf("%-15s%10d%9.3f%%%9.3f%%\n","MC mu0 is m+",  n, 100.0*n/t,100.*(o-n)/o);
  n=count_mun0;  printf("%-15s%10d%9.3f%%%9.3f%%\n","MC mu0 is m-",  n, 100.0*n/t,100.*(o-n)/o);

  for(int i=0; i<50; ++i)printf("-"); printf("\n");
  printf("No. of DIS events (without trigger cut): %10d\n",fTotNEvts);
  printf("                  (with    trigger cut): %10d\n",fTrigger);
  //cout<<"3rd check SaveHadronTree"<<endl;

  if( fIsMC ){
    cout<<"End of MC job: MC hadrons: "<<fNMChadr<<endl\
	<<"               multiple associations (~) : "<<fNMChadrMultAssoc<<endl;
  }
//cout<<"check SaveHadronTree"<<endl;
}




//--------------------------------------------------- Set methods
void LCAnalysis::SetQ2Min(const double& Q2min)
{
  fQ2min=Q2min;
  //cout<<"check SetQ2Min"<<endl;
}
void LCAnalysis::SetEbeamRange(const double& Emin,const double& Emax)
{
  fEbeamMin=Emin;
  fEbeamMax=Emax;
  //cout<<"check SetEbeamRange"<<endl;
}
void LCAnalysis::SetyRange(const double& ymin,const double& ymax)
{
  fyMin=ymin;
  fyMax=ymax;
  //cout<<"check SetyRange"<<endl;
}
void LCAnalysis::SetXX0max(const double& XX0max)
{
  fXX0max=XX0max;
  //cout<<"check SetXX0max \t "<<fXX0max<<endl;
}






//------------------------------------- Dynamic cuts
CheckMethod FindCheckMethod(std::string cutname)
{
  if( cutname == "BestPV")//  cout<<"found BestPV"<<endl;
    return &LCAnalysis::IsThereABestPV;
  if( cutname == "mu1Recons"     ) return &LCAnalysis::IsMu1Reconstructed;
  if( cutname == "intInTarget"   ) return &LCAnalysis::InteractionInTarget;
  if( cutname == "intInTarg2009" ) return &LCAnalysis::InteractionInTarget2009;
  if( cutname == "intInTarg2016" ) return &LCAnalysis::InteractionInTarget2016;
  if( cutname == "cellsCrossed"  ) return &LCAnalysis::CellsCrossed;
  if( cutname == "cellsCrossed2016"  ) return &LCAnalysis::CellsCrossed2016;
  //if( cutname == "hodos"         ) return &LCAnalysis::HodosCut;                  //cut mu' track passes trhough ineffcient slab of the hodoscope HM05 (added by Quiela)
  return 0;
  //cout<<"check FindCheckMethod"<<endl;
}



LCEventSelectionCut::LCEventSelectionCut(string name,LCAnalysis* analysis):
  fName(name),fTitle(name),fAnalysis(analysis),fNsample(0),fNpassed(0)
{
  fCheck = FindCheckMethod(name);
  //cout<<"check LCEventSelectionCut "<<endl;
}

bool LCEventSelectionCut::Check()
{
  bool check = DoCheck();
  ++fNsample;
  if(check) ++fNpassed;
  return check;
  //cout<<"check Check"<<endl;
}

bool LCEventSelectionCut::DoCheck()
{
  return (fAnalysis->*fCheck) ();
  //cout<<"check DoCheck"<<endl;
}

void LCEventSelectionCut::Output(char* output,int Ntot)
{
  sprintf(output,"%-25s%10d%9.3f%%%9.3f%%\n",
	  fTitle.data(),
	  fNpassed,
	  ((float)fNpassed)/Ntot*100.,
	  ((float)(fNsample-fNpassed))*100./fNsample);
  //cout<<"check Output"<<endl;
}



//------------------------------------- Dynamic range cuts
GetValueMethod FindGetValueMethod(std::string cutname)
{
  if( cutname == "Q2"    ) return &LCAnalysis::GetQ2;
  if( cutname == "Ebeam" ) return &LCAnalysis::GetEbeam;
  if( cutname == "y"     ) return &LCAnalysis::Gety;
  if( cutname == "xBj"   ) return &LCAnalysis::GetxBj;
  return 0;
  //cout<<"check FindGetValueMethod "<<endl;
}

LCEventSelectionRangeCut::LCEventSelectionRangeCut(string name,LCAnalysis* analysis,
						   double minVal,double maxVal):
  LCEventSelectionCut(name,analysis),fMinVal(minVal),fMaxVal(maxVal)
{
  fGetValue = FindGetValueMethod(name);

  if( fMaxVal < fMinVal ) fMaxVal = FLT_MAX;

  char title[512];
  if(fMaxVal < FLT_MAX )
    sprintf(title,"%-s (%.3f,%.3f)",fName.data(),fMinVal,fMaxVal);
  else
    sprintf(title,"%-s > %.3f",fName.data(),fMinVal);
  fTitle = title;
  //cout<<"check LCEventSelectionRangeCut"<<endl;
}

bool LCEventSelectionRangeCut::DoCheck()
{
  fVal = (fAnalysis->*fGetValue)();
  return fVal > fMinVal && fVal < fMaxVal;
  //cout<<"check LCEventSelectionRangeCut::DoCheck"<<endl;
}
