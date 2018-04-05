#include "PaEvent.h"

#include <vector>
#include <string>

#include "DISEventData.h"

class TFile;
class TTree;
class TNtuple;
class TEventList;
class PaPid;
class PaCalorimeter;

class LCEventSelectionCut;
class LCEventSelectionRangeCut;

//class DISEventData;
//class DISEventMCData;
class HadronPairData;

class TargetCell;

class LCAnalysis{

  //--- type of analysis from config file
  int fAnalysisType;

  //--- event selection utility flags
  bool fEventAccepted;
  bool fSaveEvents;

  //--- event selection stats.
  int fTotNEvts;
  int fTrigger;

  //--- vector of dynamic event selection cuts
  std::vector<LCEventSelectionCut*> fCutVector;

  //--- event selection variables
  double fQ2min;              // default = 1.0
  double fEbeamMin,fEbeamMax; // default = 140 ... 180 GeV
  double fyMin,fyMax;         // default = 0.1 ... 0.9
  double fxBjMin,fxBjMax;     // default = 0.004 ... 0.7

  //--- hadron selection stats
  int fTotNHadr;
  int fNsec;
  int fNzlt1;
  int fNZfirst;
  int fNZlast;
  int fNXX0;
  int fNTrigHad;
  //--- for MC
  int fNMChadr, fNMChadrMultAssoc;


  //--- hadron selection variables
  double fZref; // just one value now for comparing with Zfirst and
	       // Zlast, value fixed to default: 350
  double fXX0max;// default = 30

  //--- event infos
  PaEvent* fEvent;
  bool fIsMC;
  const PaVertex* fBPV;
  int fiBPV;
  double fzVTX;
  int fimu0,fimu1;
  bool fReconsEvent;
  bool fValidMu;
  int fTrigMask;
  bool fTrigOk;
  bool fCellsCrossed;
  bool fChi2CutFlag;
  DISEventData* fDISEvt;
  DISEventMCData* fDISMCEvt;
  std::vector<HadronData> fHadrons;
  std::vector<HadronData>* fHadronsPtr;
  std::vector<HadronMCData> fMCHadrons;
  std::vector<HadronMCData>* fMCHadronsPtr;
  //-- utility variables
  int HM04h,HM05h,HL04h,HL05h,HO03h,HO04h,HG01h,HG02h;
  double HM04x,HM04y,HM04z,HM05x,HM05y,HM05z;
  double HL04x,HL04y,HL04z,HL05x,HL05y,HL05z;
  double HO03x,HO03y,HO03z,HO04x,HO04y,HO04z;
  double HG01x,HG01y,HG01z,HG021x,HG021y,HG021z,HG022x,HG022y,HG022z;
  double MM01x,MM01y,MM01z,MM02x,MM02y,MM02z,MM03x,MM03y,MM03z;
  double Z2Ax,Z2Ay,Z2Az,Z2Bx,Z2By,Z2Bz;
  double RICHx,RICHy; //added by Quiela

  //--- kinematical quantities
  TLorentzVector fkMu0, fkMu1;
  double fEbeam,fQ2,fnu,fy,fxBj,fW2,fz;
  double fEmu1,fThmu1,fPhmu1;

  //---- track infos
  double fXX0mu1;
  double fZfirst,fZlast,fXX0;
  short fCharge;
  double fMom,fthP,fphP,fthRICH;
  //double fMom,fthP,fphP,fthRICH;
  //---Target 2012
  TargetCell* fTcell;
  //--- Likelihoods
  PaPid* fPid;
  short fLikePid,fMCPid,fTruePidMC;
  double fLpi,fLK,fLp,fLe,fLmu,fLBg;

  //--- output saving members
  TFile* fOutFile;
  TTree* fEvtTree;
  TTree* fEvtTreeMC;
  TTree* fHadrTree;
  TTree* fHadrTreeMC;
  HadronPairData* fHadrPair;

  TNtuple* fHistInfo;
  TTree* fDISEvtTree;
  TTree* fDISMCEvtTree;
  TEventList* fDISEvtList;
  TTree* fUtilityTree;

  //--- log file name
  char fLogFileName[512];

  //--- setup info
  int fdatatargetType;
  int fMCtargetType;

  //--- pointers to hadronic calorimeters
  PaCalorimeter *fHadrCal1, *fHadrCal2;

public:
  LCAnalysis();
  ~LCAnalysis();

  void DoEvent(PaEvent& ev); // General Event handler
  void SetMC(bool isMC);
  void JobEnd(); // Output


protected:
  //---------------------------------------------- Pointer to selected handler
  typedef void (LCAnalysis::*EventHandler) (PaEvent&);
  EventHandler fEvHandler;

  //---------------------------------------------- Initialisation
  void ReadCutFile();

  //---------------------------------------------- Specific handlers
  void SelectEvent(PaEvent& ev);
  void FindHadrons(PaEvent& ev);
  void SelectHadronPairs(PaEvent& ev);

  //---------------------------------------------- Utility methods
  void ResetEventInfo();
  void ReadMinimumEventInfo();
  void CheckTriggerMask(const PaEvent& ev);
  void SetMuKinematics(const PaEvent& ev,const int& iVtx,
		       const int& imu0,const int& imu1);
  void Calcz(const PaTrack& tr);
  void CopyDISEvtData(int pReconsEvent);
  void InitHadronTree();
  int FindHadronsMC(const PaMCvertex&);
  void InitHadronPairTree();
  double GetMassPid(int pid) const;
  bool parseMCgen(const vector<PaMCgen> &vMCgen,
		  TLorentzVector &gammTr);
public:
  //---------------------------------------------- Set methods
  void SetQ2Min(const double& Q2min);
  void SetEbeamRange(const double& Emin,const double& Emax);
  void SetyRange(const double& ymin,const double& ymax);
  void SetXX0max(const double& XX0max);

  //---------------------------------------------- Check methods (for dynamic cuts)
  bool IsThereABestPV();
  bool IsMu1Reconstructed();
  bool IsMu0Valid(PaTPar par);
  bool InteractionInTarget();
  bool InteractionInTarget2009();
  bool InteractionInTarget2016();
  bool CellsCrossed();
  bool CellsCrossed2016();

  //---------------------------------------------- Get methods (for dynamic cuts)
  const double& GetQ2()    const {return fQ2;};
  const double& GetEbeam() const {return fEbeam;};
  const double& Gety()     const {return fy;};
  const double& GetxBj()   const {return fxBj;};

protected:
  //---------------------------------------------- Specific end of job methods
  void PrintEventStats();
  void SaveHadronTree();

};


//------------------------------------ pointer to LCAnalysis check method
typedef bool (LCAnalysis::*CheckMethod) ();

//------------------------------------ cut name -> CheckMethod mapping function
CheckMethod FindCheckMethod(std::string cutname);

class LCEventSelectionCut{ //-------- general selection cut

protected:
  std::string fName;
  std::string fTitle;
  LCAnalysis* fAnalysis;

private:
  CheckMethod fCheck;

protected:
  int fNsample; // number of checks done in the session
  int fNpassed; // number of passed checks

public:
  LCEventSelectionCut(std::string name,LCAnalysis* analysis);
  virtual ~LCEventSelectionCut() {};

  //--- interface for LCAnalysis
  bool Check();
  const int& GetNpassed() const {return fNpassed;};

  //--- formats stats output for a session. Ntot is the total number
  //--- of events processed.
  virtual void Output(char* output,int Ntot);

  //--- Get methods
  const std::string& GetName(){return fName;};

protected:
  //--- specific implementation of the check
  virtual bool DoCheck();

};


//------------------------------------ pointer to LCAnalysis get value method
typedef const double& (LCAnalysis::*GetValueMethod) () const;

//------------------------------------ cut name -> GetValueMethod mapping function
GetValueMethod FindGetValueMethod(std::string cutname);

//---- selection cut checking if a parameter lies in a given interval
// Checks supported: 1. value within a range (min < x < max) and
// 2. value larger than a lower threshold (x > min). If the maximum
// value is left blank in the configuration file or it is smaller than
// the minimum value, it is set to FLT_MAX (no maximum check).
//
class LCEventSelectionRangeCut: virtual public LCEventSelectionCut{

private:
  GetValueMethod fGetValue;

  double fMinVal; // minimum allowed value
  double fMaxVal; // maximum allowed value
  double fVal;    // value to be checked for current event

public:
  LCEventSelectionRangeCut(std::string name,LCAnalysis* analysis,
			   double minVal,double maxVal);

protected:
  //--- specific implementation of the check
  virtual bool DoCheck();

};
