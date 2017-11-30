
#ifndef DISEventTreeManager_h 
#define DISEventTreeManager_h 1 

#include "DISEventData.h"

class TFile;
class TChain;

class DISEventHandler;

class DISEventTreeManager{
  
  TChain *tree;
  DISEventData* dis;
  std::vector<HadronData> hadrons;
  std::vector<HadronData>* hadronsPtr;

  bool isMCtree;
  bool isMCgenTree;  // contains only generator info
  DISEventMCData* MCdis;
  std::vector<HadronMCData> MChadrons;
  std::vector<HadronMCData>* MChadronsPtr;

public:
  DISEventTreeManager():tree(0),dis(0),hadronsPtr(&hadrons),
			isMCtree(false),MCdis(0),MChadronsPtr(&MChadrons),
			entry(0){};

  void OpenTree(const char *treename,const char* filename);
  void LoopOnTree(DISEventHandler* handler,int NmaxEntries=-1);

private:
  long entry;

};


class DISEventHandler{
public:
  virtual ~DISEventHandler(){}
  virtual void HandleEvent(DISEventData* disPtr, std::vector<HadronData>* hdrnsPtr,
			   DISEventMCData* mcdisPtr, std::vector<HadronMCData>* mchdrnsPtr) = 0;
 
};


#endif
