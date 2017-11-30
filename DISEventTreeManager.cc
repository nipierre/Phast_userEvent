

#include "DISEventTreeManager.h"

#include "TFile.h"
#include "TChain.h"
#include <iostream>

using namespace std;



void DISEventTreeManager::OpenTree(const char *treename,const char* filename)
{

  cout<<"DISEventTreeManager: opening tree "<<treename
      <<" from chain "<<endl
      <<"   "<<filename<<endl;

  tree = new TChain(treename);
  tree->Add(filename);

  isMCgenTree = !tree->FindBranch("DISEvt");
  
  if(dis) delete dis;
  if(!isMCgenTree){
    dis = new DISEventData;
    tree->SetBranchAddress("DISEvt",&dis);
    hadrons.clear();
    tree->SetBranchAddress("Hadrons",&hadronsPtr);
  }else{
    dis = 0;
    hadronsPtr = 0;
  }

  isMCtree = tree->FindBranch("DISMCEvt");

  if(MCdis) delete MCdis;
  if(isMCtree){
    MCdis = new DISEventMCData;
    tree->SetBranchAddress("DISMCEvt",&MCdis);
    MChadrons.clear();
    MChadronsPtr = &MChadrons;
    tree->SetBranchAddress("MCHadrons",&MChadronsPtr);
  }else{
    MCdis=0;
    MChadronsPtr=0;
  }


  entry=0;
}

void DISEventTreeManager::LoopOnTree(DISEventHandler* handler,int NmaxEntries)
{
  if( NmaxEntries < 0 ) NmaxEntries=TChain::kBigNumber;

  entry=0;
  while(entry<NmaxEntries && tree->GetEntry(entry)){
    
    if(dis) dis->CalcKin();
    if(MCdis) MCdis->CalcKin();
    handler->HandleEvent(dis,hadronsPtr,MCdis,MChadronsPtr);

    //    cout<<entry<<" qua"<<endl;

    ++entry;
  }

  cout<<"DISEventTreeManager::LoopOnTree: "<<entry<<" entries read in tree "<<tree->GetName()<<endl
      <<"                                 from "<<tree->GetTreeNumber()+1<<" files"<<endl;


}
