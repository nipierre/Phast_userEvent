#include "Phast.h"
#include "PaEvent.h"
#include "PaAlgo.h"
#include "TH1D.h"

#include "LCAnalysis.h"
using namespace std;
static LCAnalysis *gAnalysis;
void UserEvent20(PaEvent& e)
{

//NLUDATA ld;
//vector<LUJET> lujets;         // create LUJET structure
  // int nparticles;               // number of particles
   //e.MCgen(nparticles,lujets);
//if(nparticles<9)cout<<"elastic"<<endl;
//for(int i=0;i<nparticles;i++){
//cout<<i<<" "<<lujets[i].k[1]<<endl;
//cout<<i<<endl;
//    cout<<"LUJET: lu2kine = "<<lujets[i].lu2kine<<endl
//	<<"  k[5] = "<<lujets[i].k[0]<<", "<<lujets[i].k[1]<<", "<<lujets[i].k[2]<<", "<<lujets[i].k[3]<<", "<<lujets[i].k[4]<<endl
//	<<"  p[5] = "<<lujets[i].p[0]<<", "<<lujets[i].p[1]<<", "<<lujets[i].p[2]<<", "<<lujets[i].p[3]<<", "<<lujets[i].p[4]<<endl;
//}
//cout<<e.MCgen(ld)<<endl;


  static bool first=true;
  if(first){
    gAnalysis = new LCAnalysis;
    first = false;
  }
  gAnalysis->DoEvent(e);
 
}

void UserJobEnd20()
{
  gAnalysis->JobEnd();;
}

