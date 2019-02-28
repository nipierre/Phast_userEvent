#include "Phast.h"
#include "PaEvent.h"
#include "TTree.h"

using namespace std;

namespace {
  Double_t zVtx;
  Double_t xPos;
  Double_t yPos;

  TTree *tree;

  void
  init()
  {
    TTree* t = tree = new TTree("GhostTracks", "Ghost Tracks Study");

    t->Branch("zVtx", "Double_t", &zVtx);

    t->Branch("xPos", "Double_t", &xPos);
    t->Branch("yPos", "Double_t", &yPos);
  }
}

void UserEvent1992(PaEvent& e)
{

  const PaVertex& v = ev.vVertex(ev.iBestPrimaryVertex());
  set<int> tracklist;

  zVtx = v.Z();

  for(int i=0; i<v.NOutParticles(); i++)
  {
    const PaParticle& pa = ev.vParticle(v.iOutParticle(i));
    if(pa.iTrack()!=-1)
      tracklist.insert(pa.iTrack());
    else
      cout << "Particle track index does not exist !" << endl;
  }

  for(int i=0; i<ev.NTrack(); i++)
  {
    const PaTrack& tr = ev.vTrack(i);
    if(tracklist.find(i) != tracklist.end())
    {
      tracklist.erase(tracklist.find(i));
    }
    else
    {
      PaTPar parH;
      tr.Extrapolate(v.Z(),parH);
      xPos = parH(1)-v.X();
      yPos = parH(2)-v.Y();
      tree->Fill();
    }
  }
}
