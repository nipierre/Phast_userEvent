#include "Phast.h"
#include "PaEvent.h"
#include "TTree.h"

using namespace std;

namespace {
  double xVtx;
  double yVtx;
  double zVtx;
  double xPos;
  double yPos;

  TTree *tree;

  void
  init()
  {
    TTree* t = tree = new TTree("GhostTracks", "Ghost Tracks Study");

    t->Branch("xVtx", &zVtx, "zVtx/D");
    t->Branch("yVtx", &zVtx, "zVtx/D");
    t->Branch("zVtx", &zVtx, "zVtx/D");

    t->Branch("xPos", &xPos, "xPos/D");
    t->Branch("yPos", &yPos, "yPos/D");
  }
}

void UserEvent1992(PaEvent& ev)
{
  static bool first = true;
  if (first)
  {
    first = false;
    init();
  }

  const PaVertex& v = ev.vVertex(ev.iBestPrimaryVertex());
  set<int> tracklist;

  xVtx = v.X();
  yVtx = v.Y();
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
