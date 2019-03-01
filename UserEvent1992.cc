#include "Phast.h"
#include "PaEvent.h"
#include "PaAlgo.h"
#include "TTree.h"
#include <fstream>
#include <sstream>

#include "G3part.h"

using namespace std;

static const double M_mu = G3partMass[5];

namespace {
  int RunNumber;
  int ChunkNumbr;
  int SpillNumber;
  int EvSpill;
  double mQ2;
  double x;
  double y;
  double mW;
  int trigMask;
  double p;
  double z;
  double xVtx;
  double yVtx;
  double zVtx;
  double xPos;
  double yPos;
  int PotGhost;

  TTree *tree;

  void
  init()
  {
    TTree* t = tree = new TTree("GhostTracks", "Ghost Tracks Study");

    t->Branch("RunNumber", &RunNumber, "RunNumber/I");
    t->Branch("ChunkNumbr", &ChunkNumbr, "ChunkNumbr/I");
    t->Branch("SpillNumber", &SpillNumber, "SpillNumber/I");
    t->Branch("EvSpill", &EvSpill, "EvSpill/I");

    t->Branch("mQ2", &mQ2, "mQ2/D");
    t->Branch("x", &x, "x/D");
    t->Branch("y", &y, "y/D");
    t->Branch("mW", &mW, "mW/D");
    t->Branch("trigMask", &trigMask, "trigMask/I");
    t->Branch("p", &p, "p/D");
    t->Branch("z", &z, "z/D");

    t->Branch("xVtx", &xVtx, "xVtx/D");
    t->Branch("yVtx", &yVtx, "yVtx/D");
    t->Branch("zVtx", &zVtx, "zVtx/D");

    t->Branch("xPos", &xPos, "xPos/D");
    t->Branch("yPos", &yPos, "yPos/D");

    t->Branch("PotGhost", &PotGhost, "PotGhost/I");
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

  RunNumber = ev.RunNum();
  ChunkNumbr = ev.ChunkNumber();
  SpillNumber = ev.SpillNum();
  EvSpill = ev.EvInSpill();

  int fiBPV = ev.iBestPrimaryVertex();

  const PaVertex& v = ev.vVertex(fiBPV);
  set<int> tracklist;

  xVtx = v.X();
  yVtx = v.Y();
  zVtx = v.Z();

  int imu0=-1, imu1=-1;
  imu0 = v.InParticle();
  imu1 = v.iMuPrim(0,1,1,1,30);

  if(fiBPV==-1 || imu1==-1 || imu0==-1) return;

  const PaParticle& Mu0   = ev.vParticle(imu0); // the beam muon
  const PaParticle& Mu1   = ev.vParticle(imu1); // the scattered muon
  const PaTPar& ParamMu0  = Mu0.ParInVtx(fiBPV); // fitted mu  parameters in the primary vertex
  const PaTPar& ParamMu1  = Mu1.ParInVtx(fiBPV); // fitted mu' parameters in the primary vertex
  TLorentzVector kMu0 = ParamMu0.LzVec(M_mu); // beam      mu Lorentz vector
  TLorentzVector kMu1 = ParamMu1.LzVec(M_mu); // scattered mu Lorentz vector
  int itr0 = Mu0.iTrack();
  int itr1 = Mu1.iTrack();
  const PaTrack& track0 = ev.vTrack(itr0);
  const PaTrack& track1 = ev.vTrack(itr1);

  double Ebeam = kMu0.E();
  int fEbeam = (140<Ebeam && Ebeam<180) ? 1 : 0;
  mQ2 = PaAlgo::Q2(kMu0, kMu1);
  int fQ2    = mQ2>1 ? 1 : 0;
  double nu    = (Ebeam - kMu1.E());
  y = nu/Ebeam;
  int fy     = (0.1<y && y<0.7) ? 1 : 0;
  x = PaAlgo::xbj(kMu0,kMu1);
  int fxBj   = (0.004<x && x<0.4) ? 1 : 0;
  mW = sqrt(PaAlgo::W2(kMu0,kMu1));
  int fW2    = (5<mW && mW<17) ? 1 : 0;
  trigMask = ev.TrigMask();
  int TrigOk = trigMask&2 || trigMask&4 || trigMask&8 || trigMask&512;

  int fBMS = track0.NHitsFoundInDetect("BM")>3 ? 1 : 0;
  int fChi2beam = track0.Chi2tot()/float(track0.Ndf())<10 ? 1 : 0;
  int fChi2muprim = track1.Chi2tot()/float(track1.Ndf())<10 ? 1 : 0;
  int fMZfirst = track1.ZFirst()<350 ? 1 : 0;
  int fInTarget = PaAlgo::InTarget(ParamMu0,'O',ev.RunNum(),1.9,1.2,-325,-71,1.9);
  int fCellsCrossed = PaAlgo::CrossCells(ParamMu0,ev.RunNum(),1.9,1.2,-325,-71,1.9);

  if(!(fQ2 && fy && fxBj && fW2 && TrigOk && fBMS && fChi2beam && fChi2muprim
      && fMZfirst && fInTarget && fCellsCrossed && fEbeam)) return;

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
    PaTPar param;
    tr.Extrapolate(v.Z(),param);
    p = param.Mom();
    PaTPar tParRich;
    tr.Extrapolate(615.6, tParRich);
    double RICHx=tParRich.Pos(0);
    double RICHy=tParRich.Pos(1);
    z = sqrt(pow(p,2)+pow(0.13957018,2))/nu;

    int fhchi2 = tr.Chi2tot()/float(tr.Ndf())<10 ? 1 : 0;
    int fXX0 = tr.XX0()<15 ? 1 : 0;
    int fZfirst = tr.ZFirst()<350 ? 1 : 0;
    int fZlast = tr.ZLast()>350 ? 1 : 0;
    int fp = (12 < p && p < 40) ? 1 : 0;
    int fthRICH = (0.01<tParRich.Theta(false) && tParRich.Theta(false)<0.12) ? 1 : 0;
    int fpipe = (pow(RICHx,2)+pow(RICHy,2))>25 ? 1 : 0;
    int fz = (0.2<z && z<0.85) ? 1 : 0;

    if(!(fhchi2 && fXX0 && fZfirst && fZlast && fp && fthRICH && fpipe
        && fz )) continue;

    xPos = param(1)-v.X();
    yPos = param(2)-v.Y();

    if(tracklist.find(i) != tracklist.end())
    {
      tracklist.erase(tracklist.find(i));
      PotGhost = 0;
    }
    else
    {
      PotGhost = 1;
      if(-110<zVtx && zVtx<-70)
      {
        ofstream ofs_debug("/sps/compass/npierre/Ghost/debug.txt", ofstream::out | ofstream::app);
        ofs_debug << ev.RunNum() << " " << ev.SpillNum() << " " << ev.EvInSpill() << " " << ev.ChunkNumber() << endl;
        ofs_debug.close();
      }
    }

    tree->Fill();
  }
}
