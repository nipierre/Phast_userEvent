#include <iomanip>
#include <unistd.h>

#include "TargetCell.h"

TargetCell* TargetCell::ptr = NULL;

TargetCell& TargetCell::Instance()
{
  if( ptr == NULL ) ptr = new TargetCell();
  return *ptr;
}


TargetCell::TargetCell(): year(-1)
{

}


void TargetCell::Init()
{
  // Check if already initialized
  // if( year == e.Year() ) return;

  //char tstr[500];
  std::ifstream fin;
  std::string tstr;
  year = 2016;
  // if(e.IsMC()) tstr = "/sps/compass/npierre/PHAST/user/Target/target-mc.dat";
  else {
    if(year==2012) tstr = "/sps/compass/npierre/PHAST/user/Target/target-107924-109081.dat";
    if(year==2016) tstr = "/sps/compass/npierre/PHAST/user/Target/target-274508-274901-1.dat";
    if(year==2017) tstr = "/sps/compass/npierre/PHAST/user/Target/target-278473-278706-0.dat";
  }
  cout<<"Opening target cell description: "<<tstr<<"..."<<endl;
  fin.open(tstr.c_str());
  while(fin.is_open() && !fin.eof()) {
    float z, x, y, dummy;
    fin >> z >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy
        >> x >> y;
    zv.push_back(z);
    xv.push_back(x);
    yv.push_back(y);
    cout<<"  z="<<z<<"  x="<<x<<"  y="<<y<<endl;
  }
  cout<<"Target cell description loaded"<<endl;
  /*
  for(unsigned int i = 0; i < zv.size(); i++) {
    cout << zv[i] << "  "  << xv[i] << "  "  << yv[i] << endl;
  }*/
}


/*
  The total path length of a given beam track through the target cell is computed
  by summing up the individual contributions of each 10cm segment into which
  it has been divided.
  For each segment, we compute the values of the Z coordinate for which the distance
  between the two lines defined by the cell axis and the beam track is >= R.
  This reduces to solving a 2nd order equation of the form ax^2 + bx + c = 0.
  In this case, the three coefficients are as follows:
  Delta(dx/dz) = (dx/dz)_beam - (dx/dz)_cell_axis
  Delta(dy/dz) = (dy/dz)_beam - (dy/dz)_cell_axis
  Delta(x) = x_beam - x_cell_axis
  Delta(y) = y_beam - y_cell_axis
  a = Delta(dx/dz)^2 + Delta(dy/dz)^2
  b = 2*(Delta(dx/dz)*Delta(x) + Delta(dy/dz)*Delta(y))
  c = Delta(x)^2 + Delta(y)^2 - R^2
  Different cases are possible:
  - D = b*b - a*c*4.0f < 0
    in this case the two lines are never closer than R, or in other words the beam
    track is outside of the cell segment
 */
float TargetCell::PathLength(PaTrack* trk, float zmin, float zmax, float R)
{
  if( zv.size() != xv.size() ) return 0;
  if( zv.size() != yv.size() ) return 0;
  if( zv.size() < 2 ) return 0;

  double length = 0;
  PaTPar hout1, hout2;
  //bool stop = false;
  for( unsigned int i = 0; i < zv.size()-1; i++ ) {
    double z1 = zv[i];
    double z2 = zv[i+1];
    if( z2 <= zmin ) continue;
    if( z1 >= zmax ) continue;
    //double dz = z2-z1;

    double xc1 = xv[i];
    double xc2 = xv[i+1];

    double yc1 = yv[i];
    double yc2 = yv[i+1];

    double dxcdz = (xc2-xc1)/(z2-z1);
    double dycdz = (yc2-yc1)/(z2-z1);

    if( z1 < zmin ) {
      z1 = zmin;
      xc1 = xv[i] + dxcdz*(z1-zv[i]);
      yc1 = yv[i] + dycdz*(z1-zv[i]);
    }

    if( z2 > zmax ) {
      z2 = zmax;
      xc2 = xv[i] + dxcdz*(z2-zv[i]);
      yc2 = yv[i] + dycdz*(z2-zv[i]);
    }

    trk->Extrap(z1, hout1);
    double xt1 = hout1.X();
    double yt1 = hout1.Y();
    double dxtdz = hout1.dXdZ();
    double dytdz = hout1.dYdZ();

    trk->Extrap(z2, hout2);
    double xt2 = hout2.X();
    double yt2 = hout2.Y();

    double Dx = xt1 - xc1;
    double Dy = yt1 - yc1;
    double Dx2 = xt2 - xc2;
    double Dy2 = yt2 - yc2;
    double ddx = dxtdz - dxcdz;
    double ddy = dytdz - dycdz;

    double a = ddx*ddx + ddy*ddy;
    double b = (ddx*Dx + ddy*Dy)*2.0f;
    double c = Dx*Dx + Dy*Dy - R*R;
    double D = b*b - a*c*4.0f;

    /*
    cout<<endl<<"z1="<<z1<<"  z2="<<z2<<endl;
    cout<<"xc1="<<xc1<<"  yc1="<<yc1<<endl;
    cout<<"xc2="<<xc2<<"  yc2="<<yc2<<endl;
    cout<<"R="<<R<<endl;
    cout<<"xt1="<<xt1<<"  yt1="<<yt1<<endl;
    cout<<"xt2="<<xt2<<"  yt2="<<yt2<<endl;
    cout<<"Dx="<<Dx<<"  Dy="<<Dy<<endl;
    cout<<"Dx2="<<Dx2<<"  Dy2="<<Dy2<<endl;
    cout<<"ddx="<<ddx<<"  ddy="<<ddy<<endl;
    cout<<"a="<<a<<"  b="<<b<<"  c="<<c<<"  D="<<D<<endl;
    */

    //double Rbegin = sqrt(pow(xt1-xc1,2)+pow(yt1-yc1,2));
    //double Rend = sqrt(pow(xt2-xc2,2)+pow(yt2-yc2,2));

//     if( D < 0 && (Rbegin<R || Rend<R))
//       {
//  cout << "D < 0:" << endl;
//  cout<<"z1="<<z1<<"  z2="<<z2<<endl;
//  cout<<"xc1="<<xc1<<"  yc1="<<yc1<<endl;
//  cout<<"xc2="<<xc2<<"  yc2="<<yc2<<endl;
//  cout<<"R="<<R<<endl;
//  cout<<"xt1="<<xt1<<"  yt1="<<yt1<<endl;
//  cout<<"xt2="<<xt2<<"  yt2="<<yt2<<endl;
//  cout<<"Dx="<<Dx<<"  Dy="<<Dy<<endl;
//  cout<<"Dx2="<<Dx2<<"  Dy2="<<Dy2<<endl;
//  cout<<"ddx="<<ddx<<"  ddy="<<ddy<<endl;
//  cout<<"a="<<a<<"  b="<<b<<"  c="<<c<<"  D="<<D<<endl;
//       }

    if( D < 0 ) continue;

    if(ddx==0.0 && ddy==0.0)// track parallel to cell axis
      {
  if(c>0){
    continue;
  }else{
    double dZ = z2 - z1;
    double dx = dxtdz*dZ;
    double dy = dytdz*dZ;
    double r = TMath::Sqrt(dx*dx + dy*dy);
    double L = TMath::Sqrt(r*r + dZ*dZ);
    length += L;

    continue;
  }
      }

    double Z1 = (-1.0f*b - TMath::Sqrt(D))/(a*2.0f) + z1;
    double Z2 = (-1.0f*b + TMath::Sqrt(D))/(a*2.0f) + z1;

    //cout<<"Z1="<<Z1<<"  Z2="<<Z2<<endl;

    if( (Z1<z1) && (Z2<z1) ) continue;
    if( (Z1>z2) && (Z2>z2) ) continue;

    if( Z1 < z1 ) Z1 = z1;
    if( Z2 > z2 ) Z2 = z2;

    //cout<<"Z1="<<Z1<<"  Z2="<<Z2<<endl;

    double dZ = Z2 - Z1;
    // The track crosses the whole segment
    double dx = dxtdz*dZ;
    double dy = dytdz*dZ;
    double r = TMath::Sqrt(dx*dx + dy*dy);

    double L = TMath::Sqrt(r*r + dZ*dZ);
    //cout<<"r="<<r<<"  dZ="<<dZ<<"    L="<<L<<endl<<endl;
    length += L;

    //DEBUG
//     if(L>0 && Rbegin>R && Rend>R)
//       {
//  cout << "L > 0" << endl;
//  cout<<"z1="<<z1<<"  z2="<<z2<<endl;
//  cout<<"xc1="<<xc1<<"  yc1="<<yc1<<endl;
//  cout<<"xc2="<<xc2<<"  yc2="<<yc2<<endl;
//  cout<<"R="<<R<<endl;
//  cout<<"xt1="<<xt1<<"  yt1="<<yt1<<endl;
//  cout<<"xt2="<<xt2<<"  yt2="<<yt2<<endl;
//  cout<<"Dx="<<Dx<<"  Dy="<<Dy<<endl;
//  cout<<"Dx2="<<Dx2<<"  Dy2="<<Dy2<<endl;
//  cout<<"ddx="<<ddx<<"  ddy="<<ddy<<endl;
//  cout<<"a="<<a<<"  b="<<b<<"  c="<<c<<"  D="<<D<<endl;

//  cout<<"Z1="<<Z1<<"  Z2="<<Z2<<endl;

//  cout<<"r="<<r<<"  dZ="<<dZ<<"    L="<<L<<endl<<endl;
//  cout<<"Rbegin="<<Rbegin<<"   Rend="<<Rend<<endl;
//       }

//     if(L>=10.0 && ((Rbegin<R && Rend>R) || (Rbegin>R && Rend<R)))
//     if(L==0.0 && ((Rbegin<R && Rend>R) || (Rbegin>R && Rend<R)))
//       {
//  cout << "L == 0" << endl;
//  cout<<"z1="<<z1<<"  z2="<<z2<<endl;
//  cout<<"xc1="<<xc1<<"  yc1="<<yc1<<endl;
//  cout<<"xc2="<<xc2<<"  yc2="<<yc2<<endl;
//  cout<<"R="<<R<<endl;
//  cout<<"xt1="<<xt1<<"  yt1="<<yt1<<endl;
//  cout<<"xt2="<<xt2<<"  yt2="<<yt2<<endl;
//  cout<<"Dx="<<Dx<<"  Dy="<<Dy<<endl;
//  cout<<"Dx2="<<Dx2<<"  Dy2="<<Dy2<<endl;
//  cout<<"ddx="<<ddx<<"  ddy="<<ddy<<endl;
//  cout<<"a="<<a<<"  b="<<b<<"  c="<<c<<"  D="<<D<<endl;

//  cout<<"Z1="<<Z1<<"  Z2="<<Z2<<endl;

//  cout<<"r="<<r<<"  dZ="<<dZ<<"    L="<<L<<endl<<endl;
//    cout<<"Rbegin="<<Rbegin<<"   Rend="<<Rend<<endl;
//       }
  }
  //cout<<"length="<<length<<endl<<endl;
  //if(stop) getchar();
  return length;
}


bool TargetCell::InTarget(const PaVertex& vtx, float R)
{
  float xvtx = vtx.X();
  float yvtx = vtx.Y();
  float zvtx = vtx.Z();

  float xc, yc;
  CellCenter(zvtx, xc, yc);
  double dx = xvtx-xc;
  double dy = yvtx-yc;
  double r = TMath::Sqrt(dx*dx + dy*dy);

  //if(r>R)cout << "Z = "<< zvtx << " R = " << r << endl;
  /*
  cout<<"TargetCell::InTarget()"<<endl
      <<"  vertex: "<<xvtx<<"  "<<yvtx<<"  "<<zvtx<<endl
      <<"  cell:   "<<xc<<"  "<<yc<<endl
      <<"  dist:   "<<dx<<"  "<<dy<<"  "<<r<<endl;
  */
  return( r <= R && yvtx < 1.2);
}


bool TargetCell::CrossCells(const PaTrack& trk, float zmin, float zmax, float R, float Y)
{
  if( zv.size() != xv.size() ) return 0;
  if( zv.size() != yv.size() ) return 0;
  if( zv.size() < 2 ) return 0;

  if( zmin < zv[0] ) zmin = zv[0];
  if( zmax > zv[zv.size()-1] ) zmax = zv[zv.size()-1];

  double R2 = R*R;
  PaTPar hout1;

  for( unsigned int i = 1; i < zv.size(); i++ ) {
    double z0 = zv[i-1];
    double z1 = zv[i];
    if( z1 <= zmin ) continue;
    if( z0 >= zmax ) continue;



    double xc0 = xv[i-1];
    double yc0 = yv[i-1];
    double xc1 = xv[i];
    double yc1 = yv[i];


    if( zmin > z0 && zmin < z1) {
      // Zmin is inside the current segment, so we need to interpolate
      float xc, yc;
      CellCenter( zmin, xc, yc );

      trk.Extrap(zmin, hout1);
      double xt1 = hout1.X();
      double yt1 = hout1.Y();

      double dx = xt1 - xc;
      double dy = yt1 - yc;
      double r2 = dx*dx + dy*dy;
      if(yt1>=Y) return false;
      if( r2 > R2 ) return false;
    }

    if( zmax > z0 && zmax < z1) {
      // Zmax is inside the current segment, so we need to interpolate
      float xc, yc;
      CellCenter( zmax, xc, yc );

      trk.Extrap(zmax, hout1);
      double xt1 = hout1.X();
      double yt1 = hout1.Y();

      double dx = xt1 - xc;
      double dy = yt1 - yc;
      double r2 = dx*dx + dy*dy;
      if(yt1>=Y) return false;
      if( r2 > R2 ) return false;
    }

    // Extrapolate to the current segment
    trk.Extrap(z1, hout1);
    double xt1 = hout1.X();
    double yt1 = hout1.Y();

    double dx = xt1 - xc1;
    double dy = yt1 - yc1;
    double r2 = dx*dx + dy*dy;
    if( r2 > R2 || yt1>=Y ) return false;
  }
  return true;
}


void TargetCell::CellCenter(float z, float& xc, float& yc)
{
  xc = 1000000;
  yc = 1000000;
  for( int i = 0; i < zv.size()-1; i++ ) {
    double z1 = zv[i];
    double z2 = zv[i+1];
    if( z2 < z ) continue;
    if( z1 > z ) continue;

    double xc1 = xv[i];
    double xc2 = xv[i+1];

    double yc1 = yv[i];
    double yc2 = yv[i+1];

    double dxcdz = (xc2-xc1)/(z2-z1);
    double dycdz = (yc2-yc1)/(z2-z1);

    double dz = z-z1;
    xc = xc1 + dxcdz*dz;
    yc = yc1 + dycdz*dz;

    break;
  }
}
