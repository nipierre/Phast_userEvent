#include <iostream>
#include <iomanip>

#include "TSpline.h"
#include "TRandom.h"

#include "TargetCell.h"
#include "Phast.h"
#include "PaAlgo.h"



/* \brief Gives the target location in space: shift and tilting.
  Returns false if no information for the given year.

  \param run  number of run ("-1" = MC)
  \param x   the shift of the target along x
  \param y   the shift of the target along y
  \param z_1 z position of the cell
  \param z_2 z position of the cell
  \param R    the radial cut
  \param yCUT the cut for upper part of the cell
*/
bool TargetCell::GetTargetLocation( int run,
		       double &x, double &y, double &z_1, double &z_2,
		       double &R, double &yCUT)
{
  if (run == -1)
  { // MC
    x = 0; y = 0;  z_1 = -325; z_2 = -71;
    R    = 1.9;
    yCUT = 1.2;
  }
  else if( run > 0 )
  { // 2016 DVCS - lh2
    x = 0; y = 0;  z_1 = -325; z_2 = -71;
    R    = 1.9;
    yCUT = 1.2;
  }
  else { // Not known
    cout<<"TargetCell::GetTargetLocation ==> target location for the run "<<run<<" is not known."<<endl;
    exit(1);
  }

  return true;
}


/* \brief The requirement that the muon beam trajectory crosses
  entirely the target cell. It is used in order to equalize fluxes
  through both cells.

  \param par    the beam track parameters in the primary vertex
  \param run    the run number ("-1" = MC)
  \param R_U    the user defined radial cut (R<1.5cm)
  \param yCUT_U the user defined vertical cut (y<1cm)
*/
bool TargetCell::CrossCells( PaTPar par, int run )
{
  PaTPar parE;
  double x,y,z_1,z_2,R,yCUT;

  if( run > 0 || run==-1)
  { // 1 cell
    if( !GetTargetLocation(run,x,y,z_1,z_2,R,yCUT) ) {
      cout<<"CrossCells PROBLEM: no info for the run "<<run<<" (1 cell)"<<endl;
      return false;
    }
  }

  //if( R_U    != -9999 ) R    = R_U;
  //if( yCUT_U != -9999 ) yCUT = yCUT_U;

  par.Extrapolate(z_1,parE,0);
  double xL = parE(1);
  double yL = parE(2);
  double rL = sqrt( (xL-x)*(xL-x) + (yL-y)*(yL-y) );
  if( rL > R ) return false;
  if( yL-y > yCUT ) return false;

  par.Extrapolate(z_2,parE,0);
  double xR = parE(1);
  double yR = parE(2);
  double rR = sqrt( (xR-x)*(xR-x) + (yR-y)*(yR-y) );
  if( rR > R ) return false;
  if( yR-y > yCUT ) return false;

  return true;
}



/* \brief The check for the primary vertex to be in one of the target cells.
  \param par    the beam track parameters in the primary vertex
  \param run    the run number ("-1" = MC)
  \param R_U    the user defined radial cut (R<1.4cm)
  \param yCUT_U the user defined vertical cut (y<1cm)
*/
bool TargetCell::InTarget( PaTPar par, int run )
{
  double x_1,y_1,z_1,z_2,R,yCUT;

  if( run > 0 && run==-1) {      // 2 cells
    if( !GetTargetLocation(run,x_1,y_1,z_1,z_2,R,yCUT) ) {
      cout<<"PaAlgo::InTarget() PROBLEM: no info for the run "<<run<<" (2 cells)"<<endl;
      return false;
    }
  }

  //if( R_U    != -9999 ) R    = R_U;
  //if( yCUT_U != -9999 ) yCUT = yCUT_U;

  double x = par(1);
  double y = par(2);
  double z = par(0);

  double r = sqrt( (x-x_1)*(x-x_1) + (y-y_1)*(y-y_1) );

  if( z < z_1 || z > z_2 ) return false;
  if(         r > R        ) return false;
  if(      y-y_1 > yCUT     ) return false;

  return true;
}
