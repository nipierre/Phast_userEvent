#ifndef TARGETCELL_h
#define TARGETCELL_h


/*!
  \class TargetCell
  \brief Target position and checks if interaction's done in target

*/

#include <cmath>
#include <TVector3.h>
#include <TLorentzVector.h>
#include "PaEvent.h"

using namespace std;


class TargetCell

{

public:

/*! \brief Gives the target location in space: shift and tilting.
  Returns false if no information for the given year.

  \param run  number of run ("-1" = MC)
  \param x   the shift of the target along x
  \param y   the shift of the target along y
  \param z_1 z position of the cell
  \param z_2 z position of the cell
  \param R    the radial cut
  \param yCUT the cut for upper part of the cell
*/
static bool GetTargetLocation(int run,
                       double &x, double &y, double &z_1, double &z_2,
                       double &R, double &yCUT);


/*! \brief The requirement that the muon beam trajectory crosses
  entirely the target cell. It is used in order to equalize fluxes
  through both cells.

  \param par    the beam track parameters in the primary vertex
  \param run    the run number ("-1" = MC)
  \param R_U    the user defined radial cut (R<1.5cm)
  \param yCUT_U the user defined vertical cut (y<1cm)
*/
 static bool CrossCells( PaTPar par, int run);



 /*! \brief The check for the primary vertex to be in one of the target cells.
   \param par    the beam track parameters in the primary vertex
   \param run    the run number ("-1" = MC)
   \param R_U    the user defined radial cut (R<1.4cm)
   \param yCUT_U the user defined vertical cut (y<1cm)
 */
 static bool InTarget( PaTPar par, int run);

};

#endif
