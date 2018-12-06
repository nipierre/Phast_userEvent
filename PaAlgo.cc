#include <iostream>
#include <iomanip>

#include "TSpline.h"
#include "TRandom.h"

#include "Phast.h"
#include "PaAlgo.h"

vector<double> PaAlgo::xv = std::vector<double>();	// initiaize empty vectors
vector<double> PaAlgo::yv = std::vector<double>();	// (ctarget cell centres)
vector<double> PaAlgo::zv = std::vector<double>();
double PaAlgo::xMC = 0;
double PaAlgo::phiMC = 0;
double PaAlgo::yMC = 0;
double PaAlgo::thetaMC = 0;
double PaAlgo::zMC = 0;

/* \brief Gives the target location in space: shift and tilting.
  Returns false if no information for the given year.
  Can be used for two cells configuration of the target
  (up to the end of 2004 or run 44108).

  \param run the run number ("-2" = 2 Cell MC)
  \param xU   the shift of the target most upstream part along x
  \param yU   the shift of the target most upstream part along y
  \param zU_1 z position of the upstream cell most upstream part
  \param zU_2 z position of the upstream cell most downstream part
  \param xD   the shift of the target most downstream part along x
  \param yD   the shift of the target most downstream part along y
  \param zD_1 z position of the downstream cell most upstream part
  \param zD_2 z position of the downstream cell most downstream part
  \param R    the radial cut
  \param yCUT the cut for upper part of the cell
  \author Alexandre.Korzenev@cern.ch
*/
bool PaAlgo::GetTargetLocation(int run,
		       double &xU, double &yU, double &zU_1, double &zU_2,
		       double &xD, double &yD, double &zD_1, double &zD_2,
		       double &R, double &yCUT)
{
  if( run > 0 && run < 10000 ) {
    cout<<"PaAlgo::GetTargetLocation ==> target location for the run "<<run
	<<" is not known. Most probably you entered the year instead of the run number."<<endl;
    exit(1);
  }

  if( run == -1 || run == -2 ) {  // for MC
    xU =    0; yU = 0;      zU_1 = -100; zU_2 = -40;
    xD =    0; yD = 0;      zD_1 = -30 ; zD_2 =  30;
    R    = 1.5;
    yCUT =   R;
  } else if( run < 15992 ) { // 2001 and before
    cout<<"PaAlgo::GetTargetLocation ==> target location for the run "<<run
	<<" is not known. The first run of 2002 is 15992."<<endl;
    exit(1);
  } else if( run < 24277 ) { // 2002
    xU = -0.2; yU = 0.1;    zU_1 = -100; zU_2 = -40;
    xD = -0.3; yD = -0.15;  zD_1 = -30 ; zD_2 =  30;
    R    = 1.4;
    yCUT = 1;
  } else if( run < 32833 ) { // 2003
    xU =  0.0; yU = -0.1;   zU_1 = -100; zU_2 = -40;
    xD = -0.1; yD = -0.25;  zD_1 = -30 ; zD_2 =  30;
    R    = 1.4;
    yCUT = 1;
  } else if( run < 44108 ) { // 2004
    // averaged over Marcin & Sebastien values
    xU =  0.0; yU = -0.1;   zU_1 = -100; zU_2 = -40;
    xD =  0.0; yD = -0.3;   zD_1 = -30 ; zD_2 =  30;
    R    = 1.4;
    yCUT = 1;
    if( run > 40000 ) { // after target reinforcement 06 September 2004.
      yU += -0.03;
      yD +=  0.1;
    }
  } else if( run > 253100 && run < 255230 ) { // DY2014
    xU =  0.0; yU = 0.0;   zU_1 = -294.6; zU_2 = -239.7;
    xD =  0.0; yD = 0.0;   zD_1 = -220.0; zD_2 = -164.3;
    R    = 2.0;
    yCUT = R;
  } else if( run > 256000 && run < 264860 ) { // DY2015
    xU =  0.0; yU = 0.0;   zU_1 = -294.5; zU_2 = -239.3;
    xD =  0.0; yD = 0.0;   zD_1 = -219.5; zD_2 = -164.3;
    R    = 2.0;
    yCUT = R;
  } else { // 2005 and after
    cout<<"PaAlgo::GetTargetLocation ==> target location for the run "<<run
	<<" is not known. The last run of 2004 is 43425."<<endl;
    exit(1);
  }

  // from transversity group ---- start

  if(         run > 21177 && run <  21880 ) { // 2002 transversity - p2b/p2c
    xU = -0.23; yU =  0.06;  zU_1 = -100; zU_2 = -40;
    xD = -0.31; yD = -0.14;  zD_1 = -30 ; zD_2 =  30;
    R    = 1.3;
    yCUT = 1.4;
  } else if ( run > 23480 && run < 23840 ) { // 2002 transversity - p2h
    xU = -0.23; yU =  0.06;  zU_1 = -100; zU_2 = -40;
    xD = -0.31; yD = -0.14;  zD_1 = -30 ; zD_2 =  30;
    R    = 1.3;
    yCUT = 1.4;
  } else if ( run > 30770 && run < 31530 ) { // 2003 transversity - p1g/p1h
    xU = -0.02; yU = -0.01;  zU_1 = -100; zU_2 = -40;
    xD = -0.09; yD = -0.19;  zD_1 = -30 ; zD_2 =  30;
    R    = 1.3;
    yCUT = 1.4;
  } else if ( run > 38990 && run < 39547 ) { // 2004 transversity - w33/w34
    xU =  0.03; yU = -0.09;  zU_1 = -100; zU_2 = -40;
    xD =  0.01; yD = -0.35;  zD_1 = -30 ; zD_2 =  30;
    R    = 1.3;
    yCUT = 1.4;
  } else if ( run > 39546 && run < 39990 ) { // 2004 transversity - w34/w35
    xU =  0.03; yU = -0.09;  zU_1 = -100; zU_2 = -40;
    xD =  0.01; yD = -0.35;  zD_1 = -30 ; zD_2 =  30;
    R    = 1.3;
    yCUT = 1.4;
  }

  // from transversity group ---- end

  return true;
}





/* \brief Gives the target location in space: shift and tilting.
  Returns false if no information for the given year.
  Can be used for three cells configuration of the target
  (from the beginning of year 2006).

  \param run the run number ("-3" = 3 Cell MC LiD, "-4" = 3 Cell MC NH3)
  \param xU    shift in x of the most upstream part of the target
  \param yU    shift in y of the most upstream part of the target
  \param zU_1  z position of the most upstream part of the upstream cell
  \param zU_2  z position of the most downstream part of the upstream cell
  \param zC_1  z position of the most upstream part of the central cell
  \param zC_2  z position of the most downstream part of the central cell
  \param xD    shift in x of the most downstream part of the target
  \param yD    shift in y of the most downstream part of the target
  \param zD_1  z position of the most upstream part of the downstream cell
  \param zD_2  z position of the most downstream part of the downstream cell
  \param R     the radial cut
  \param yCUT  the cut for upper part of the cell
  \author Alexandre.Korzenev@cern.ch
*/
//Addition of the 2011 target location determined by Artem (Vincent)
bool PaAlgo::GetTargetLocation(int run,
		       double &xU, double &yU, double &zU_1, double &zU_2,
		         double &zC_1, double &zC_2,
		       double &xD, double &yD, double &zD_1, double &zD_2,
		       double &R, double &yCUT)
{
  if( run > 0 && run < 45000 ) {
    cout<<"PaAlgo::GetTargetLocation ==> target location for the run "<<run
	<<" is not known. For runs < 45000 use 2-cells configuration."<<endl;
    exit(1);
  }

  if( run == -3 ) {  // for MC LiD
    double dz = 2;
    xU = 0; yU = 0;      zU_1 = -65+dz; zU_2 = -35+dz;
                         zC_1 = -30+dz; zC_2 =  30+dz;
    xD = 0; yD = 0;      zD_1 =  35+dz; zD_2 =  65+dz;
    R    = 1.4;
    yCUT = R;
  } else if( run == -4 ) {  // for MC NH3
    double dz = 2;
    xU = 0; yU = 0;      zU_1 = -65+dz; zU_2 = -35+dz;
                         zC_1 = -30+dz; zC_2 =  30+dz;
    xD = 0; yD = 0;      zD_1 =  35+dz; zD_2 =  65+dz;
    R    = 1.9;
    yCUT = R;
  } else if( run < 55456 ) { // 2006
    double dz = 2;
    xU = -0.1;  yU = 0.33;      zU_1 = -65+dz+4; zU_2 = -35+dz;
                                zC_1 = -30+dz+8; zC_2 =  30+dz;
    xD = -0.07; yD = 0.33;      zD_1 =  35+dz+2; zD_2 =  65+dz;
    R    = 1.4;
    yCUT = R;
  } else if( run < 64435 ) { // 2007
    if( (run>=57058 && run<=60332) || (run>=62650 && run<=63816) ) { // transv
      xU = -0.15;  yU = 0.39;      zU_1 = -62.5; zU_2 = -32.5;
                                   zC_1 = -27.5; zC_2 =  32.5;
      xD =  0.02;  yD = 0.23;      zD_1 =  37.5; zD_2 =  67.5;
      R    = 1.9;
      yCUT = R;
    } else { // longit
      xU = -0.15;  yU = 0.39;      zU_1 = -62.5+2; zU_2 = -32.5;
                                   zC_1 = -27.5+2; zC_2 =  32.5;
      xD =  0.02;  yD = 0.23;      zD_1 =  37.5+0; zD_2 =  67.5;
      R    = 1.9;
      yCUT = R;
    }
  }  else  if( run>=82363 && run<= 89392) { //2010, transv
      xU = -0.2;  yU = 0.02;      zU_1 = -62.5; zU_2 = -32.5;
                                  zC_1 = -27.5; zC_2 =  32.5;
      xD = -0.2;  yD = 0.02;      zD_1 =  37.5; zD_2 =  67.5;
      R    = 1.9;
      yCUT = R;
  }   else if(run>=89393 && run<= 96223) {//2011, longit
      double dz = -1.2;
      xU = -0.0073;  yU = -0.1680;       zU_1 = -61.7+dz; zU_2 = -32.4+dz;
                                         zC_1 = -26.5; zC_2 =  32.1;
      xD =  0.1261;  yD = -0.1796;       zD_1 =  37.2; zD_2 =  66.5;
    R  = 1.9;
    yCUT = R;
  } else { // 2008 and 2009
     cout<<"PaAlgo::GetTargetLocation ==> target location for the run "<<run
	 <<" is not known. "<<endl;
     exit(1);
   }

  return true;
}

/* \brief ives the x and y coordinates of the target center for a given z.
  Returns false if no information for the given year.
  Can be used for one cell configuration of the target
  (from the years 2012/2016/2017).

  \param run the run number
  \param xC    x(z) of the centre of the target
  \param yC    y(z) of the centre of the target
  \param z     z position in the target (input parameter)
  \param R     the recommended radial cut
  \param yCUT  the recommended 'hydrogen level cut' (y < yCUT)
  \author antoine.vidon@cern.ch, nicolas.pierre@cern.ch, karolina.juraskova@cern.ch, jan.matousek@cern.ch
*/

bool PaAlgo::GetTargetLocation(int run, double &xC, double &yC, double &xCmc, double &yCmc, double z, double &R, double &RMC, double &yCUT) // !!NEW!!
{
  xC = 1000000;
  yC = 1000000;
 if( !(xv.size() && yv.size() && zv.size()) )  // Check if already initialized
 {
    std::ifstream fin, finmc;
    std::string tstr, tstrmc;

    if( 96224 <= run && run <= 109125 )
		{
			tstr = "/afs/cern.ch/compass/dvcs/Production/Analysis/Target/target-107924-109081.dat"; // 2012
			tstrmc = "/afs/cern.ch/compass/dvcs/Production/Analysis/Target/target-mc-2012.dat"; // 2012
		}
    else if( 264860 <= run && run <= 276879 )
		{
			tstr = //"/afs/cern.ch/compass/dvcs/Production/Analysis/Target/target-274508-274901-1.dat"; // 2016
			"/eos/user/j/jmatouse/analysis/Sidis-2016/realdata/TargetCell_cut/target-274508-274901-2.dat";
			tstrmc = "/afs/cern.ch/compass/dvcs/Production/Analysis/Target/target-mc-2016.dat"; // 2016
		}
    else if( 276880 <= run && run <= 281775 )
		{
			tstr = "/afs/cern.ch/compass/dvcs/Production/Analysis/Target/target-278473-278706-0.dat"; // 2017
			tstrmc = "/afs/cern.ch/compass/dvcs/Production/Analysis/Target/target-mc-2017.dat"; // 2017
		}
    else return false; //check, otherwise segmentation fault
    fin.open(tstr.c_str());
    while(fin.is_open() && !fin.eof()) {
      double z, x, y, dummy;
      fin >> z >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy
          >> x >> y;
      zv.push_back(z);
      xv.push_back(x);
      yv.push_back(y);
      if (z >=0) //check if file is empty or with wrong numbers
			{
      	cout << "!!! PaAlgo::GetTargetLocation() PROBLEM: EMPTY or WRONG TARGET GEOMETRY FILE !!!" << endl;
      	return false;
    	}
    }
		fin.close();
		finmc.open(tstrmc.c_str());
		while(finmc.is_open() && !finmc.eof()) {
      fin >> xMC >> phiMC >> yMC >> thetaMC >> zMC;
      if (zMC >=0) //check if file is empty or with wrong numbers
			{
      	cout << "!!! PaAlgo::GetTargetLocation() PROBLEM: EMPTY or WRONG MC TARGET GEOMETRY FILE !!!" << endl;
      	return false;
    	}
    }
		finmc.close();
  }

  if( !(xv.size() && yv.size() && zv.size()) ) {
  	cout << "!!! PaAlgo::GetTargetLocation() PROBLEM: NONEXISTING TARGET GEOMETRY FILE !!!" << endl;
  	return false; //check, otherwise segmentation fault if there is nothing in vectors xv, yv, zv
	}

  R=1.9;    //will be set if R is not defined by user in CrossCells or in InTarget
	RMC=2;    //will be set if R is not defined by user in CrossCells or in InTarget
	yCUT=1.2; //will be set if yCUT is not defined by user in CrossCells or in InTarget

  for(unsigned int i = 0; i < zv.size()-1; i++ ) {
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
    xC = xc1 + dxcdz*dz;
    yC = yc1 + dycdz*dz;
    break;
  }

	xCmc = xMC+sin(phiMC)*(zMC-z);
	yCmc = yMC+sin(thetaMC)*(zMC-z);

  return true;
}

/* \brief The requirement that the muon beam trajectory crosses
  entirely two target cells. It is used in order to equalize fluxes
  through both cells.

  \param par    the beam track parameters in the primary vertex
  \param run    the run number ("-2" = 2 Cell MC, "-3" = 3 Cell MC LiD, "-4" = 3 Cell MC NH3)
  \param R_U    the user defined radial cut (R<1.4cm), for 2012/2016/2017: use (R<1.9cm), if R_U is not set by user it is set according to GetTargetLocation.
  \param yCUT_U the user defined vertical cut (y<1cm), for 2012/2016/2017: use (y<1.2cm)), if yCUT_U is not set by user it is set according to GetTargetLocation.
  \param zmin_U the user defined zmin of the target - now available only for 2012/2016/2017. If zmin_U is not set by user it is set according to GetTargetLocation.
  \param zmax_U the user defined zmax of the target - now available only for 2012/2016/2017. If zmax_U is not set by user it is set according to GetTargetLocation.
*/
bool PaAlgo::CrossCells( PaTPar par, int run, double R_U, double yCUT_U, double zmin_U, double zmax_U, double RMC_U ) // Added zmin/zmax in case some people want to have stricter cuts than in target file
{
  PaTPar parE;
  double xU,yU,zU_1,zU_2, xD,yD,zD_1,zD_2, R, RMC, yCUT, xC, yC, yCmc, yCmc, zmin, zmax; // !!NEW!! xC, yC, zmin, zmax declaration for 2012/2016/2017

  // !!NEW!!
  if( (96224 <= run && run <= 109125) || (264860<= run && run <= 281775) ) // 1 cell   2012/2016/2017
  {
    double z = par.Z();  //function GetTargetLocation(run,xC,yC,z) is called latter with different z arguments, here it is just check
    if( !GetTargetLocation(run,xC,yC,xCmc,yCmc,z, R, RMC, yCUT) ) {
      cout<<"PaAlgo::CrossCells PROBLEM: no info for the run "<<run<<" (1 cell)"<<endl;
      return false;
    }
  }
  // !!NEW!!
  else if( run < 45000 && run>=-2) {      // 2 cells
    if( !GetTargetLocation(run, xU,yU,zU_1,zU_2, xD,yD,zD_1,zD_2, R,yCUT) ) {
      cout<<"PaAlgo::CrossCells PROBLEM: no info for the run "<<run<<" (2 cells)"<<endl;
      return false;
    }
  } else {                 // 3 cells
    double zC_1, zC_2; // dummy
    if( !GetTargetLocation(run, xU,yU,zU_1,zU_2,zC_1,zC_2, xD,yD,zD_1,zD_2, R,yCUT) ) {
      cout<<"PaAlgo::CrossCells PROBLEM: no info for the run "<<run<<" (3 cells)"<<endl;
      return false;
    }
  }

  if( R_U    != -9999 ) R    = R_U;
	if( RMC_U  != -9999 ) RMC  = RMC_U;
  if( yCUT_U != -9999 ) yCUT = yCUT_U;

  // !!NEW!!
  if( (96224 <= run && run <= 109125) || (264860 <= run && run <= 281775) )  //2012/2016/2017
  {
    if( zmin_U == -9999) {
      zmin = zv[0];   // if zmin_U is not set by user it is set according to the target file in GetTargetLocation
      cout <<"PaAlgo::CrossCells WARNING: zmin was not defined by user and it was set automatically to zmin= "<< zmin << "cm" << endl;
    }
    else zmin = zmin_U;

    if( zmax_U == -9999) {
    zmax = zv[zv.size()-1];  // if zmax_U is not set by user it is set according to the target file in GetTargetLocation
    cout <<"PaAlgo::CrossCells WARNING: zmax was not defined by user and it was set automatically to zmax= "<< zmax << "cm" << endl;
    }
    else zmax = zmax_U;

    for( unsigned int i = 1; i < zv.size(); i++ )
    {
      double z0 = zv[i-1];
      double z1 = zv[i];
      if( z1 <= zmin ) continue;
      if( z0 >= zmax ) continue;

      if( zmin > z0 && zmin < z1) {
        // Zmin is inside the current segment, so we need to interpolate
        par.Extrapolate(zmin, parE, false);
        if (!InTarget(parE,'O',run,R,yCUT,zmin,zmax,RMC)) return false;
        }
      if( zmax > z0 && zmax < z1) {
        // Zmax is inside the current segment, so we need to interpolate
        par.Extrapolate(zmax, parE, false);
        if (!InTarget(parE,'O',run,R,yCUT,zmin,zmax,RMC)) return false;
		}
	  if (zmax > z1) {
        // Extrapolate to the downstream end of the current segment
        par.Extrapolate(z1, parE, false);
        if (!InTarget(parE,'O',run,R,yCUT,zmin,zmax,RMC)) return false;
	  }
    }
  }
  // !!NEW!!
  else
  {
    par.Extrapolate(zU_1,parE,0);
    double xL = parE(1);
    double yL = parE(2);

    par.Extrapolate(zD_2,parE,0);
    double xR = parE(1);
    double yR = parE(2);

    double rL = sqrt( (xL-xU)*(xL-xU) + (yL-yU)*(yL-yU) );
    if( rL > R ) return false;
    if( yL-yU > yCUT ) return false;

    double rR = sqrt( (xR-xD)*(xR-xD) + (yR-yD)*(yR-yD) );
    if( rR > R ) return false;
    if( yR-yD > yCUT ) return false;
  }

  return true;

}



/* \brief The check for the primary vertex to be in one of the target cells.
  \param par    the beam track parameters in the primary vertex
  \param Cell   the one of cells (if 'U' - upstream cell, if 'D' - downstream, if 'C' - center, if 'O' - one cell)
  \param run    the run number ("-2" = 2 Cell MC, "-3" = 3 Cell MC LiD, "-4" = 3 Cell MC NH3)
  \param R_U    the user defined radial cut (R<1.4cm), for 2012/2016/2017: use (R<1.9cm), if R_U is not set by user it is set according to the target file in GetTargetLocation.
  \param yCUT_U the user defined vertical cut (y<1cm), for 2012/2016/2017: use (y<1.2cm)), if yCUT_U is not set by user it is set according to the target file in GetTargetLocation.
  \param zmin_U the user defined zmin of the target - now available only for 2012/2016/2017. If zmin_U is not set by user it is set according to the target file in GetTargetLocation.
  \param zmax_U the user defined zmax of the target - now available only for 2012/2016/2017. If zmax_U is not set by user it is set according to the target file in GetTargetLocation.
*/
bool PaAlgo::InTarget( double x, double y, double z, char Cell, int run, double R_U, double yCUT_U, double zmin_U, double zmax_U, double RMC_U ) // Added zmin/zmax in case some people want to have stricter cuts than in target file
{
  double xU,yU, zU_1, zU_2, zC_1, zC_2, xD, yD, zD_1, zD_2, R, RMC, yCUT, xC, yC, xCmc, yCmc, r, rMC;
	// xC, yC, r defined here because it is used bellow for different numberb of cells
  double zmin = 0;	// !!NEW!! zmin, zmax declaration for 2012/2016/2017
  double zmax = 0;

  // !!NEW!!
  if( (96224 <= run && run <= 109125) || (264860 <= run && run <= 281775) ) // DVCS 2012/2016/2017
  {
    if( Cell == 'U' || Cell == 'C' || Cell == 'D' )
    {
      cout<<"PaAlgo::InTarget() PROBLEM: target for the run "<<run<<" has 1 cell!"<<endl;
      return false;
    }
    if( !GetTargetLocation(run,xC,yC,xCmc,yCmc,z, R, RMC, yCUT) ) {
      cout<<"PaAlgo::InTarget() PROBLEM: no info for the run "<<run<<" (1 cell)"<<endl;
      return false;
    }
  }
  // !!NEW!!
  else if( (run < 45000 && run>=-2) || (255235 <= run && run <= 264859) || (281776 <= run && run <= 287559) ) // 2 cells
  {
    if( Cell == 'C' || Cell == 'O' ) {
      cout<<"PaAlgo::InTarget() PROBLEM: target for the run "<<run<<" has 2 cells!"<<endl;
      return false;
    }
    if( (run < 45000 && run>=-2) )
    {
      if( !GetTargetLocation(run, xU,yU,zU_1,zU_2, xD,yD,zD_1,zD_2, R,yCUT) ) {
        cout<<"PaAlgo::InTarget() PROBLEM: no info for the run "<<run<<" (2 cells)"<<endl;
        return false;
      }
    }
    else // DY 2015/2018 TBD
    {

    }
  }
  else // 3 cells
  {
    if( Cell == 'O' ) {
      cout<<"PaAlgo::InTarget() PROBLEM: target for the run "<<run<<" has 3 cells!"<<endl;
      return false;
    }
    if( !GetTargetLocation(run, xU,yU,zU_1,zU_2,zC_1,zC_2, xD,yD,zD_1,zD_2, R,yCUT) ) {
      cout<<"PaAlgo::InTarget() PROBLEM: no info for the run "<<run<<" (3 cells)"<<endl;
      return false;
    }
  }

  if( R_U    != -9999 ) R    = R_U;
	if( RMC_U  != -9999 ) RMC  = RMC_U;
  if( yCUT_U != -9999 ) yCUT = yCUT_U;

  if(Cell != 'O')
  {
    xC = (xD-xU) * (zU_1-z) / (zU_1-zD_2) + xU;
    yC = (yD-yU) * (zU_1-z) / (zU_1-zD_2) + yU;
    r = sqrt( (x-xC)*(x-xC) + (y-yC)*(y-yC) );
  }
  // !!NEW!!
  else  //2012/2016/2017
  {
    r = sqrt( (x-xC)*(x-xC) + (y-yC)*(y-yC) );
		rMC = sqrt( (x-xCmc)*(x-xCmc) + (y-yCmc)*(y-yCmc) );

		if( R_U    == -9999 ) cout <<"PaAlgo::InTarget WARNING: R_cut was not defined by user and it was set automatically to R= "<< R << "cm" << endl; //set according to the target file in GetTargetLocation
		if( RMC_U    == -9999 ) cout <<"PaAlgo::InTarget WARNING: RMC_cut was not defined by user and it was set automatically to R= "<< RMC << "cm" << endl; //set according to the MC target file in GetTargetLocation
    if( yCUT_U == -9999 ) cout <<"PaAlgo::InTarget WARNING: Y_cut was not defined by user and it was set automatically to Y= "<< yCUT << "cm" << endl;  //set according to the target file in GetTargetLocation

    if( zmin_U == -9999) {
      zmin = zv[0];   // if zmin_U is not set by user it is set according to the target file in GetTargetLocation
      cout <<"PaAlgo::InTarget WARNING: zmin was not defined by user and it was set automatically to zmin= "<< zmin << "cm" << endl;
    }
    else zmin = zmin_U;

    if( zmax_U == -9999) {
    zmax = zv[zv.size()-1];  // if zmax_U is not set by user it is set according to the target file in GetTargetLocation
    cout <<"PaAlgo::InTarget WARNING: zmax was not defined by user and it was set automatically to zmax= "<< zmax << "cm" << endl;
    }
    else zmax = zmax_U;

    if( zmin < zv[0] || zmin > 0) {
      zmin = zv[0];
      cout << "PaAlgo::InTarget WARNING: user defined zmin out of range, zmin was set according to the target file" << endl;
      cout << "If you need longer target, create different target file and call it instead of defautl one in GetTargetLocation() for 2012/2016/2017 " << endl;
    }

    if( zmax > zv[zv.size()-1] ) {
      zmax = zv[zv.size()-1];
      cout << "PaAlgo::InTarget WARNING: user defined zmax out of range, zmax was set according to the target file" << endl;
      cout << "If you need longer target, create different target file and call it instead of defautl one in GetTargetLocation() for 2012/2016/2017 " << endl;
    }
  }
  // !!NEW!!

  if(        Cell == 'U' ) {
    if( z < zU_1 || z > zU_2 ) return false;
  } else if( Cell == 'C' ) {
    if( z < zC_1 || z > zC_2 ) return false;
  } else if( Cell == 'D' ) {
    if( z < zD_1 || z > zD_2 ) return false;
  } else if( Cell == 'O' ) {
    if( z < zmin || z > zmax ) return false; // !!NEW!!
  } else {
    cout<<"inTarget PROBLEM: no info for cell "<<Cell<<endl;
    return false;
  }

  if(         r > R        ) return false;
  if(Cell != 'O')
  {
    if(      y-yC > yCUT     ) return false;
  }
  else
  {
    if(      y > yCUT || rMC > RMC     ) return false; // !!NEW!!
  }

  return true;
}


// by Alexandre.Korzenev@cern.ch
// mods for performance improvement by  Roland Kuhn <rkuhn@e18.physik.tu-muenchen.de>

// Addition of the 2011 beam polarisation and correction of the sign of the polarisation (Vincent)
// Correction of the range of the TSpline evaluation for year < 2011 (Vincent)

// Improvments for 2011 beam polarization
// (mail of Ana Sofia Nunes 20 April 2015)

double PaAlgo::GetBeamPol( float mom, int year )
{
  static const int n = 14;
//  static const int n2011 = 12;
  static const int n2011 = 25;

  static double x_03[n]     = {  120,  130,  145,  152,  155,   157,   159,   161,   163,   165,   168,  175,  185,   200};
  static double x_04[n]     = {  120,  135,  145,  152,  155,   157,   159,   161,   163,   165,   168,  175,  185,   200};
//  static double x_05[n2011] = {  185,  191,  193,  195,  197,   199,   201,   203,   205,   207,   209,  215};
  static double x_05[n2011] = {165,172.5,177.5,181,183,185,187,189,191,193,195,197,199,201,203,205,207,209,211,213,215,217,219,222.5,227.5};

  static double y_03[n]     = {-0.55 ,-0.57 ,-0.582,-0.630,-0.684,-0.723 ,-0.756 ,-0.782 ,-0.803 ,-0.823 ,-0.844,-0.855,-0.8551,-0.8552};
  static double y_04[n]     = {-0.574,-0.575,-0.61 ,-0.70 ,-0.757,-0.7836,-0.8040,-0.8251,-0.8443,-0.8588,-0.8708,-0.88 ,-0.883,-0.8835};
//  static double y_05[n2011] = {-0.6491, -0.7143, -0.7471, -0.7720, -0.7936, -0.8063, -0.8273, -0.8438, -0.8601, -0.8722, -0.8802, -0.8862}; //2011 beam good muon polarisation provided by Lau, statistical error
  // well below 1 permille and the systematics are estimated to about 3%
  static double y_05[n2011] = {-0.6703060915,-0.6534957096,-0.6313457782,-0.6302088022,-0.627365988,-0.6325876395,-0.658481391,-0.6920352993,-0.7268828016,-0.7577037021,-0.781735445,-0.800817675,-0.819343704,-0.8296627879,-0.8426811817,-0.855950396,-0.8697056884,-0.8818718939,-0.8896650166,-0.8911495472,-0.8882634868,-0.8856977802,-0.8890861625,-0.9070184047,-0.9400768204};//2011 beam good muon polarisation provided by Lau, statistical error

  static TSpline3 *spl03 = NULL; static TSpline3 *spl04 = NULL; static TSpline3 *spl05 = NULL;
  if (spl03 == NULL) {
    spl03 = new TSpline3("spl03", x_03, y_03, n);
    spl04 = new TSpline3("spl04", x_04, y_04, n);
    spl05 = new TSpline3("spl05", x_05, y_05, n2011);
  }

  if( year == 2002 || year == 2003 ) {
    if( mom <= x_03[0] )   return y_03[0];// put 0 index instead of 1 (Vincent)
    if( mom >= x_03[n-1] ) return y_03[n-1];
    return spl03->Eval(mom);
  } else if( year == 2004 || year == 2006 || year == 2007 ) {
    if( mom <= x_04[0] )   return y_04[0];// put 0 index instead of 1 (Vincent)
    if( mom >= x_04[n-1] ) return y_04[n-1];
    return spl04->Eval(mom);
  } else if( year == 2011 ) {
    if( mom <= x_05[0] )   return y_05[0];
    if( mom >= x_05[n2011-1] ) return y_05[n2011-1];
    return spl05->Eval(mom);
  } else {
    cout<<"PaAlgo::GetBeamPol() ==> ERROR: No parametrization for year "<<year<<endl;
    return 0;
  }
}



// the code in this function is based on the A1 analysis code which was made
// available by A. Korzenev
void PaAlgo::GetDepolarizationFactor(double q2, double xBj, double y,
                                     double R, double dR,
                                     double&D, double&dD, bool do_err)
{
  const double M = G3partMass[14]; // proton  0.938272
  const double m = G3partMass[5];  // muon    0.105658357;

  Double_t g2 = 4*M*M*xBj*xBj/q2;

  D = y * ( (1+g2*y/2)*(2-y) - 2*y*y*m*m/q2 ) /
    ( y*y* (1-2*m*m/q2) * (1+g2) + 2*(1+R) * (1-y-g2*y*y/4) );

  if (do_err) {
    double a = y * ( (1+g2*y/2)*(2-y) - 2*y*y*m*m/q2 );
    double b = y*y* (1-2*m*m/q2) * (1+g2) + 2 * (1-y-g2*y*y/4);
    double c = 2 * (1-y-g2*y*y/4);

    dD = dR * a*c / pow( (b+c*R), 2 );
  } else {
    dD = -1;
  }
}

#ifndef NO_FORTRAN

typedef struct rcpoin {
	float hydrc;
	float herc;
	float deutrc;
	float carbrc;
	float oxrc;
	float curc;
	float nirc;
	float deutjm;
	float hydjm;
	float nkrc;
} rcpoin;

extern "C" {

  void dfdeut_(float*,float*,int*,float*,float*);
  void df2007_(float*,float*,int*,float*,float*);
  void df2011_(float*,float*,int*,float*,float*);
  void r1990_(float*,float*,float*,float*);
  void r1998_(float*,float*,float*,float*);
  rcpoin rcpoin_;
  void rcinlx_(float*, float*, float*);
  void rcinhd_(float*, float*, float*);
  void rcinlx11_(float*, float*, float*);
  void rcinhd11_(float*, float*, float*);
}

// contributed by Alexander.Korzenev@cern.ch

// Dilution factor for 2010 data computed as if it was 2007 data, df2010.F has to be created and fill
// with the values available at http://compass02.cern.ch/elog/target_polar/127
// Dilution factor for 2011 data included and computed in the class df2011.F (Vincent)
void PaAlgo::GetDilutionFactor(float xBj, float y, char Cell, int run, int flag,
                               float &f, float &df)
{
  float dftot[12], err_dftot[12];

  if( Cell != 'U' && Cell != 'C' && Cell != 'D' ) {
    cout<<"PaAlgo::GetDilutionFactor() PROBLEM: Unclear cell = "<<Cell<<endl;
    exit(0);
  }

  if( flag != 1 && flag != 2 ) {
    cout<<"PaAlgo::GetDilutionFactor() PROBLEM: Unclear event type = "
	<<flag<<". Should be 1 for inclusive events or 2 for hadronic events."<<endl;
    exit(0);
  }

  if( run < 45000 ) {      // 2 cells
    if( Cell == 'C' ) {
      cout<<"PaAlgo::GetDilutionFactor() PROBLEM: target for the run "<<run<<" has 2 cells!"<<endl;
      exit(0);
    }
    if(      flag == 1 ) flag = 3; // inclusive events
    else if( flag == 2 ) flag = 4; // hadronic events
    dfdeut_(&xBj, &y, &flag, dftot, err_dftot);
    if( Cell == 'U' ) { f = dftot[1]; df = err_dftot[1]; }
    else              { f = dftot[2]; df = err_dftot[2]; }

  } else {                 // 3 cells

    if( run < 55456 ) { // deuteron target, 2006
      if(      flag == 1 ) flag = 5; // inclusive events
      else if( flag == 2 ) flag = 6; // hadronic events
      dfdeut_(&xBj, &y, &flag, dftot, err_dftot);
      if     ( Cell == 'U' ) { f = dftot[1]; df = err_dftot[1]; }
      else if( Cell == 'C' ) { f = dftot[2]; df = err_dftot[2]; }
      else                   { f = dftot[3]; df = err_dftot[3]; }

    } else if( run < 64434 || ( run > 82362 && run < 89393))  { // proton target, 2007 = run<64434 or 2010 run > 82362 && run < 89393
      df2007_(&xBj, &y, &flag, dftot, err_dftot);
      if     ( Cell == 'U' ) { f = dftot[1]; df = err_dftot[1]; }
      else if( Cell == 'C' ) { f = dftot[2]; df = err_dftot[2]; }
      else                   { f = dftot[3]; df = err_dftot[3]; }
      //printf("2007\n");

    } else if( run > 89392 && run < 96224){ // proton target, 2011 at 200 GeV
      df2011_(&xBj, &y, &flag, dftot, err_dftot);
      if     ( Cell == 'U' ) { f = dftot[1]; df = err_dftot[1]; }
      else if( Cell == 'C' ) { f = dftot[2]; df = err_dftot[2]; }
      else                   { f = dftot[3]; df = err_dftot[3]; }
    } else { //2008, 2009, 2012
      cout<<"PaAlgo::GetDilutionFactor ==> Dilution factor for the run "<<run<<" is not known. "<<endl;
      f = 0.;
      df = 0.;
    }
  } // if 2 or 3 cells
}



// contributed by Roland Kuhn
// <roland.kuhn@cern.ch>
void PaAlgo::GetDilutionFactor(float xBj, float y, int flag, float*dftot,
                               float*err_dftot)
{
  //float flag_ = flag;
  dfdeut_(&xBj, &y, &flag, dftot, err_dftot);
}

void PaAlgo::GetR(float q2, float xBj, float&R, float&dR, bool do_err)
{
  if (!do_err) {
    dR = -1;
  } else {
    dR = 0;
  }
  // r1990_(&xBj, &q2, &R, &dR);
  r1998_(&xBj, &q2, &R, &dR);
}

void PaAlgo::GetDepolarizationFactor(double q2, double xBj, double y,
                                     double&D, double&dD, bool do_err)
{
  float R, dR;
  PaAlgo::GetR(q2, xBj, R, dR, do_err);
  PaAlgo::GetDepolarizationFactor(q2, xBj, y, R, dR, D, dD, do_err);
}

// contributed by Konrad Klimaszewski
// <Konrad.Klimaszewski@cern.ch>
// Take into account the year of the data taking to use the correct table for 2011 data at 200 GeV.
float PaAlgo::GetRadiativeWeight(float xBj, float y, int flag)
{
  if( flag != 1 && flag != 2 && flag != 11 && flag != 12 ) {
    cerr<<"PaAlgo::GetRadiativeWeight() PROBLEM: Unclear event type = "<<flag<<"."<<endl;
    cerr<<"Should be 1 for inclusive events or 2 for hadronic events for deuteron target."<<endl;
    cerr<<"Should be 11 for inclusive events or 12 for hadronic events for hydrogen target."<<endl;
    exit(0);
  }
  const PaEvent& e = Phast::Ref().event;//Get the year of the event to the use the tables corresponding to the right beam energy.
  // According to Barbara Badelek tables are calculated for 160 GeV energy. No point to put variable one here.
  float energy = 0;
  float rc;

  if(e.Year() == 2011){
    energy = 200.0;
    if(flag==11) // inclusive events
      rcinlx11_(&xBj, &y, &energy);
    else if (flag==12)// hadronic events
      rcinhd11_(&xBj, &y, &energy);
    else {
      cerr<<"PaAlgo::GetRadiativeWeight() PROBLEM: Unclear event type = "<<flag<<"."<<endl;
      cerr<<"Should be either 11 for inclusive or 12 for hadronic events for 2011 data"<<endl;
      exit(0);
    }
  }
  else {
    energy = 160.0;
    if(flag==1 || flag == 11) // inclusive events
      rcinlx_(&xBj, &y, &energy);
    else // hadronic events
      rcinhd_(&xBj, &y, &energy);
  }

  if(flag==1||flag==2) // deuteron target
      rc = 1.0/rcpoin_.deutrc;
  else // hydrogen target
      rc = 1.0/rcpoin_.hydrc;

  return rc;
}

void PaAlgo::GetRadiativeWeight(float xBj, float y, int flag, float& rc)
{
	cerr<<"The function 'void PaAlgo::GetRadiativeWeight(xBj, y, flag, rc)' is deprecated!"<<endl;
	cerr<<"Please use the new one: 'float PaAlgo::GetRadiativeWeight(xBj, y, flag)'."<<endl;
	cerr<<"Also please carefuly study the documentation as the meaning of returned value has changed."<<endl;
	exit(0);
}


#endif // NO_FORTRAN

// by Andrea Ferrero <aferrero@to.infn.it>
bool  PaAlgo::DoMW1ID(const PaTrack& t, int* nma1, int* nma2)
{
  int nma01 = t.NHitsFoundInDetect("MA01");
  int nma02 = t.NHitsFoundInDetect("MA02");
  if(nma1) *nma1 = nma01;
  if(nma2) *nma2 = nma02;
  if(nma01 < 4) return false;
  if(nma02 < 6) return false;
  return true;
}

// by Andrea Ferrero <aferrero@to.infn.it>
int  PaAlgo::GetMW1ScatMuon(const PaEvent& e)
{
  if((e.TrigMask()&0x3FFFFFF) != 16) return -1; // mask out onl filter bits

  int iprim = e.iBestPrimaryVertex();
  if(iprim == -1) return -1;

  const PaVertex& pv = e.vVertex(iprim);

  map<int,int> mup_id;
  for(int ti = 0; ti < pv.NOutParticles(); ++ti) {
    int pj = pv.iOutParticle(ti);
    const PaParticle& p = e.vParticle(pj);
    int tj = p.iTrack();
    if(tj < 0) continue;
    const PaTrack& t = e.vTrack(tj);
    int nma1, nma2;
    if(DoMW1ID(t,&nma1,&nma2)) {
      mup_id.insert(make_pair(nma1+10*nma2,tj));
    }
  }

  if(mup_id.size() == 0) return -1;
  return (*mup_id.rbegin()).second;
}
