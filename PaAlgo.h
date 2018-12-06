#ifndef PaAlgo_h
#define PaAlgo_h


/*!
  \class PaAlgo
  \brief Miscellaneous functions

  Collection of independent,"self-contained"
  algorithms and useful physics-related functions.

*/

#include <cmath>
#include <TVector3.h>
#include <TLorentzVector.h>
#include "PaEvent.h"
#include <vector>

using namespace std;


class PaAlgo

{

protected:
  static vector<double> xv, yv, zv;
  static double xMC, phiMC, yMC, thetaMC, zMC;

public:


/*! \brief Gives the target location in space: shift and tilting.
  Returns false if no information for the given year.

  \param run  the run number ("-2" = 2 cell MC)
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
static bool GetTargetLocation(int run,
                       double &xU, double &yU, double &zU_1, double &zU_2,
                       double &xD, double &yD, double &zD_1, double &zD_2,
                       double &R, double &yCUT);



/* \brief Gives the target location in space: shift and tilting.
  Returns false if no information for the given year.
  Can be used for three cells configuration of the target
  (from the beginning of year 2006).

  \param run the run number ("-3" = 3 cell MC LiD, "-4" = 3 cell MC NH3)
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
static bool GetTargetLocation(int run,
		       double &xU, double &yU, double &zU_1, double &zU_2,
		                               double &zC_1, double &zC_2,
		       double &xD, double &yD, double &zD_1, double &zD_2,
		       double &R, double &yCUT);



/*! \brief Gives the x and y coordinates of the target center for a given z.
  Returns false if no information for the given year.
  Can be used for one cell configuration of the target (2012/2016/2017).

  \param run the run number
  \param xC    x(z) of the centre of the target
  \param yC    y(z) of the centre of the target
  \param z     z position in the target (input parameter)
  \param R     the recommended radial cut
  \param yCUT  the recommended 'hydrogen level cut' (y < yCUT)
  \author karolina.juraskova@cern.ch, antoine.vidon@cern.ch, nicolas.pierre@cern.ch, jan.matousek@cern.ch
*/
static bool GetTargetLocation(int run,
           double &xC, double &yC, double &xCmc, double &yCmc, double z, double &R, double &RMC, double &yCUT); // !!NEW!!


/*! \brief The requirement that the muon beam trajectory crosses
  entirely two target cells. It is used in order to equalize fluxes
  through both cells.

  \param par    the beam track parameters in the primary vertex
  \param run    the run number ("-2" = 2 cell MC, "-3" = 3 cell MC LiD, "-4" = 3 cell MC NH3)
  \param R_U    the user defined radial cut (R<1.4cm), for 2012/2016/2017: use (R<1.9cm), if R_U is not set by user it is set according GetTargetLocation.
  \param yCUT_U the user defined vertical cut (y<1cm), for 2012/2016/2017: use (y<1.2cm)), if yCUT_U is not set by user it is set according to GetTargetLocation.
  \param zmin_U the user defined zmin of the target - only for 2012/2016/2017. If zmin_U is not set by user it is set according to GetTargetLocation.
  \param zmax_U the user defined zmax of the target - only for 2012/2016/2017. If zmax_U is not set by user it is set according to GetTargetLocation.
  \author Alexandre.Korzenev@cern.ch
  \author 2012/16/17 update: karolina.juraskova@cern.ch, antoine.vidon@cern.ch, nicolas.pierre@cern.ch, jan.matousek@cern.ch
*/
 static bool CrossCells( PaTPar par, int run, double R_U = -9999, double yCUT_U = -9999, double zmin_U = -9999, double zmax_U = -9999, double RMC_U = -9999 );



 /*! \brief The check for the primary vertex to be in one of the target cells.

 \param par  The beam track parameters in the primary vertex
 \param Cell The one of cells (if 'U' - upstream cell, if 'D' - downstream, if 'C' - central, 'O' - one cell)
 \param run  The run number ("-2" = 2 cell MC, "-3" = 3 cell MC LiD, "-4" = 3 cell MC NH3)
 \param R_U    the user defined radial cut (R<1.4cm), for 2012/2016/2017: use (R<1.9cm), if R_U is not set by user it is set according to the target file in GetTargetLocation.
 \param yCUT_U the user defined vertical cut (y<1cm), for 2012/2016/2017: use (y<1.2cm)), if yCUT_U is not set by user it is set according to the target file in GetTargetLocation.
  // !!NEW!!
 \param zmin_U the user defined zmin of the target - now available only for 2012/2016/2017. If zmin_U is not set by user it is set according to the target file in GetTargetLocation.
 \param zmax_U the user defined zmax of the target - now available only for 2012/2016/2017. If zmax_U is not set by user it is set according to the target file in GetTargetLocation.
 \author Alexandre.Korzenev@cern.ch
 \author 2012/16/17 update: karolina.juraskova@cern.ch, antoine.vidon@cern.ch, nicolas.pierre@cern.ch, jan.matousek@cern.ch
 */
 static bool InTarget( PaTPar par, char Cell, int run, double R_U = -9999, double yCUT_U = -9999, double zmin_U = -9999, double zmax_U = -9999, double RMC_U = -9999 )
 { return InTarget(par.X(), par.Y(), par.Z(), Cell, run, R_U, yCUT_U, zmin_U, zmax_U, RMC_U); }

 static bool InTarget( double X, double Y, double Z, char Cell, int run, double R_U = -9999, double yCUT_U = -9999, double zmin_U = -9999, double zmax_U = -9999 , double RMC_U = -9999 );

 /*! \brief Returns the average muon beam polarization

  \param mom the momentum of the beam muon track
  \param year the year of data taking (2002, 2003 or 2004).

  */
  static double GetBeamPol( float mom, int year );

  /*! \brief Returns the depolarization factor

  \param q2 squared invariant mass of the virtual photon
  \param xBj Bjorken x
  \param y energy fraction carried by the virtual photon
  \param R R value (e.g. from r1990.F)
  \param dR error on R
  \param D resulting depolarization factor
  \param dD calculated error on D (if do_err==true)
  \param do_err flag which controls whether the error calculation should be
  done (defaults to 'yes')

  The code in this function is taken from the A1 analysis code (thanks to
  A. Korzenev). This function does not invoke any external FORTRAN code, which
  means that if you have your own source for the R values you can use it even
  when compiling Phast with <code>NO_FORTRAN=1</code>.

  \author Roland.Kuhn@cern.ch
  */
  static void GetDepolarizationFactor( double q2, double xBj, double y,
                                       double R, double dR,
                                       double&D, double&dD,
                                       bool do_err = true);

#ifndef NO_FORTRAN

  /*! \brief Returns the depolarization factor
  \param q2 squared invariant mass of the virtual photon
  \param xBj Bjorken x
  \param y energy fraction carried by the virtual photon
  \param D resulting depolarization factor
  \param dD calculated error on D (if do_err==true)
  \param do_err flag which controls whether the error calculation should be
  done (defaults to 'yes')

  <b>This function is not available when Phast is compiled with
  <code>NO_FORTRAN=1</code></b>

  This function is a convenience wrapper for the first variant which takes the
  value of R+-dR from PaAlgo::GetR().

  \author Roland.Kuhn@cern.ch
  */
  static void GetDepolarizationFactor( double q2, double xBj, double y,
                                       double&D, double&dD,
                                       bool do_err = true);

  /*! \brief Returns the dilution factor and its error

  \param xBj Bjorken x
  \param y energy fraction carried by the virtual photon
  \param Cell target cell (U=upstream, C=central, D=downstream)
  \param run run number.
  \param flag processing flag (1=inclusive, 2=hadron_tagged)
  \param f dilution factor
  \param df error for the dilution factor

  <b>This function is not available when Phast is compiled with
  <code>NO_FORTRAN=1</code></b>

  This function calls the FORTRAN routine dfdeut from
  <code>./fortran/dilut_main.F</code> to get the dilution factor for the various
  regions of the target.

  \b Beware: The FORTRAN routine gives only <code>float</code> values!

  \author Alexandre.Korzenev@cern.ch
  */
  static void GetDilutionFactor( float xBj, float y, char Cell, int run, int flag,
                                 float &f, float &df);


  /* \brief Returns the dilution factor

  If you are not sure what you need, use better the function
  static void PaAlgo::GetDilutionFactor( float xBj, float y, char Cell, int run, int flag, float &f, float &df);

  \param x Bjorken x
  \param y energy fraction carried by the virtual photon
  \param flag processing flag (1=inclusive, 2=hadron_tagged)
  \param dftot pointer to an array of 11 \b floats to hold the result
  \param err_dftot pointer to an array of 11 \b floats to hold the errors

  <b>This function is not available when Phast is compiled with
  <code>NO_FORTRAN=1</code></b>

  This function calls the FORTRAN routine dfdeut from
  <code>./fortran/dfdeut.F</code> to get the dilution factor for the various
  regions of the target.

  For the 2-cells target the following information is returned:

  \arg dftot[0] overall dilution factor for entire target
  \arg dftot[1] overall dilution factor for upstream target cell
  \arg dftot[2] overall dilution factor for downstream target cell
  \arg dftot[3] 0<= \a r <= 8mm dilution factor for upstream target cell
  \arg dftot[4] 0<= \a r <= 8mm dilution factor for downstream target cell
  \arg dftot[5] 8mm < \a r <= 15mm dilution factor for upstream target cell
  \arg dftot[6] 8mm < \a r <= 15mm dilution factor for downstream target cell
  \arg dftot[7] 15mm < \a r <= 16mm dilution factor for upstream target cell
  \arg dftot[8] 15mm < \a r <= 16mm dilution factor for downstream target cell
  \arg dftot[9] 16mm < \a r <= 35mm dilution factor for upstream target cell
  \arg dftot[10] 16mm < \a r <= 35mm dilution factor for downstream target cell

  For the 3-cells target:

  \arg dftot[0] overall dilution factor for entire target
  \arg dftot[1] overall dilution factor for upstream target cell
  \arg dftot[2] overall dilution factor for central target cell
  \arg dftot[3] overall dilution factor for downstream target cell


  See comments in FORTRAN code for further details.

  \b Beware: The FORTRAN routine gives only <code>float</code> values!

  \author Roland.Kuhn@cern.ch
  */
  static void GetDilutionFactor( float x, float y, int flag, float*dftot,
                                 float*err_dftot);


  /*! \brief Returns the R1990 value

  \param q2 squared invariant mass of the virtual photon
  \param xBj Bjorken x
  \param R result (\b float!)
  \param dR error on R (if do_err==true)
  \param do_err flag which controls whether the error calculation should be
  done (defaults to 'yes')

  <b>This function is not available when Phast is compiled with
  <code>NO_FORTRAN=1</code></b>

  This function calls the FORTRAN routing r1990 from
  <code>./fortran/r1990.F</code> to get the R parameterization.

  \b Beware: The FORTRAN routine gives only <code>float</code> values!

  \author Roland Kuhn <roland.kuhn@cern.ch>
  */
  static void GetR( float q2, float xBj, float&R, float&dR, bool do_err = true);

  /*! \brief Returns the radiative correction weight (DEPRECATED)

  \deprecated This function will exit PHAST with an error message. Please use instead the function:
  static float PaAlgo::GetRadiativeWeight(float xBj, float y, int flag)

  \param xBj Bjorken x
  \param y energy fraction carried by the virtual photon
  \param flag processing flag (1=inclusive(deuter), 2=hadron_tagged(deuter),
  11=inclusive(hydrogen), 12=hadron_tagged(hydrogen))
  \param rc resulting radiative correction weight it is defined as:
  \f[
  rc = 1/\eta (x,y) = \frac{\sigma_{measured}}{\sigma_{1\gamma}}
  \f]

  <b>This function is not available when Phast is compiled with
  <code>NO_FORTRAN=1</code></b>

  This function calls the FORTRAN routine rcinlx or rcinhd from
  <code>./fortran/dfdeut.F</code> to get the radiative correction weights from
  precalculated tables. Tables that are used were generated for dilution factor
  calculations.

  \b Beware: The FORTRAN routine gives only <code>float</code> values!

  \author Konrad.Klimaszewski@cern.ch
  */
  static void GetRadiativeWeight(float xBj, float y, int flag, float& rc);

  /*! \brief Returns the radiative correction weight

  \param xBj Bjorken x
  \param y energy fraction carried by the virtual photon
  \param flag processing flag (1=inclusive(deuter), 2=hadron_tagged(deuter),
  11=inclusive(hydrogen), 12=hadron_tagged(hydrogen))
  \return The radiative correction factor defined as:
  \f[
  \eta (x,y) = \frac{\sigma_{1\gamma}}{\sigma_{measured}}
  \f]
  This factor can be _both_ smaller and larger than 1.
  For the reference in case of related questions, please look
  \verbatim
  =========================================================================
  B. Badelek, D. Bardin, K. Kurek and K. Scholz, Z. Phys. C 66 (1995) 591
  =========================================================================
  and in particular the Eq (4) and Fig 5.
  \endverbatim

  <b>This function is not available when Phast is compiled with
  <code>NO_FORTRAN=1</code></b>

  This function calls the FORTRAN routine rcinlx or rcinhd from
  <code>./fortran/dfdeut.F</code> to get the radiative correction weights from
  precalculated tables. Tables that are used were generated for dilution factor
  calculations.

  \b Beware: The FORTRAN routine gives only <code>float</code> values!

  \author Konrad.Klimaszewski@cern.ch
  */
  static float GetRadiativeWeight(float xBj, float y, int flag);

#endif // NO_FORTRAN

  /*!
    \brief return "true" is track "t" is identified in MW1

    The number of clusters associated to the track before and after the absorber is stored in
    nma1 and nma2 respectively if they are not null pointers.
    \author Andrea Ferrero <aferrero@to.infn.it>
  */
  static bool     DoMW1ID(const PaTrack& t, int* nma1=0, int* nma2=0);

  /*!
    \brief returns the index of the candidate scattered muon identified in MW1

    Returns the best scattered muon candidate as identified by MW1.
    The track is selected among those coming out of the "best primary vertex",
    as given by iBestPrimaryVertex() function. The identification is performed
    only in the case of calorimetric trigger events.
    \author Andrea Ferrero <aferrero@to.infn.it>
  */
  static int      GetMW1ScatMuon(const PaEvent& e);

  /*! \brief Q<sup>2</sup> via momenta and Cos
    \param pmu0 the absolute value of the beam momentum
    \param pmu the absolute value of the scattered lepton momentum
    \param Cos the cos of angle between beam and scattered lepton
    \param mmu the mass of the lepton (default is muon)
    \author Alexandre.Korzenev@cern.ch
  */

  static double Q2( double pmu0, double pmu, double Cos, double mmu = 0.105658357 );

  /*! \brief Q<sup>2</sup> via 3-vectors
    \param pmu0 the beam vector
    \param pmu the scattered lepton vector
    \param mmu the mass of the lepton (default is muon)
    \author Alexandre.Korzenev@cern.ch
 */
  static double Q2( const TVector3& pmu0, const TVector3& pmu, double mmu = 0.105658357 );

  /*! \brief Q<sup>2</sup> via Lorentz vectors
    \param pmu0 the beam Lorentz vector
    \param pmu the scattered lepton Lorentz vector
    \author Alexandre.Korzenev@cern.ch
 */
  static double Q2( const TLorentzVector& pmu0, const TLorentzVector& pmu );


  /*! \brief x<sub>Bjorken</sub> via momenta and Cos
    \param pmu0 the absolute value of the beam momentum
    \param pmu the absolute value of the scattered lepton momentum
    \param Cos the cos of angle between beam and scattered lepton
    \param m_P the mass of the target particle (default is proton)
    \param mmu the mass of the scattering particle (default is muon)
    \author Alexandre.Korzenev@cern.ch
  */
  static double xbj( double pmu0, double pmu, double Cos, double m_P = 0.93827231, double mmu = 0.105658357 );

  /*! \brief x<sub>Bjorken</sub> via 3-vectors
    \param pmu0 the beam vector
    \param pmu the scattered particle vector
    \param m_P the mass of the target particle (default is proton)
    \param mmu the mass of the scattering particle (default is muon)
    \author Alexandre.Korzenev@cern.ch
  */
    static double xbj( const TVector3& pmu0, const TVector3& pmu, double m_P = 0.93827231, double mmu = 0.105658357 );

  /*! \brief x<sub>Bjorken</sub> via Lorentz vectors
    \param pmu0 the beam Lorentz vector
    \param pmu the scattered particle Lorentz vector
    \param m_P the mass of the target particle (default is proton)
    \author Alexandre.Korzenev@cern.ch
  */
  static double xbj( const TLorentzVector& pmu0, const TLorentzVector& pmu, double m_P = 0.93827231 );


  /*! \brief W<sup>2</sup> via momenta and Cos
    \param pmu0 the absolute value of the beam momentum
    \param pmu the absolute value of the scattered lepton momentum
    \param Cos the cos of angle between beam and scattered lepton
    \param m_P the mass of the target particle (default is proton)
    \param mmu the mass of the scattering particle (default is muon)
    \author Alexandre.Korzenev@cern.ch
  */
  static double W2( double pmu0, double pmu, double Cos, double m_P = 0.93827231, double mmu = 0.105658357 );

  /*! \brief W<sup>2</sup> via  3-vectors
    \param pmu0 the beam vector
    \param pmu the scattered particle vector
    \param m_P the mass of the target particle (default is proton)
    \param mmu the mass of the scattering particle (default is muon)
    \author Alexandre.Korzenev@cern.ch
  */
  static double W2( const TVector3& pmu0, const TVector3& pmu, double m_P = 0.93827231, double mmu = 0.105658357 );

  /*! \brief W<sup>2</sup> via  Lorentz vectors
    \param pmu0 the beam Lorentz vector
    \param pmu the scattered particle Lorentz vector
    \param m_P the mass of the target particle (default is proton)
    \author Alexandre.Korzenev@cern.ch
  */
  static double W2( const TLorentzVector& pmu0, const TLorentzVector& pmu, double m_P = 0.93827231);


  /*! \brief boost via 3-vectors

    The Lorentz transformation of the particle {\a v1, \a m1} to the rest frame of
    the particle {\a v0, \a m0}.
    \param v0 the momentum vector of particle <b>0</b> in laboratory frame
    \param m0 the mass of the particle <b>0</b>
    \param v1 the momentum vector of particle <b>1</b> in laboratory frame
    \param m1 the mass of the particle <b>1</b>
    \param v10 the momentum vector of the particle <b>1</b> in the rest frame of <b>0</b>
    \author Alexandre.Korzenev@cern.ch
  */
  static void Boost( const TVector3& v0, double m0, const TVector3& v1, double m1, TVector3& v10 );

  /*! \brief boost via Lorentz vectors

    The Lorentz transformation of the particle \a v1 to the rest frame of
    the particle \a v0.
    \param v0 the four-vector of particle <b>0</b> in laboratory frame
    \param v1 the four-vector of particle <b>1</b> in laboratory frame
    \param v10 the momentum vector of the particle <b>1</b> in the rest frame of <b>0</b>
    \author Alexandre.Korzenev@cern.ch
 */
  static void Boost( const TLorentzVector& v0, const TLorentzVector& v1, TVector3& v10 );

  /*! \brief boost via Lorentz vectors

    The Lorentz transformation of the particle \a v1 to the rest frame of
    the particle \a v0.
    \param v0 the four-vector of particle <b>0</b> in laboratory frame
    \param v1 the four-vector of particle <b>1</b> in laboratory frame
    \param v10 the four-vector of the particle <b>1</b> in the rest frame of <b>0</b>
    \author Alexandre.Korzenev@cern.ch
 */
  static void Boost( const TLorentzVector& v0, const TLorentzVector& v1, TLorentzVector& v10 );


  /*! \brief x<sub>F</sub> via 3-vectors
    \param pmu0 the beam vector
    \param pmu the scattered lepton vector
    \param plh the hadron vector
    \param mh  the hadron mass
    \param m_P the target paticle mass (default is proton)
    \param mmu the lepton mass
    \author Alexandre.Korzenev@cern.ch
 */
  static double Xf( const TVector3& pmu0, const TVector3& pmu, const TVector3& plh,
		    double mh = 0.13957018, double m_P = 0.93827231, double mmu = 0.105658357 );

  /*! \brief x<sub>F</sub> via  Lorentz vectors
    \param pmu0 the beam four-vector
    \param pmu the scattered lepton four-vector
    \param plh the hadron four-vector
    \param m_P the target particle mass (default is proton)
    \author Alexandre.Korzenev@cern.ch
 */
  static double Xf( const TLorentzVector& pmu0, const TLorentzVector& pmu, const TLorentzVector& plh,
		    double m_P = 0.93827231 );

  //! Energy via mass and 3-vector
  static double E( double m, const TVector3& p ) { return sqrt(m*m+p.Mag2()); }

  //! Mass via energy and 3-vector
  static double M( double e, const TVector3& p ) { return sqrt(e*e-p.Mag2()); }

  //! Trajectory propagation in magnetic field by Runge-Kutta method with Jacobian calculations.
  static bool RKutta (double* SU, double* VO, double& Path);

  //! Change the energy part of the Lorentz vector "lv" to corresponds exactly to mass "mass".
  static void RedefMass(TLorentzVector& lv, double mass);

  /*
    \brief Inversion of a positively definite symmetric 5x5 matrix,
    represented by it's lower triangle

    \param a input matrix (array of 15)
    \param b inverted matrix (array of 15)

    \note for use in Kalman fit
  */
  static int Inv5(double*, double*);

  // vertex finder
  static bool FindVtx(const vector<PaTPar>& par, double zmin, double zmax,  double& zout, double& ezout, double& dist, int& niter);

  /*!
    \brief "M out of N" combinator

    If to call in the loop (while it returns "true"), on every call this function put in array "i[M]"
    indexes [0-N] for "M out of N" combinations.
    Could be used in combination of tracks for vertex candidates,
    in combination of tracks or calorimeter clusters for invariant mass calculation etc.
  */
  static bool CombMofN(int N, int M, int i[]);

  /*!
    \brief "N nested loops" combinator

    If to call in the loop (while it returns "true"), on every call this function put in array "ind[Nloops]
    running indexes of N nested loops where every loop do Niter cycles.
  */
  static bool CombNloops(int Nloops, int Niter, int ind[]);

  /*!
    \brief return interpolated pol.
    \param ev PaEvent object
    \param z  z coordinate where you want to calculate pol. at
    \details only pol. in 2015 is available
    \author Nukazuka, Genki [genki@quark.kj.yamagata-u.ac.jp]
  */
  static double GetDYtargetPolarization(const PaEvent& ev, double z );


  /*!
    \brief return a boolean and the dilution factor as well as its associated error by reference
    \param x2        Proton Bjorken x
    \param qT        Transverse momentum of the virtual photon
    \param Q2        Di-muon invariant mass SQUARED
    \param cell      Target cell: 'U', 'D' for upstream, downstream respectively
    \param dilution  Dilution factor
    \param edilution Error on the dilution factor
    \author vincent.andrieux@cern.ch
  */
  static bool GetDYdilutionFactor(const PaEvent& ev, double x2, double qT, double Q2, char cell, double& dilution, double& edilution);

}; // end of class definition

// inlined functions for kinematic variables calculations

inline double PaAlgo::Q2( double pmu0, double pmu, double Cos, double mmu) {
  double pp = pmu0*pmu*Cos;
  double E0 = sqrt( pmu0*pmu0 + mmu*mmu );
  double E  = sqrt( pmu *pmu  + mmu*mmu );
  return 2 * ( - mmu * mmu - pp + E0 * E );
}

inline double PaAlgo::Q2( const TVector3& pmu0, const TVector3& pmu, double mmu) {
  double pp = pmu0 * pmu;
  double E0 = PaAlgo::E( mmu, pmu0 );
  double E  = PaAlgo::E( mmu, pmu  );
  return 2 * ( - mmu * mmu - pp + E0 * E );
}

inline double PaAlgo::Q2( const TLorentzVector& pmu0, const TLorentzVector& pmu ) {
  return PaAlgo::Q2( pmu0.Vect(), pmu.Vect(), pmu0.M() );
}



inline double PaAlgo::xbj( double pmu0, double pmu, double Cos, double m_P, double mmu) {
  double q2 = PaAlgo::Q2( pmu0, pmu, Cos, mmu );
  double E0 = sqrt( pmu0*pmu0 + mmu*mmu );
  double E  = sqrt( pmu *pmu  + mmu*mmu );
  return q2 / ( 2 * m_P * ( E0 - E ) );
}

inline double PaAlgo::xbj( const TVector3& pmu0, const TVector3& pmu, double m_P, double mmu) {
  double q2 = PaAlgo::Q2( pmu0, pmu, mmu );
  double E0 = PaAlgo::E( mmu, pmu0 );
  double E  = PaAlgo::E( mmu, pmu  );
  return q2 / ( 2 * m_P * ( E0 - E ) );
}

inline double PaAlgo::xbj( const TLorentzVector& pmu0, const TLorentzVector& pmu, double m_P) {
  return PaAlgo::xbj( pmu0.Vect(), pmu.Vect(), m_P, pmu0.M() );
}


inline double PaAlgo::W2( double pmu0, double pmu, double Cos, double m_P, double mmu) {
  double q2 = PaAlgo::Q2( pmu0, pmu, Cos, mmu );
  double E0 = sqrt( pmu0*pmu0 + mmu*mmu );
  double E  = sqrt( pmu *pmu  + mmu*mmu );
  return m_P*m_P + 2*m_P*(E0-E) - q2;
}

inline double PaAlgo::W2( const TVector3& pmu0, const TVector3& pmu, double m_P, double mmu) {
  double q2 = PaAlgo::Q2( pmu0, pmu , mmu );
  double E0 = PaAlgo::E( mmu, pmu0 );
  double E  = PaAlgo::E( mmu, pmu  );
  return m_P*m_P + 2*m_P*(E0-E) - q2;
}

inline double PaAlgo::W2( const TLorentzVector& pmu0, const TLorentzVector& pmu, double m_P) {
  return PaAlgo::W2( pmu0.Vect(), pmu.Vect(), m_P, pmu0.M() );
}



inline void PaAlgo::Boost( const TVector3& v0, double m0, const TVector3& v1, double m1, TVector3& v10 ) {
  double E0 = PaAlgo::E( m0, v0 );
  double E1 = PaAlgo::E( m1, v1 );
  double E10=(E0*E1-v0*v1)/m0;
  v10=v1-(E1+E10)/(E0+m0)*v0;
}

inline void PaAlgo::Boost( const TLorentzVector& v0, const TLorentzVector& v1, TVector3& v10 ) {
  PaAlgo::Boost(v0.Vect(),v0.M(), v1.Vect(),v1.M(), v10 );
}

inline void PaAlgo::Boost( const TLorentzVector& v0, const TLorentzVector& v1, TLorentzVector& v10 ) {
  TVector3 v;
  PaAlgo::Boost(v0.Vect(),v0.M(), v1.Vect(),v1.M(), v );
  v10.SetVectM(v,v1.M());
}



inline double PaAlgo::Xf( const TVector3& pmu0, const TVector3& pmu, const TVector3& plh,
			  double mh, double m_P, double mmu) {

  double xf;
  double ElCM,mCM, Eg,Elg, Eh,Elh;
  TVector3 pCM, pg,plg, ph;

  plg = pmu0 - pmu;
  Elg = PaAlgo::E(mmu,pmu0) - PaAlgo::E(mmu,pmu);

  const TVector3& plCM = plg;
  ElCM = Elg + m_P;
  mCM = PaAlgo::M( ElCM, plCM );

  Eg = ( ElCM * Elg - plCM * plg ) / mCM;
  pg = plg - ( Elg + Eg )/( ElCM + mCM ) * plCM;

  Elh = PaAlgo::E( mh, plh );
  Eh = (ElCM * Elh - plh * plCM ) / mCM;
  ph = plh - ( Elh + Eh )/( ElCM + mCM ) * plCM;

  xf = 2 * pg.Unit() * ph / mCM;

  // if( fabs(xf) > 2 ) cout<<" |xf|>2, xf="<<xf<<endl;

  return xf;
}

inline double PaAlgo::Xf( const TLorentzVector& pmu0, const TLorentzVector& pmu,
			  const TLorentzVector& plh, double m_P) {
  return PaAlgo::Xf( pmu0.Vect(),pmu.Vect(),plh.Vect(), plh.M(),m_P,pmu0.M() );
}

inline void PaAlgo::RedefMass(TLorentzVector& lv, double mass){

  TLorentzVector lvtmp(lv.Vect(), sqrt((lv.Vect()).Mag2() + mass*mass));
  lv = lvtmp;
}

#endif
