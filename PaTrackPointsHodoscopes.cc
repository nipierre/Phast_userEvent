#include "Phast.h"
#include "PaEvent.h"
#include "PaTrack.h"


#define MU_PRIM_HODO 0x1                 // Was PointsHodoscope called
#define MU_PRIM_HODO_RESULT 0x2          // PointsHodoscope result
#define MU_PRIM_ABS 0x4                  // Was CanBeMuon called
#define MU_PRIM_ABS_RESULT 0x8           // CanBeMuon result
#define MU_PRIM_SM2 0x10                 // Was crossYokeSM2(false) called
#define MU_PRIM_SM2_RESULT 0x20          // crossYokeSM2(false) result
#define MU_PRIM_SM2_STRICT 0x40          // Was crossYokeSM2(true) called
#define MU_PRIM_SM2_STRICT_RESULT 0x80   // crossYokeSM2(true) result
// Konrad.Klimaszewski@cern.ch


//LargeQ2 and LAST hodoscopes taken into account for 2007 and after 2007 respectively (Konrad)
//Assert expression added for 2012 Primakov data taking to treat muon and hadron data in the same way (Vincent)
bool PaTrack::PointsHodoscopes() const
// Konrad.Klimaszewski@cern.ch
{
	// If the result of the check was stored use it.
	if((muPrimInfo&MU_PRIM_HODO) != 0) {
		return ((muPrimInfo&MU_PRIM_HODO_RESULT)!=0);
	} else {
		muPrimInfo |= MU_PRIM_HODO;
		muPrimInfo |= MU_PRIM_HODO_RESULT;
	}

	// Extract information about the Hodoscope planes once and store it in std::map for fast access
	// If run number changes refresh the information
	static bool first = 1;
	static int run_number = 0;
	static int known_trigger_bits = 0; // Bits of known triggers
	static int outer_trigger_bits = 0; // Bits for OT and cOT if exists
	static int calo_trigger_bits = 0;  // Bits for CT, RPD and other hodoscope less triggers
	typedef map<string, PaDetect> HodosMap;
	typedef HodosMap::iterator HodosIter;
	static HodosMap Hodos;
	const PaEvent& e = Phast::Ref().event;
	int TriggerMask = e.TrigMask();
	if(first || run_number != e.RunNum()) {
		run_number = e.RunNum();

		// Most common bit masks
		known_trigger_bits = 0x31f;
		outer_trigger_bits = 0x8;
		calo_trigger_bits = 0x10;

		// Per year variations
		if(e.Year() == 2003 || e.Year() == 2004) {
			// In 2003 and 2004 there was a caloOT x200
			outer_trigger_bits = 0x208;
		} else if (e.Year() >= 2007 || e.IsMC()) {
			// x200 corresponds to LargeQ2 in 2007, LAST in >2010, RPD in 2008-09
			// Hadron Runs (treat all triggers as hodoscope less ones)
			if (e.Year() == 2008 || e.Year() == 2009) {
				known_trigger_bits = 0x21f;
				outer_trigger_bits = 0;
				calo_trigger_bits = 0x21f;
			}
			// DVCS 2009
			if (e.Year() == 2009 && e.RunNum() >= 79584 && e.RunNum() <= 79964) {
				known_trigger_bits = 0x21f;
				outer_trigger_bits = 0x8;
				calo_trigger_bits = 0x201;
			}
			// Primakoff and DVCS 2012
			if (e.Year() == 2012) {
			  assert(e.RunNum()< 102893 || e.RunNum() > 105692);// Primakov data taking : No checks done to treat muon and hadron data in the same way
			  known_trigger_bits = 0x20e;
			  outer_trigger_bits = 0x8;
			  calo_trigger_bits  = 0x10;
			}
			// DVCS 2016/2017
			if (e.Year() == 2016 || e.Year() == 2017) {
				known_trigger_bits = 0x20e;
				outer_trigger_bits = 0x8;
				calo_trigger_bits  = 0x10;
			}
		} else {
			known_trigger_bits = 0x11f;
		}

		static string HodoNames[21] =
		{
			"HO03Y1_m", // 0
			"HI04X1_u", // 1
			"HI04X1_d", // 2
			"HO04Y1_m", // 3
			"HO04Y2_m", // 4
			"HM04Y1_d", // 5
			"HM04Y1_u", // 6
			"HM04X1_d", // 7
			"HM04X1_u", // 8
			"HL04X1_m", // 9
			"HM05X1_d", // 10
			"HM05X1_u", // 11
			"HM05Y1_d", // 12
			"HM05Y1_u", // 13
			"HL05X1_m", // 14
			"HI05X1_d", // 15
			"HI05X1_u", // 16
			"HG01Y1__", // 17
			"HG02Y1__", // 18
			"HG02Y2__", // 19
			"HQ01Y1_m", // 20
		};
		int idet;

		for(int i=0; i<21; i++) {
		  if (HodoNames[i].find("HI") != string::npos && e.Year()>=2008 && e.Year()<=2009 && e.Year() == 2012) continue;
			else if (HodoNames[i].find("HQ") != string::npos && e.Year() != 2007) continue;
			else if (HodoNames[i].find("HG") != string::npos && e.RunNum() < 86202 && !e.IsMC()) continue;

			idet = PaSetup::Ref().iDetector(HodoNames[i]);
			if(idet>0) {
				Hodos[HodoNames[i]] = PaSetup::Ref().Detector(idet);
			} else {
				cerr<<"PaTrack::PointsHodoscopes Warning: Detector "<<HodoNames[i]<<" not found!"<<endl;
			}
		}
		first = 0;
	}
	TriggerMask &= known_trigger_bits;

	// No known trigger fired
	if(TriggerMask == 0)
		return 0;

	// Pure CT case no checking needed
	if((TriggerMask&(~calo_trigger_bits)) == 0)
		return 1;

	int Npars = NTPar();
	if(Npars == 0) {
		cerr<<"PaTrack::PointsHodoscopes ERROR: track does not have associated helixes!"<<endl;

		muPrimInfo &= ~MU_PRIM_HODO_RESULT;
		return 0;
	}

	PaTPar partr_last = vTPar(Npars-1); //Track parameters in last measured point
	PaTPar partr_first = vTPar(0);      //Track parameters just after SM2
	PaTPar partr_zero = vTPar(0);       //Track parameters in first measured point
	PaTPar partr_between_magnets = vTPar(0); //Track parameters between SM1 and SM2
	PaTPar partr;                       //Track parameters after extrapolation
	double zSM1 = PaSetup::Ref().PtrMagField()->getMagInfo()[1].zcm / 10.0; // SM1 position
	double zSM2 = PaSetup::Ref().PtrMagField()->getMagInfo()[2].zcm / 10.0; // SM2 position
	double zMF2 = 3800; // MF2 position
	bool MaterialMaps = false; // Should we use material maps for the extrapolations

	// In case of OT and IT events make a preextrapolation from the first measured point
	// to just after SM2 to speed things up. If the last helix is before the MF2 use it
	// instead of the first one.
	// For all triggers if the last helix is before SM2 preextrapolate it to just after SM2.
	if(partr_last(0) < zMF2) {
		if(partr_last(0) > zSM2) {
			partr_first = partr_last;
		} else {
			partr_last.Extrapolate(zSM2 + 250, partr, MaterialMaps);
			partr_first = partr;
			partr_last = partr;
		}
	} else if((TriggerMask&(~(outer_trigger_bits+0x1)))==0) {
		partr_first.Extrapolate(zSM2 + 250, partr, MaterialMaps);
		partr_first = partr;
	}
	if((e.IsMC() || e.Year()>=2010) && (TriggerMask&0x200)!=0) { // LAST is positioned between magnets
		partr_zero.Extrapolate(zSM1 + 220, partr, MaterialMaps);
		partr_between_magnets = partr;
	}

	string name;  // Hodo plane name
	HodosIter it; // Position in the Hodo detectors map.
	// In the trigger logic usually two planes have to give a hit for a trigger to fire.
	// h1 - did the first plane from a "trigger pair" fire
	// h1 - did the second plane from a "trigger pair" fire
	// fatal - set to true if a detector is missing
	bool h1, h2, fatal;

	// IT
	if((TriggerMask&0x1)!=0)
	{
		h1 = h2 = fatal = 0;
		name = "HI04X1_u";
		it = Hodos.find(name);
		if(it == Hodos.end()) {fatal = 1;}
		if(!fatal && partr_first.Extrapolate(it->second.Z(), partr, MaterialMaps)) {
			h1 = true;

			// HI04 is not perpendicular to the beam direction.
			// The InActive method gives wrong results for this plane.
			// Do the check by hand.
			double y = partr(2) - it->second.Y();
			if(fabs(y) > (it->second.YSiz()*it->second.Ca())/2.0) h1 = false;

			double uwrs = partr(1) - it->second.Uorig() + it->second.Pitch()/2;
			double srange = it->second.Pitch()*it->second.Nwires();
			if ( uwrs < min(0.,srange) || uwrs > max(0.,srange) ) h1 = false;
		}
		name = "HI05X1_u";
		it = Hodos.find(name);
		if(it == Hodos.end()) {fatal = 1;}
		if(!fatal && h1 && partr_last.Extrapolate(it->second.Z(), partr, MaterialMaps)) {
			h2 = it->second.InActive(partr(1), partr(2));
		}
		if(h1 && h2) return 1;

		h1 = h2 = 0;
		name = "HI04X1_d";
		it = Hodos.find(name);
		if(it == Hodos.end()) {fatal = 1;}
		if(!fatal && partr_first.Extrapolate(it->second.Z(), partr, MaterialMaps)) {
			h1 = true;

			// HI04 is not perpendicular to the beam direction.
			// The InActive method gives wrong results for this plane.
			// Do the check by hand.
			double y = partr(2) - it->second.Y();
			if(fabs(y) > (it->second.YSiz()*it->second.Ca())/2.0) h1 = false;

			double uwrs = partr(1) - it->second.Uorig() + it->second.Pitch()/2;
			double srange = it->second.Pitch()*it->second.Nwires();
			if ( uwrs < min(0.,srange) || uwrs > max(0.,srange) ) h1 = false;
		}
		name = "HI05X1_d";
		it = Hodos.find(name);
		if(it == Hodos.end()) {fatal = 1;}
		if(!fatal && h1 && partr_last.Extrapolate(it->second.Z(), partr, MaterialMaps)) {
			h2 = it->second.InActive(partr(1), partr(2));
		}
		if(h1 && h2) return 1;
	}

	// MT
	if((TriggerMask&0x102)!=0)
	{
		fatal = h1 = h2 = 0;
		if(e.Year() == 2016 || e.Year() == 2017)
		{
			cout << "pouet" << endl;
			name = "HM04Y1_d";
			it = Hodos.find(name);
			if(it == Hodos.end()) {fatal = 1;}
			if(!fatal && partr_last.Extrapolate(it->second.Z(), partr, MaterialMaps)) {
				h1 = it->second.InActive(partr(1), partr(2));
			}

			if(h1)
			{
				name = "HM05Y1_d";
				it = Hodos.find(name);
				if(it == Hodos.end()) {fatal = 1;}
				if(!fatal && h1 && partr_last.Extrapolate(it->second.Z(), partr, MaterialMaps)) {
					h1 = it->second.InActive(partr(1), partr(2));
				}
			}
			else { h1=0; }

			name = "HM04Y1_u";
			it = Hodos.find(name);
			if(it == Hodos.end()) {fatal = 1;}
			if(!fatal && !h1 && partr_last.Extrapolate(it->second.Z(), partr, MaterialMaps)) {
				h2 = it->second.InActive(partr(1), partr(2));
			}

			if(h2)
			{
				name = "HM05Y1_u";
				it = Hodos.find(name);
				if(it == Hodos.end()) {fatal = 1;}
				if(!fatal && h2 && partr_last.Extrapolate(it->second.Z(), partr, MaterialMaps)) {
					h2 = it->second.InActive(partr(1), partr(2));
				}
			}
			else { h2=0; }
			if(h1 || h2) return 1;
		}
		else
		{
			name = "HM04X1_u";
			it = Hodos.find(name);
			if(it == Hodos.end()) {fatal = 1;}
			if(!fatal && partr_last.Extrapolate(it->second.Z(), partr, MaterialMaps)) {
				h1 = it->second.InActive(partr(1), partr(2));
			}

			if(h1)
			{
				name = "HM05X1_u";
				it = Hodos.find(name);
				if(it == Hodos.end()) {fatal = 1;}
				if(!fatal && h1 && partr_last.Extrapolate(it->second.Z(), partr, MaterialMaps)) {
					h1 = it->second.InActive(partr(1), partr(2));
				}
			}
			else { h1=0; }

			name = "HM04X1_d";
			it = Hodos.find(name);
			if(it == Hodos.end()) {fatal = 1;}
			if(!fatal && !h1 && partr_last.Extrapolate(it->second.Z(), partr, MaterialMaps)) {
				h2 = it->second.InActive(partr(1), partr(2));
			}

			if(h2)
			{
				name = "HM05X1_d";
				it = Hodos.find(name);
				if(it == Hodos.end()) {fatal = 1;}
				if(!fatal && h2 && partr_last.Extrapolate(it->second.Z(), partr, MaterialMaps)) {
					h2 = it->second.InActive(partr(1), partr(2));
				}
			}
			else { h2=0; }

			if(h1 || h2) {
				h1 = h2 = 0;
				name = "HM04Y1_d";
				it = Hodos.find(name);
				if(it == Hodos.end()) {fatal = 1;}
				if(!fatal && partr_last.Extrapolate(it->second.Z(), partr, MaterialMaps)) {
					h1 = it->second.InActive(partr(1), partr(2));
				}

				if(h1)
				{
					name = "HM05Y1_d";
					it = Hodos.find(name);
					if(it == Hodos.end()) {fatal = 1;}
					if(!fatal && h1 && partr_last.Extrapolate(it->second.Z(), partr, MaterialMaps)) {
						h1 = it->second.InActive(partr(1), partr(2));
					}
				}
				else { h1=0; }

				name = "HM04Y1_u";
				it = Hodos.find(name);
				if(it == Hodos.end()) {fatal = 1;}
				if(!fatal && !h1 && partr_last.Extrapolate(it->second.Z(), partr, MaterialMaps)) {
					h2 = it->second.InActive(partr(1), partr(2));
				}

				if(h2)
				{
					name = "HM05Y1_u";
					it = Hodos.find(name);
					if(it == Hodos.end()) {fatal = 1;}
					if(!fatal && h2 && partr_last.Extrapolate(it->second.Z(), partr, MaterialMaps)) {
						h2 = it->second.InActive(partr(1), partr(2));
					}
				}
				else { h2=0; }
				if(h1 || h2) return 1;
			}
		}
	}

	// LT
	if((TriggerMask&0x4)!=0)
	{
		fatal = h1 = h2 = 0;
		name = "HL04X1_m";
		it = Hodos.find(name);
		if(it == Hodos.end()) {fatal = 1;}
		if(!fatal && partr_last.Extrapolate(it->second.Z(), partr, MaterialMaps)) {
			h1 = it->second.InActive(partr(1), partr(2));
		}
		name = "HL05X1_m";
		it = Hodos.find(name);
		if(it == Hodos.end()) {fatal = 1;}
		if(!fatal && h1 && partr_last.Extrapolate(it->second.Z(), partr, MaterialMaps)) {
			h2 = it->second.InActive(partr(1), partr(2));
		}
		if(h1 && h2) return 1;
	}

	// OT
	if((TriggerMask&outer_trigger_bits)!=0)
	{
		fatal = h1 = h2 = 0;
		name = "HO03Y1_m";
		it = Hodos.find(name);
		if(it == Hodos.end()) {fatal = 1;}
		if(!fatal && partr_first.Extrapolate(it->second.Z(), partr, MaterialMaps)) {
			h1 = it->second.InActive(partr(1), partr(2));
			// For HO03 InActive properly describes larger rectangles.
			// Take into account two horizontal missing strips.
			// Pictures in documentation are correct only for 2002.
			if(partr(2) < it->second.Y() + 6.5 && partr(2) > it->second.Y() - 6.5)
				h1 = false;
		}
		// For HO04 InActive properly describes larger rectangles.
		// Take into account two horizontal missing strips. (See pictures ...)
		name = "HO04Y1_m";
		it = Hodos.find(name);
		if(it == Hodos.end()) {fatal = 1;}
		if(!fatal && h1 && partr_last.Extrapolate(it->second.Z(), partr, MaterialMaps)) {
			h2 = it->second.InActive(partr(1), partr(2));
			if(partr(2) < it->second.Y() + 14 && partr(2) > it->second.Y() - 14 &&
					partr(1) < it->second.X() + 190)
				h2 = false;
		}
		name = "HO04Y2_m";
		it = Hodos.find(name);
		if(it == Hodos.end()) {fatal = 1;}
		if(!fatal && h1 && !h2 && partr_last.Extrapolate(it->second.Z(), partr, MaterialMaps)) {
			h2 = it->second.InActive(partr(1), partr(2));
			if(partr(2) < it->second.Y() + 14 && partr(2) > it->second.Y() - 14)
				h2 = false;
		}
		if(h1 && h2) return 1;
	}

	// HighQ2
	if((TriggerMask&0x200)!=0 and e.Year() == 2007)
	{
		fatal = h1 = 0;
		name = "HQ01Y1_m";
		it = Hodos.find(name);
		if(it == Hodos.end()) {fatal = 1;}
		if(!fatal && partr_zero.Extrapolate(it->second.Z(), partr, MaterialMaps)) {
			h1 = it->second.InActive(partr(1), partr(2));
		}
		if(h1) return 1;
	}

	// LAST
	if((TriggerMask&0x200)!=0 and (e.Year() >= 2010 || e.IsMC()))
	{
		fatal = h1 = h2 = 0;
		name = "HG01Y1__";
		it = Hodos.find(name);
		if(it == Hodos.end()) {fatal = 1;}
		if(!fatal && partr_between_magnets.Extrapolate(it->second.Z(), partr, MaterialMaps)) {
			h1 = it->second.InActive(partr(1), partr(2));
		}
		name = "HG02Y1__";
		it = Hodos.find(name);
		if(it == Hodos.end()) {fatal = 1;}
		if(!fatal && h1 && partr_between_magnets.Extrapolate(it->second.Z(), partr, MaterialMaps)) {
			h2 = it->second.InActive(partr(1), partr(2));
		}
		if(!h2)
		{
			name = "HG02Y2__";
			it = Hodos.find(name);
			if(it == Hodos.end()) {fatal = 1;}
			if(!fatal && h1 && partr_between_magnets.Extrapolate(it->second.Z(), partr, MaterialMaps)) {
				h2 = it->second.InActive(partr(1), partr(2));
			}
		}
		if(h1 && h2) return 1;
	}

	// Store the result
	muPrimInfo &= ~MU_PRIM_HODO_RESULT;
	return 0;
}


bool PaTrack::CanBeMuon(void) const
// Konrad.Klimaszewski@cern.ch
{
	// Return stored result if it extists
	if((muPrimInfo&MU_PRIM_ABS) != 0) {
		return ((muPrimInfo&MU_PRIM_ABS_RESULT)!=0);
	} else {
		muPrimInfo |= MU_PRIM_ABS;
		muPrimInfo &= ~MU_PRIM_ABS_RESULT;
	}

	// Extract information about the Hodoscope planes once and store it in std::map for fast access
	// If run number changes refresh the information
	static bool first = true;
	static PaDetect det_d;
	static PaDetect det_u;
	static bool noHI = false;
	static int muon_sign;
	static int run_number = 0;
	const PaEvent& e = Phast::Ref().event;
	if(first || run_number != e.RunNum()) {
		run_number = e.RunNum();
		const PaSetup& setup = PaSetup::Ref();
		int idet_d = setup.iDetector("HI05X1_d");
		int idet_u = setup.iDetector("HI05X1_u");
		if( idet_d == -1 || idet_u == -1 ) noHI = true; // there no such detectors

		det_d = setup.Detector(idet_d);
		det_u = setup.Detector(idet_u);

		// Beam particle sign is determined based on SM2 field scale
		if(PaSetup::Ref().PtrMagField()->getMagInfo()[2].fsc < 0)
			muon_sign = -1;
		else
			muon_sign = 1;

		first = false;
	}

	double ZmuF2 = 3700; // entrance to muF2
	double ZHI05 = det_d.Z(); // 51 m
	int Npar = NTPar();
	PaTPar Hout;

	if( Q() != muon_sign ) return false; // for positive particles only

	if(Npar==0) {
		cerr<<"PaTrack::CanBeMuon ERROR: there are no helices in the track!"<<endl;
		return false;  // no helices in the track
	}

	const PaTPar& par = vTPar(Npar-1);

	if( ZLast() > 5000 ) goto can_be_muon_true; // hit in HI05 => muon

	if( par.Mom() < 5 ) return false; // too soft to be a subject of h-mu misidentification

	if( !par.Extrapolate ( ZmuF2, Hout, false) ) return false; // out of the hole region

	if( Hout(1) > 13 && Hout(1) < 52 && fabs(Hout(2)) < 18.5 ) { // entrance hole in muF2, numbers from COMGEANT
		// Check whether the particle in active zone of HI05
		if( noHI ) goto can_be_muon_true; // there no such detectors

		if( !par.Extrapolate ( ZHI05, Hout, false) ) goto can_be_muon_true; // out of HI05 region

		// stopped before HI05 => hadron
		if(!det_d.InActive(Hout(1),Hout(2)) && !det_u.InActive(Hout(1),Hout(2)))
			goto can_be_muon_true;
	}

	return false; // out of the hole region or stopped before HI05

	// Store positive result
can_be_muon_true:
	muPrimInfo |= MU_PRIM_ABS_RESULT;
	return true;
}



bool PaTrack::CrossYokeSM2(bool strict) const
{
	// Return stored result if it extists
	if(strict) {
		if((muPrimInfo&MU_PRIM_SM2_STRICT) != 0) {
			return ((muPrimInfo&MU_PRIM_SM2_STRICT_RESULT)!=0);
		} else {
			muPrimInfo |= MU_PRIM_SM2_STRICT;
			muPrimInfo &= ~MU_PRIM_SM2_STRICT_RESULT;
		}
	} else {
		if((muPrimInfo&MU_PRIM_SM2) != 0) {
			return ((muPrimInfo&MU_PRIM_SM2_RESULT)!=0);
		} else {
			muPrimInfo |= MU_PRIM_SM2;
			muPrimInfo &= ~MU_PRIM_SM2_RESULT;
		}
	}

	const PaSetup& setup = PaSetup::Ref();
	PaField *field = setup.PtrMagField();
	PaMagInfo *magInfo = field->getMagInfo();
	double xSM2 = magInfo[2].xcm / 10.; // in cm
	double ySM2 = magInfo[2].ycm / 10.;
	double zSM2 = magInfo[2].zcm / 10.;

	if( ZLast() < zSM2 ) return false;

	int NPAR = NTPar();
	if(NPAR==0) {
		cerr<<"PaTrack::CrossYokeSM2 ERROR: there are no helices in the track!"<<endl;
		return false;
	}
	// Select the last helix of a track. If it was indeed stored then the
	// extrapolation at least for downstream edge of the SM2 will not go
	// through magnetic field and as such will not be affected by wrong
	// momentum asignment.
	const PaTPar &par = vTPar(NPAR-1); // the last helix of the track

	double zEntr = zSM2 - 200;
	double zExit = zSM2 + 200;

	// In case where our last helix is indeed after the SM2 we will extrapolate
	// to the downstream face of SM2 first.
	if(par(0) > zSM2) {
		zEntr = zSM2 + 200;
		zExit = zSM2 - 200;
	}

	double xJura = xSM2 - 100;
	double xSale = xSM2 + 100;

	double yUp = ySM2 + 50;
	double yDw = ySM2 - 50;

	PaTPar parEntr;
	PaTPar parExit;
	if( !par.Extrapolate( zEntr, parEntr, strict) ) return false;

	double xTr = parEntr(1);
	double yTr = parEntr(2);
	double ExTr;
	double EyTr;

	if(strict) {
		ExTr = 2*sqrt( parEntr(1,1) );
		EyTr = 2*sqrt( parEntr(2,2) );
		if( xTr < xJura+ExTr || xTr > xSale-ExTr ) goto sm2_yoke_true_strict;
		if( yTr < yDw  +EyTr || yTr > yUp  -EyTr ) goto sm2_yoke_true_strict;
	} else {
		if( xTr < xJura || xTr > xSale ) goto sm2_yoke_true;
		if( yTr < yDw   || yTr > yUp   ) goto sm2_yoke_true;
	}

	if( !parEntr.Extrapolate( zExit, parExit, strict) ) {
		if(strict) goto sm2_yoke_true_strict;
		else goto sm2_yoke_true;
	}

	xTr = parExit(1);
	yTr = parExit(2);

	if(strict) {
		ExTr = 2*sqrt( parExit(1,1) );
		EyTr = 2*sqrt( parExit(2,2) );
		if( xTr < xJura+ExTr || xTr > xSale-ExTr ) goto sm2_yoke_true_strict;
		if( yTr < yDw  +EyTr || yTr > yUp  -EyTr ) goto sm2_yoke_true_strict;
	} else {
		if( xTr < xJura || xTr > xSale ) goto sm2_yoke_true;
		if( yTr < yDw   || yTr > yUp   ) goto sm2_yoke_true;
	}

	return false;

sm2_yoke_true_strict:
	muPrimInfo |= MU_PRIM_SM2_STRICT_RESULT;
	return true;

sm2_yoke_true:
	muPrimInfo |= MU_PRIM_SM2_RESULT;
	return true;
}
