#include <iostream>
#include <iomanip>

#include "TSpline.h"
#include "TRandom.h"

#include "Phast.h"
#include "PaAlgo.h"

double target_file_2012[25][9] =
{{-311.19, -0.360342, 0.00112758, 142.793, 0.00353858, 1.99989, 0.0012368, -0.04, -0.45},
{-301.19, -0.360342, 0.00112758, 142.793, 0.00353858, 1.99989, 0.0012368, -0.05, -0.4549},
{-291.19, -0.466043, 0.000601126, 142.799, 0.00168855, 1.99989, 0.000651646, -0.0665722, -0.461263},
{-281.19, -0.47561, 0.000614169, 142.798, 0.00156694, 1.99945, 0.000624622, -0.0685443, -0.470645},
{-271.19, -0.487868, 0.00064998, 142.821, 0.00155072, 2.00121, 0.000648421, -0.0589272, -0.484296},
{-261.19, -0.500966, 0.000666676, 142.849, 0.00153378, 2.00063, 0.000659744, -0.0469298, -0.498763},
{-251.19, -0.513325, 0.000788929, 142.882, 0.00160906, 1.99941, 0.000759828, -0.0312396, -0.512374},
{-241.19, -0.523261, 0.000799581, 142.916, 0.00155487, 2.00013, 0.000752633, -0.0136532, -0.523082},
{-231.19, -0.531869, 0.000864788, 142.952, 0.00155671, 1.99959, 0.00082212, 0.00480739, -0.531847},
{-221.19, -0.537759, 0.000901252, 142.981, 0.00153163, 1.99795, 0.000847752, 0.0208353, -0.537355},
{-211.19, -0.548344, 0.00084799, 143.01, 0.00140758, 2.00034, 0.000793529, 0.0372393, -0.547078},
{-201.19, -0.562294, 0.0010845, 143.031, 0.00149944, 2.00317, 0.00103106, 0.0499375, -0.560072},
{-191.19, -0.563077, 0.00113537, 143.037, 0.00150214, 1.99867, 0.00104078, 0.0533941, -0.560539},
{-181.19, -0.570993, 0.00102211, 143.052, 0.00148746, 2.0042, 0.000980157, 0.0622739, -0.567587},
{-171.19, -0.571367, 0.00111599, 143.054, 0.00155968, 2.00487, 0.00107006, 0.0634556, -0.567833},
{-161.19, -0.564802, 0.00104564, 143.047, 0.0013972, 1.99993, 0.000981339, 0.0588035, -0.561733},
{-151.19, -0.564016, 0.00109004, 143.039, 0.00158343, 1.99756, 0.00102487, 0.0543487, -0.561392},
{-141.19, -0.57844, 0.00107463, 143.042, 0.00146923, 2.00342, 0.00101544, 0.0574397, -0.575581},
{-131.19, -0.598211, 0.00114126, 143.045, 0.00168688, 2.00762, 0.00114181, 0.0613223, -0.59506},
{-121.19, -0.617062, 0.00117336, 143.038, 0.00157445, 2.00502, 0.001107, 0.0591593, -0.614219},
{-111.19, -0.639078, 0.0013255, 143.035, 0.00168094, 2.00466, 0.00127579, 0.059164, -0.636333},
{-101.19, -0.67542, 0.0014279, 143.016, 0.0017744, 2.00143, 0.00141485, 0.049698, -0.673589},
{-91.19, -0.714508, 0.00160412, 143.001, 0.00188937, 2.00333, 0.00163327, 0.0420778, -0.713268},
{-81.19, -0.760857, 0.00154153, 142.973, 0.00190479, 2.00511, 0.0015849, 0.0233096, -0.7605},
{-71.19, -0.795688, 0.00185853, 142.937, 0.00195403, 1.99213, 0.00170385, -0.0040598, -0.795677}};

double target_file_2016[29][9] =
{{-500.00, -0.137288, 0.000122812, -0.780854, 0.00120411, 1.97243, 9.98977e-05, 0.0975173, -0.0966352},
{-323.79, -0.137288, 0.000122812, -0.780854, 0.00120411, 1.97243, 9.98977e-05, 0.0975173, -0.0966352},
{-313.79, -0.143987, 0.000384912, -0.893519, 0.00325604, 1.97459, 0.000313727, 0.0902332, -0.112207},
{-303.79, -0.191253, 0.000544815, -1.16598, 0.00325718, 1.97785, 0.000448898, 0.0753255, -0.175795},
{-293.79, 0.24925, 0.00049329, 1.65587, 0.00204219, 1.99941, 0.000399223, 0.0211781, -0.248349},
{-283.79, 0.288479, 0.000476493, 1.55461, 0.0013603, 2.00425, 0.000333477, -0.00466937, -0.288442},
{-273.79, 0.327547, 0.00047984, 1.54251, 0.00113361, 2.00373, 0.000324566, -0.00926407, -0.327416},
{-263.79, 0.361748, 0.000457651, 1.55214, 0.000975106, 2.00548, 0.000297845, -0.00674825, -0.361685},
{-253.79, 0.391426, 0.000430639, 1.58079, 0.000943403, 2.00424, 0.000317716, 0.00391214, -0.391407},
{-243.79, 0.415539, 0.000419331, 1.61008, 0.000839509, 2.00367, 0.000305715, 0.0163214, -0.415219},
{-233.79, 0.436593, 0.000395061, 1.64343, 0.000816587, 2.00286, 0.000309984, 0.0316818, -0.435442},
{-223.79, 0.454741, 0.00041294, 309.547, 0.000760231, 2.0044, 0.000301802, 0.0453748, -0.452471},
{-213.79, 0.473328, 0.000425352, 1.70183, 0.000774725, 2.00344, 0.000323519, 0.0618448, -0.46927},
{-203.79, 0.490203, 0.00042692, 1.72634, 0.000743467, 2.00393, 0.000334269, 0.0759427, -0.484284},
{-193.79, 0.500946, 0.000388136, 1.74979, 0.000732282, 2.0017, 0.000318957, 0.0891882, -0.492943},
{-183.79, 0.505313, 0.000384736, 1.76893, 0.000732011, 2.00204, 0.00031878, 0.0994651, -0.495427},
{-173.79, -0.506375, 0.00040591, -1.36103, 0.000707547, 2.00261, 0.000336325, 0.105445, -0.495275},
{-163.79, -0.500812, 0.000406612, -1.35612, 0.000726149, 1.99944, 0.000333942, 0.10669, -0.489315},
{-153.79, -0.506437, 0.00042259, -51.6226, 0.000731611, 2.00289, 0.000343524, 0.107384, -0.494921},
{-143.79, 0.512317, 0.000421799, 1.77558, 0.00075051, 2.00167, 0.00034176, 0.104184, -0.501612},
{-133.79, -0.526122, 0.000462565, -1.37833, 0.000790483, 2.00257, 0.00036895, 0.100635, -0.516408},
{-123.79, -0.539493, 0.00044704, -1.39133, 0.000797124, 2.00144, 0.000366697, 0.0963038, -0.530828},
{-113.79, 0.551925, 0.000465019, 1.73143, 0.00074705, 2.00073, 0.000386386, 0.0882784, -0.544819},
{-103.79, 0.557702, 0.00049837, 1.71261, 0.000785211, 2.00007, 0.000425565, 0.0788264, -0.552103},
{-93.79, 0.557139, 0.000446873, 1.67956, 0.000806354, 2.00084, 0.000398517, 0.0604751, -0.553847},
{-83.79, 0.554022, 0.000512909, 1.63832, 0.000867301, 2.00302, 0.000456136, 0.0373791, -0.55276},
{-73.79, 0.548578, 0.000571343, 1.60366, 0.000906274, 2.00107, 0.000492774, 0.0180271, -0.548282},
{-70.00, 0.548578, 0.000571343, 1.60366, 0.000906274, 2.00107, 0.000492774, 0.0180271, -0.548282},
{-60.00, 0.548578, 0.000571343, 1.60366, 0.000906274, 2.00107, 0.000492774, 0.0180271, -0.548282}};

double target_file_2017[29][9] =
{{-500.00, 0.229376, 0.000370072, 0.301833, 0.00122449, 1.99753, 0.000287844, -0.219007, -0.0681867},
{-323.79, 0.229376, 0.000370072, 0.301833, 0.00122449, 1.99753, 0.000287844, -0.219007, -0.0681867},
{-313.79, 0.25106, 0.000729378, 0.224252, 0.00182379, 1.98916, 0.000583362, -0.244773, -0.05583},
{-303.79, 0.32493, 0.00133207, 0.319432, 0.00209533, 2.00523, 0.000775463, -0.308493, -0.102037},
{-293.79, 0.418489, 0.00131345, 0.417468, 0.00173224, 2.01256, 0.000848791, -0.382548, -0.169675},
{-283.79, 0.476534, 0.00124386, 0.46923, 0.00145491, 2.00976, 0.00082504, -0.425029, -0.215488},
{-273.79, 0.544642, 0.00152373, 0.516044, 0.00128627, 2.01046, 0.00108234, -0.473718, -0.26875},
{-263.79, 0.613247, 0.00160229, 0.553393, 0.00111616, 2.01432, 0.00122245, -0.521717, -0.322308},
{-253.79, 0.659311, 0.001469, 0.593512, 0.000990419, 2.00953, 0.00121947, -0.546557, -0.368737},
{-243.79, 0.694273, 0.0014252, 0.621761, 0.000851804, 2.00885, 0.0012413, -0.564343, -0.404392},
{-233.79, 0.719013, 0.00134434, 0.640164, 0.000755602, 2.01093, 0.00121023, -0.576647, -0.429486},
{-223.79, 0.72836, 0.00135553, 0.654241, 0.000712771, 2.00832, 0.00120105, -0.577961, -0.443249},
{-213.79, 0.737463, 0.00129198, 0.664886, 0.000682412, 2.00821, 0.00116151, -0.580374, -0.454992},
{-203.79, 0.744318, 0.00133629, 0.670825, 0.000687942, 2.00829, 0.00119207, -0.583031, -0.462692},
{-193.79, 0.745385, 0.00130927, 0.674474, 0.000691416, 2.00446, 0.0011625, -0.582172, -0.465483},
{-183.79, 0.760192, 0.00155123, 0.674049, 0.000776324, 2.00997, 0.0013847, -0.593938, -0.474477},
{-173.79, 0.775427, 0.00167067, 0.670818, 0.000803443, 2.01143, 0.00149919, -0.607402, -0.482026},
{-163.79, 0.785504, 0.00160032, 0.681309, 0.000772989, 2.00953, 0.00144564, -0.610139, -0.494719},
{-153.79, 0.799322, 0.00161653, 0.696316, 0.000769732, 2.01386, 0.00147667, -0.613248, -0.512682},
{-143.79, 0.819085, 0.00170752, 0.708409, 0.000749264, 2.02108, 0.00152867, -0.622012, -0.532919},
{-133.79, 0.820403, 0.00158103, 0.722297, 0.000734502, 2.00679, 0.00139682, -0.61554, -0.542376},
{-123.79, 0.826571, 0.00141132, 0.738858, 0.000686824, 1.99616, 0.00123041, -0.611033, -0.55665},
{-113.79, 0.856104, 0.00178234, 0.74203, 0.000744214, 2.00622, 0.00158856, -0.631033, -0.578542},
{-103.79, 0.875317, 0.00181212, 0.73541, 0.000727104, 2.00627, 0.0015757, -0.649096, -0.587242},
{-93.79, 0.896213, 0.00180573, 0.720407, 0.000736007, 2.00496, 0.00160658, -0.673538, -0.591223},
{-83.79, 0.925352, 0.00196474, 0.690233, 0.000765448, 2.01039, 0.0017562, -0.713536, -0.589187},
{-73.79, 0.956502, 0.00220583, 0.647156, 0.000843521, 2.02035, 0.00200773, -0.763099, -0.576694},
{-63.79, 0.770585, 0.00288422, 0.780979, 0.0015848, 1.86346, 0.00235879, -0.547288, -0.542473},
{-60.00, 0.770585, 0.00288422, 0.780979, 0.0015848, 1.86346, 0.00235879, -0.547288, -0.542473}};

double target_mc_2012[5] = {0.0, 0.0, -0.5, -0.00109, -198.5};

double target_mc_2016[5] = {0.0, 0.0, -0.5, -0.00161, -202.0};

double target_mc_2017[5] = {0.0, 0.0, -0.5, -0.00117, -202.5};

bool PaAlgo::GetTargetLocationCenter(int run, double &xC, double &yC, double &xCmc, double &yCmc, double z, double &R, double &RMC, double &yCUT)
{
  xC = 1000000;
  yC = 1000000;

  vector<double> xv, yv, zv;
  double xMC, phiMC, yMC, thetaMC, zMC;

  if( !(xv.size() && yv.size() && zv.size()) )  // Check if already initialized
  {
    std::ifstream fin, finmc;
    std::string tstr, tstrmc;

    if( 96224 <= run && run <= 109125 )
		{
			for(int i=0; i<25; i++)
      {
        zv.push_back(target_file_2012[i][0]);
        xv.push_back(target_file_2012[i][7]);
        yv.push_back(target_file_2012[i][8]);
      }
      xMC = target_mc_2012[0]; phiMC = target_mc_2012[1]; yMC = target_mc_2012[2];
      thetaMC = target_mc_2012[3]; zMC = target_mc_2012[4];
		}
    else if( 264860 <= run && run <= 276879 )
		{
      for(int i=0; i<29; i++)
      {
        zv.push_back(target_file_2016[i][0]);
        xv.push_back(target_file_2016[i][7]);
        yv.push_back(target_file_2016[i][8]);
      }
      xMC = target_mc_2016[0]; phiMC = target_mc_2016[1]; yMC = target_mc_2016[2];
      thetaMC = target_mc_2016[3]; zMC = target_mc_2016[4];
		}
    else if( 276880 <= run && run <= 281775 )
		{
      for(int i=0; i<29; i++)
      {
        zv.push_back(target_file_2017[i][0]);
        xv.push_back(target_file_2017[i][7]);
        yv.push_back(target_file_2017[i][8]);
      }
      xMC = target_mc_2017[0]; phiMC = target_mc_2017[1]; yMC = target_mc_2017[2];
      thetaMC = target_mc_2017[3]; zMC = target_mc_2017[4];
		}
    else return false; //check, otherwise segmentation fault
    cout << "PaAlgo::GetTargetLocationCenter(): Loaded RD/MC target file " << tstr << endl;
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
