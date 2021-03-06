#include <Phast.h>


class TargetCell
{
  std::vector<float> xv, yv, zv;
  std::vector<float> xMCv, yMCv, zMCv;
  int year;

 public:

  TargetCell();

  void Init();

  float PathLength(PaTrack* track, float zmin, float zmax, float R);

  bool InTarget(const PaVertex& vtx, float R);
  bool InTargetMC(const PaMCvertex& vtx, float R);

  bool CrossCells(const PaTrack& track, float zmin=-325, float zmax=-71, float R=1.9, float Y=1.2);
  bool CrossCellsMC(const PaMCtrack& track, float zmin=-325, float zmax=-71, float R=1.9, float Y=1.2);

  void CellCenter(float z, float& xc, float& yc);
  void CellCenterMC(float z, float& xc, float& yc);
};
