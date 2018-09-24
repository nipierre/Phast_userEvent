#include <Phast.h>


class TargetCell
{
  std::vector<float> xv, yv, zv;
  int year;
  TargetCell();

 public:
  void Init();

  float PathLength(PaTrack* track, float zmin, float zmax, float R);

  bool InTarget(const PaVertex& vtx, float R);

  bool CrossCells(const PaTrack& track, float zmin=-325, float zmax=-71, float R=1.9, float Y=1.2);

  void CellCenter(float z, float& xc, float& yc);
};
