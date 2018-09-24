#include <Phast.h>


class TargetCell
{
  std::vector<float> xv, yv, zv;
  int year;

  static TargetCell* ptr;

  TargetCell();
 public:
  static TargetCell& Instance();

  void Init(const PaEvent& e);

  float PathLength(PaTrack* track, float zmin, float zmax, float R);

  bool InTarget(const PaVertex& vtx, float R);

  bool InTarget(const PaTrack& track, float zmin, float zmax, float R, float Y);

  void CellCenter(float z, float& xc, float& yc);
};
