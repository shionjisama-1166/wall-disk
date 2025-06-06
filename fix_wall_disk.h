/* -*- c++ -*- ------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
------------------------------------------------------------------------- */

#ifndef LMP_FIX_WALL_DISK_H           // ❶ new include-guard
#define LMP_FIX_WALL_DISK_H

#include "fix_wall.h"

#if !defined(FIX_CLASS)  
namespace LAMMPS_NS {

/* ---------------------------------------------------------------------- */

class FixWallDisk : public FixWall {  // ❷ new class name
 public:
  FixWallDisk(class LAMMPS *, int, char **);
  ~FixWallDisk() override;

  void init() override;
  void precompute(int) override;
  void wall_particle(int, int, double) override;

 private:
  double coeff1[6], coeff2[6], coeff3[6], coeff4[6];
};

/* ---------------------------------------------------------------------- */

} // namespace LAMMPS_NS
#endif    

#endif /* LMP_FIX_WALL_DISK_H */

/* ---- register keyword so it appears in style_fix.h ------------------- */
#ifdef FIX_CLASS
// clang-format off
FixStyle(wall/disk,FixWallDisk);      // ❸ new keyword + class
// clang-format on
#endif