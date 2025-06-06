/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Binghan Liu (Virginia Tech)
------------------------------------------------------------------------- */

#include "fix_wall_disk.h"

#include "atom.h"
#include "error.h"
#include "math_special.h"

#include <cmath>
#include "update.h"

namespace LAMMPS_NS {

FixWallDisk::FixWallDisk(LAMMPS *lmp, int narg, char **arg) : FixWall(lmp, narg, arg) {
}

FixWallDisk::~FixWallDisk() {
} 

void FixWallDisk::init()
{
  if (!atom->radius_flag) error->all(FLERR, "Fix wall/disk requires atom attribute radius");

  FixWall::init();
}

/* ---------------------------------------------------------------------- */

void FixWallDisk::precompute(int m)
{
  double sigma6 = LAMMPS_NS::MathSpecial::powint(sigma[m], 6); 
  coeff1[m] = 4.0 * epsilon[m] * sigma6;
  coeff2[m] = 4.0 * epsilon[m] * sigma6 * sigma6;
}

} 

using namespace LAMMPS_NS;
using MathSpecial::powint;

void FixWallDisk::wall_particle(int m, int which, double coord)
{

    double delta, delta2, fwall, UA, UR, Uwall;
    double rad, rad2, rad4, rad6, rad8, new_coeff1, new_coeff2, delta2_rad2;
    double vn;
    double delta4, delta6, delta8; 
    double denoA, denoR, denoR2, fwallA, fwallR; 
    const double pi = M_PI;
    double **x = atom->x;
    double **f = atom->f;
    double *radius = atom->radius;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    int dim = which / 2;
    int side = which % 2;
    double Lambda_a = 1.0; // hardcoded

    if (side == 0) side = -1;
    int onflag = 0;

    for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) {
            if (side < 0)
                delta = x[i][dim] - coord;
            else
                delta = coord - x[i][dim];

            if (delta >= cutoff[m]) continue;

            rad = radius[i];
            
            // Point-Wall Potential 

            if (rad == 0.0) {
                if (delta <= 0.0) 
                {           
                onflag = 1;                
                continue;
                }
                double delta2 = delta * delta;
                double delta4 = delta2 * delta2;
                double delta5 = delta4 * delta;      
                double delta10 = delta5 * delta5;    

                double UA = -3.0 * pi * coeff1[m] / (32.0 * delta4);
                double UR =  63.0 * pi * coeff2[m] / (2560.0 * delta10);
                Uwall = UA + UR;

                double cutoff_val  = cutoff[m];
                double cutoff2_val = cutoff_val * cutoff_val;
                double cutoff4_val = cutoff2_val * cutoff2_val;
                double cutoff10_val = cutoff4_val * cutoff4_val * cutoff2_val;

                double UA_c = -3.0 * pi * coeff1[m] / (32.0 * cutoff4_val);
                double UR_c =  63.0 * pi * coeff2[m] / (2560.0 * cutoff10_val);
                Uwall -= (UA_c + UR_c);             // shifted energy

                double FA = -3.0 * pi * coeff1[m] / ( 8.0 * delta5);
                double FR = 63.0 * pi * coeff2[m] / (256.0 * delta10 * delta);

                fwall = side * (FA + FR);

                f[i][dim] -= fwall;
                ewall[0]  += Uwall;
                ewall[m+1]+= fwall;    
                
                if (evflag) {
                    double vn = (side < 0 ? -fwall : fwall) * delta;
                    v_tally(dim, i, vn);
                }
                continue;                          
                }

            // Disk-Wall Potential 

            if (rad >= delta) { 
                onflag = 1;
                if (rad >= delta) {
                char buf[128];
                sprintf(buf,"atom %d  rad=%g  delta=%g",i,rad,delta);
                error->warning(FLERR,buf);
                onflag = 1;
                continue;
            }}
            rad2 = rad * rad;
            rad4 = rad2 * rad2;
            rad6 = rad4 * rad2;
            rad8   = rad4 * rad4; 

            delta2 = delta * delta;
            delta2_rad2 = delta2 - rad2; 
            delta4 = delta2 * delta2;
            delta6  = delta2 * delta4;
            delta8 = delta4 * delta4;  
          

            new_coeff1 = coeff1[m] * pi * pi * Lambda_a * rad2;
            new_coeff2 = coeff2[m] * pi * pi * Lambda_a * rad2;

            denoA = 32.0 * pow(delta2_rad2, 5.0 / 2.0);
            denoR = 163840.0 * pow(delta2_rad2, 17.0 / 2.0);
            denoR2 = 32768 * pow(delta2_rad2, 19.0 / 2.0);
            double P6_val_current = (35.0*rad6 + 280.0*rad4*delta2 + 336.0*rad2*delta4 + 64.0*delta6);
            UA = -(3.0 * new_coeff1 * delta) / denoA;
            UR = (63.0 * new_coeff2 * delta * P6_val_current) / denoR;
            Uwall = UA + UR; 

            // --- Calculate Potential Energy at Cutoff ---
            double UA_cutoff;
            double UR_cutoff;
            double Uwall_cutoff; 

            double cutoff_val = cutoff[m];
            double cutoff2_val = cutoff_val * cutoff_val;
            double delta2_rad2_at_cutoff = cutoff2_val - rad2;

            double denoA_at_cutoff = 64.0 * pow(delta2_rad2_at_cutoff, 5.0 / 2.0);
            double denoR_at_cutoff = 163840.0 * pow(delta2_rad2_at_cutoff, 17.0 / 2.0);
            
            double cutoff4_val = cutoff2_val * cutoff2_val;
            double cutoff6_val = cutoff4_val * cutoff2_val;
            double P6_val_at_cutoff = (35.0*rad6 + 280.0*rad4*cutoff2_val + 336.0*rad2*cutoff4_val + 64.0*cutoff6_val);

            UA_cutoff = -(3.0 * new_coeff1 * cutoff_val) / denoA_at_cutoff;
            UR_cutoff = (63.0 * new_coeff2 * cutoff_val * P6_val_at_cutoff) / denoR_at_cutoff;
            Uwall_cutoff = UA_cutoff + UR_cutoff;

            ewall[0] += (Uwall - Uwall_cutoff);
            fwallA = -(3.0 * new_coeff1 * (rad2 + 4.0 * delta2)) / (denoA * delta2_rad2); 

            double fR_num_factor = ( 7.0 * rad8 + 280.0 * rad6 * delta2 + 1120.0 * rad4 * delta4 + 896.0 * rad2 * delta6 +128 * delta8);

            fwallR = 63.0 * new_coeff2 * (fR_num_factor) / (denoR2);

            // Combine and apply force
            fwall = side * (fwallA + fwallR);
            f[i][dim] -= fwall;
            
            if (evflag) {
                if (side < 0)
                    vn = -fwall * delta;
                else
                    vn = fwall * delta;
                v_tally(dim, i, vn);
            }
        } 

    if (onflag) error->one(FLERR, "Particle on or inside fix wall surface");
}
