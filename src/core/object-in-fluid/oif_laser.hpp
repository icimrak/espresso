/*
  Copyright (C) 2012-2018 The ESPResSo project

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef _OBJECT_IN_FLUID_OIF_LASER_H
#define _OBJECT_IN_FLUID_OIF_LASER_H

/** \file
 *  Routines to calculate the OIF_LASER forces
 *  for a particle triangle (one triangle).
 * \ref forces.cpp
 */

#include "bonded_interactions/bonded_interaction_data.hpp"
#include "grid.hpp"
#include "particle_data.hpp"
#include "utils/Vector.hpp"
#include "utils/math/triangle_functions.hpp"

// set parameters for local forces
int oif_laser_set_params(int bond_type, double klas, double lasX,
                                double lasY, double lasZ);

/** Computes the local laser forces (Dupin2007) and adds them
    to the particle forces.
    @param p1,p2,p3     Pointers to particles of triangle.
    @return 0
*/
inline int calc_oif_laser(Particle *p1, Particle *p2, Particle *p3,
                          Bonded_ia_parameters *iaparams,
                          double force[3], double force2[3],
                          double force3[3]) 
{

  auto const fp1 = unfolded_position(*p1);
  auto const fp2 = fp1 + get_mi_vector(p2->r.p, fp1);
  auto const fp3 = fp1 + get_mi_vector(p3->r.p, fp1);

  for (int i = 0; i < 3; i++) {
    force[i] = 0;
    force2[i] = 0;
    force3[i] = 0;
  }

  // non-linear stretching
  if (iaparams->p.oif_laser.klas > TINY_OIF_ELASTICITY_COEFFICIENT) {
    //TODO: implementation of laser forces
    // pristup k parametrom:
    // iaparams->p.oif_laser.klas;
    // iaparams->p.oif_laser.lasX;
    
    
    //auto const dx = fp2 - fp3;
    //auto const len = dx.norm();
    //auto const dr = len - iaparams->p.oif_local_forces.r0;
    //auto const lambda = 1.0 * len / iaparams->p.oif_local_forces.r0;
    //auto const fac =
        //-iaparams->p.oif_local_forces.ks * KS(lambda) * dr; // no normalization
    //for (int i = 0; i < 3; i++) {
      //force2[i] += fac * dx[i] / len;
      //force3[i] += -fac * dx[i] / len;
    //}
  }
  return 0;
}

#endif
