/*
  Copyright (C) 2010-2018 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

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
#ifndef _VISCOUS_HPP
#define _VISCOUS_HPP
/** \file
 *  Routines to calculate the VISOUC force
 *  for a particle pair.
 *  \ref forces.cpp
 */

/************************************************************/

#include "bonded_interaction_data.hpp"
#include "particle_data.hpp"
#include "utils.hpp"
#include "grid.hpp"

/// set the parameters for the harmonic potential
int viscous_set_params(int bond_type, double k, double r, double r_cut);

/** Computes the VISCOUS pair force and adds this
    force to the particle forces (see \ref bonded_interaction_data.cpp).
    @param p1        Pointer to first particle.
    @param p2        Pointer to second/middle particle.
    @param iaparams  bond type number of the angle interaction (see \ref
   bonded_interaction_data.cpp).
    @param dx        particle distance vector
    @param force     returns force of particle 1
    @return 0.
*/
inline int calc_viscous_pair_force(Particle *p1, Particle *p2,
                                    Bonded_ia_parameters *iaparams,
                                    double dx[3], double force[3]) {
  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);

  if ((iaparams->p.viscous.r_cut > 0.0) && (dist > iaparams->p.viscous.r_cut))
    return 1;
  // viscous force
    auto const fp1 = unfolded_position(*p1);
    auto const fp2 = fp1 + get_mi_vector(p2->r.p, fp1); 
     
    auto const udx = fp1 - fp2;
    auto const len2 = udx.norm2();
    auto const v_ij = p2->m.v - p1->m.v;

    // Variant A
    // Here the force is in the direction of relative velocity btw points
    if (iaparams->p.viscous.r > 0) { 
    // Code:
        for(int i=0; i<3; i++) {
            force[i] = - iaparams->p.viscous.k*v_ij[i];
            force[i] =  iaparams->p.viscous.k*v_ij[i];
        }
    }

    // Variant B
    // Here the force is the projection of relative velocity btw points onto
    // line btw the points

    // denote p vector between p2 and p3
    // denote v the velocity difference between the points p2 and p3
    // denote alpha the angle between p and v
    // denote x the projevted v onto p
    // cos alpha = |x|/|v|
    // cos alpha = (v,p)/(|v||p|)
    // together we get |x|=(v,p)/|p|
    // also, x is along p, so x = |x|.p/|p|
    // so x = p/|p| . (v,p)/|p|
    // altogether x = p . (v,p)/(|p|)^2
    // |p|^2 is stored in len2

    // Code:
    if (iaparams->p.viscous.r < 0) { 
        auto const fac = iaparams->p.viscous.k * (udx * v_ij) / len2;

        for (int i = 0; i < 3; i++) {
          force[i] = fac * udx[i];
          force[i] = fac * udx[i];
        }
    }
  


  return 0;
}

#endif
