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
#ifndef _OBJECT_IN_FLUID_OIF_GLOBAL_FORCES_H
#define _OBJECT_IN_FLUID_OIF_GLOBAL_FORCES_H
/** \file
 *  Routines to calculate the OIF_GLOBAL_FORCES energy or/and and force
 *  for a particle triple (triangle from mesh). (Dupin2007)
 *  \ref forces.cpp
 */

#include "particle_data.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"
#include "grid_based_algorithms/lb_interface.hpp"
#include "grid_based_algorithms/lb.hpp"
#include "particle_data.hpp"

#include <utils/Vector.hpp>
#include <utils/math/triangle_functions.hpp>
using Utils::angle_btw_triangles;
using Utils::area_triangle;
using Utils::get_n_triangle;
#include <utils/constants.hpp>

#include "lbnodes_variable_viscosity.hpp"

/** set parameters for the OIF_GLOBAL_FORCES potential.
 */
int oif_global_forces_set_params(int bond_type, double A0_g, double ka_g,
                                 double V0, double kv, double inner_fluid_visc);
void calc_oif_global(double *area_volume, int molType);
void add_oif_global_forces(double *area_volume, int molType);

bool calc_vectors_of_triangles(Particle &p, Utils::Vector3d &p11, Utils::Vector3d &p22, Utils::Vector3d &p33, Particle* p1, Particle* p2,
        Particle* p3, int molType, Bonded_ia_parameters* iaparams, int test);


#ifdef LB_VARIABLE_VISCOSITY
extern LBodes_variable_viscosity *lbodes_variable_viscosity;
extern bool flagging_lbnodes_var_visc;
extern bool reflagging_lbnodes_var_visc;


void flag_lbnodes_variable_visc(LBodes_variable_viscosity *lbodes_variable_visc);
void reflag_lbnodes_variable_visc();
void update_flags_variable_visc();
#endif

/************************************************************/

extern int max_oif_objects;
#endif
