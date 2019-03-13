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

/** \file
 *  Routines to calculate the OIF_LASER
 *  for a particle triangle.
 * ) \ref forces.cpp
 */

#include "oif_laser.hpp"
#include "communication.hpp"

/** set parameters for the OIF_LASER potential.*/

int oif_laser_set_params(int bond_type, double klas, double lasX,
                                double lasY, double lasZ) {
  if (bond_type < 0)
    return ES_ERROR;

  make_bond_type_exist(bond_type);

  bonded_ia_params[bond_type].p.oif_laser.klas = klas;
  bonded_ia_params[bond_type].p.oif_laser.lasX = lasX;
  bonded_ia_params[bond_type].p.oif_laser.lasY = lasY;
  bonded_ia_params[bond_type].p.oif_laser.lasZ = lasZ;

  bonded_ia_params[bond_type].type = BONDED_IA_OIF_LASER;
  bonded_ia_params[bond_type].num = 2;

  mpi_bcast_ia_params(bond_type, -1);

  return ES_OK;
}
