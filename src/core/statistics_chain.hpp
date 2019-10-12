/*
 * Copyright (C) 2010-2019 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef STATISTICS_CHAIN_H
#define STATISTICS_CHAIN_H
/** \file
 *
 *  This file contains the code for statistics on the data using the
 *  molecule information set with analyze set chains.
 */

#include "PartCfg.hpp"

/** \name Exported Variables */
/************************************************************/
/** Particles' initial positions (needed for g1(t), g2(t), g3(t)) */
/*@{*/
extern float *partCoord_g;
extern float *partCM_g;
extern int n_part_g;
extern int n_chains_g;
/*@}*/

/** data for a system consisting of chains. TBRS. */
/*@{*/
extern int chain_start;
extern int chain_n_chains;
extern int chain_length;
/*@}*/

/** \name Exported Functions */
/************************************************************/
/*@{*/

/** Calculate the end-to-end-distance.
 *  Chain information \ref chain_start etc. must be set!
 */
std::array<double, 4> calc_re(PartCfg &);

/** Calculate the radius of gyration.
 *  Chain information \ref chain_start etc. must be set!
 */
std::array<double, 4> calc_rg(PartCfg &);

/** Calculate the hydrodynamic radius (ref. Kirkwood-Zimm theory).
 *  Chain information \ref chain_start etc. must be set!
 */
std::array<double, 2> calc_rh(PartCfg &);

/*@}*/

#endif
