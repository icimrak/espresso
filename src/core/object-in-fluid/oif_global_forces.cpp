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
 *  Routines to calculate the OIF_GLOBAL_FORCES energy or/and and force
 *  for a particle triple (triangle from mesh). (Dupin2007)
 *  \ref forces.cpp
 */

#include "oif_global_forces.hpp"
#include "../../../build/myconfig.hpp"

/** set parameters for the OIF_GLOBAL_FORCES potential.
 */
int oif_global_forces_set_params(int bond_type, double A0_g, double ka_g,
                                 double V0, double kv, double inner_fluid_visc) {
    if (bond_type < 0)
        return ES_ERROR;

    make_bond_type_exist(bond_type);

    bonded_ia_params[bond_type].p.oif_global_forces.ka_g = ka_g;
    bonded_ia_params[bond_type].p.oif_global_forces.A0_g = A0_g;
    bonded_ia_params[bond_type].p.oif_global_forces.V0 = V0;
    bonded_ia_params[bond_type].p.oif_global_forces.kv = kv;
    bonded_ia_params[bond_type].p.oif_global_forces.inner_fluid_visc = inner_fluid_visc;

    bonded_ia_params[bond_type].type = BONDED_IA_OIF_GLOBAL_FORCES;
    bonded_ia_params[bond_type].num = 2;

    /* broadcast interaction parameters */
    mpi_bcast_ia_params(bond_type, -1);

    return ES_OK;
}

/** called in force_calc() from within forces.cpp
 *  calculates the global area and global volume for a cell before the forces
 * are handled
 *  sums up parts for area with mpi_reduce from local triangles
 *  synchronization with allreduce
 *
 *  !!! loop over particles from domain_decomposition !!!
 */

#ifdef LB_VARIABLE_VISCOSITY
LBodes_variable_viscosity *lbodes_variable_viscosity{nullptr};
bool flagging_lbnodes_var_visc{false};
bool reflagging_lbnodes_var_visc{false};

void flag_lbnodes_variable_visc(LBodes_variable_viscosity *lbodes_variable_visc) {
    lbodes_variable_viscosity = lbodes_variable_visc;
    lbodes_variable_viscosity->initial_algorithm();
    flagging_lbnodes_var_visc = true;
    update_flags_variable_visc();
}


void reflag_lbnodes_variable_visc() {
    lbodes_variable_viscosity->update_algorithm();
    reflagging_lbnodes_var_visc = true;
    flagging_lbnodes_var_visc = false;
    update_flags_variable_visc();
}

void update_flags_variable_visc(){
    lbodes_variable_viscosity->making_initial_algorithm = flagging_lbnodes_var_visc;
    lbodes_variable_viscosity->making_update_algorithm = reflagging_lbnodes_var_visc;
}

#endif


void calc_oif_global(double *area_volume, int molType) { // first-fold-then-the-same approach
    double partArea = 0.0;
    double part_area_volume[2]; // added

    // z volume
    double VOL_partVol = 0.;

    int test = 0;

    for (auto &p : local_cells.particles()) {
#ifdef LB_VARIABLE_VISCOSITY
        lbodes_variable_viscosity->particle_from_main_loop(p);
#endif
        Vector3d p11, p22, p33;
        Particle *p1{nullptr}, *p2{nullptr}, *p3{nullptr};
        Bonded_ia_parameters *iaparams{nullptr};
        if (calc_vectors_of_triangles(p, p11, p22, p33, p1, p2, p3, molType, iaparams, test)) {
            // unfolded positions correct
            auto const VOL_A = area_triangle(p11, p22, p33);
            partArea += VOL_A;

            auto const VOL_norm = get_n_triangle(p11, p22, p33);
            auto const VOL_dn = VOL_norm.norm();
            auto const VOL_hz = 1.0 / 3.0 * (p11[2] + p22[2] + p33[2]);
            VOL_partVol += VOL_A * -1 * VOL_norm[2] / VOL_dn * VOL_hz;
        }
    }

    part_area_volume[0] = partArea;
    part_area_volume[1] = VOL_partVol;

    MPI_Allreduce(part_area_volume, area_volume, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#ifdef LB_VARIABLE_VISCOSITY
    //odtialto zavolam marking object inside ak prebieha init algoritmus
    if(flagging_lbnodes_var_visc){
        lbodes_variable_viscosity->marking_object_inside();
    }
    reflagging_lbnodes_var_visc = false;
    flagging_lbnodes_var_visc = false;
    update_flags_variable_visc();
#endif
}

void add_oif_global_forces(double *area_volume, int molType) { // first-fold-then-the-same approach
    double area = area_volume[0];
    double VOL_volume = area_volume[1];

    int test = 0;


    for (auto &p : local_cells.particles()) {
#ifdef LB_VARIABLE_VISCOSITY
        //lbodes_variable_viscosity->particle_from_main_loop(p);
#endif
        Vector3d p11, p22, p33;
        Particle *p1{nullptr}, *p2{nullptr}, *p3{nullptr};
        Bonded_ia_parameters *iaparams{nullptr};
        if (calc_vectors_of_triangles(p, p11, p22, p33, p1, p2, p3, molType, iaparams, test)) {
            // unfolded positions correct
            /// starting code from volume force
            auto const VOL_norm = get_n_triangle(p11, p22, p33).normalize();
            auto const VOL_A = area_triangle(p11, p22, p33);
            auto const VOL_vv = (VOL_volume - iaparams->p.oif_global_forces.V0) / iaparams->p.oif_global_forces.V0;
            auto const VOL_force = (1.0 / 3.0) * iaparams->p.oif_global_forces.kv * VOL_vv * VOL_A * VOL_norm;
            p1->f.f += VOL_force;
            p2->f.f += VOL_force;
            p3->f.f += VOL_force;
            ///  ending code from volume force

            auto const h = (1. / 3.) * (p11 + p22 + p33);
            auto const deltaA = (area - iaparams->p.oif_global_forces.A0_g) / iaparams->p.oif_global_forces.A0_g;
            auto const m1 = h - p11;
            auto const m2 = h - p22;
            auto const m3 = h - p33;

            auto const m1_length = m1.norm();
            auto const m2_length = m2.norm();
            auto const m3_length = m3.norm();

            auto const fac = iaparams->p.oif_global_forces.ka_g * VOL_A * deltaA /
                             (m1_length * m1_length + m2_length * m2_length + m3_length * m3_length);

            p1->f.f += fac * m1;
            p2->f.f += fac * m2;
            p3->f.f += fac * m3;
        }
    }
#ifdef LB_VARIABLE_VISCOSITY
    //odtialto zavolam marking object inside ak prebieha init algoritmus
    if(flagging_lbnodes_var_visc){
        lbodes_variable_viscosity->marking_object_inside();
    }
    reflagging_lbnodes_var_visc = false;
    flagging_lbnodes_var_visc = false;
    update_flags_variable_visc();
#endif
}

bool calc_vectors_of_triangles(Particle &p, Vector3d &p11, Vector3d &p22, Vector3d &p33, Particle *p1, Particle *p2,
                               Particle *p3, int molType, Bonded_ia_parameters *iaparams, int test) {
    int j = 0;
    p1 = &p;
    while (j < p1->bl.n) {
        // bond type
        auto const type_num = p1->bl.e[j++];
        iaparams = &bonded_ia_params[type_num];
        auto const type = iaparams->type;
        auto const n_partners = iaparams->num;
        auto const id = p1->p.mol_id;
        if (type == BONDED_IA_OIF_GLOBAL_FORCES && id == molType) { // BONDED_IA_OIF_GLOBAL_FORCES with correct molType
            test++;
            // fetch particle 2
            p2 = local_particles[p1->bl.e[j++]];
            if (!p2) {
                runtimeErrorMsg() << "add area: bond broken between particles "
                                  << p1->p.identity << " and " << p1->bl.e[j - 1]
                                  << " (particles not stored on the same node - "
                                     "oif_globalforce2); n "
                                  << p1->bl.n << " max " << p1->bl.max;
                return false;
            }
            // fetch particle 3
            // if(n_partners>2){
            p3 = local_particles[p1->bl.e[j++]];
            if (!p3) {
                runtimeErrorMsg()
                        << "add area: bond broken between particles " << p1->p.identity
                        << ", " << p1->bl.e[j - 2] << " and " << p1->bl.e[j - 1]
                        << " (particles not stored on the same node); n " << p1->bl.n
                        << " max " << p1->bl.max;
                return false;
            }

            p11 = unfolded_position(*p1);
            p22 = p11 + get_mi_vector(p2->r.p, p11);
            p33 = p11 + get_mi_vector(p3->r.p, p11);
            printf("naslo mi to \n");
            return true;
        } else {
            j += n_partners;
        }
    }
    return false;
}


int max_oif_objects = 0;
bool oif_objects_up_to_date = false;
