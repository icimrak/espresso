//
// Created by Tibor Poštek on 2019-03-22.
//

#include "lbnodes_variable_viscosity.hpp"

void LBodes_variable_viscosity::init_data_structure() {
    for (int x = 0; x < lblattice.halo_grid[0]; ++x) {
        for (int y = 0; y < lblattice.halo_grid[1]; ++y) {
            for (int z = 0; z < lblattice.halo_grid[2]; ++z) {
                int index = get_linear_index(x, y, z, lblattice.halo_grid);
                lbfields[index].var_visc_gamma_shear = std::numeric_limits<double>::quiet_NaN();
                lbfields[index].flag = Flag::outer;
            }
        }
    }
}

void LBodes_variable_viscosity::particle_from_main_loop(Particle &p){
    if (making_initial_algorithm || making_update_algorithm) {
        int j = 0;
        p1 = &p;

        // bond type
        auto const type_num = p1->bl.e[j++];
        iaparams = &bonded_ia_params[type_num];
        auto const type = iaparams->type;
        auto const n_partners = iaparams->num;
        auto const id = p1->p.mol_id;

        // fetch particle 2
        p2 = local_particles[p1->bl.e[j++]];
        if (!p2) {
            runtimeErrorMsg() << "add area: bond broken between particles "
                              << p1->p.identity << " and " << p1->bl.e[j - 1]
                              << " (particles not stored on the same node - "
                                 "oif_globalforce2); n "
                              << p1->bl.n << " max " << p1->bl.max;
            return;
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
            return;
        }


        //budem potrebovat aj folded aj unfolded
        A = unfolded_position(*p1);
        B = A + get_mi_vector(p2->r.p, A);
        C = A + get_mi_vector(p3->r.p, A);

        std::cout << "A " << A.m_storage.m_data[0] << "; " << A.m_storage.m_data[1] << "; "
                  << A.m_storage.m_data[2] << std::endl;

        if (making_update_algorithm){
            //tu zas inak zareagujem na update algorithm

        } else if (making_initial_algorithm){
            //tu nejako zareagujem ak robim initial algorithm

            //normalovy vektor
            const Vector3d normal_vector = get_n_triangle(A, B, C).normalize();
            //parameter d z rovnice roviny
            double d = -(normal_vector[0] * A[0] + normal_vector[1] * A[1] + normal_vector[2] * A[2]);
            int minY = ceil(std::min({A[1],B[1],C[1]}));
            int maxY = floor(std::max({A[1], B[1], C[1]}));
        }








    }
}


void LBodes_variable_viscosity::reset_algorithm_parameters(){
    min_Py = DBL_MAX;
    max_Py = DBL_MIN;
    min_Pz = DBL_MAX;
    max_Pz = DBL_MIN;
}

void LBodes_variable_viscosity::initial_algorithm() {
    init_data_structure();
    reset_algorithm_parameters();




// initializes the variable viscosity fields, all the fields will be constant with viscosity values given by lbfluid.
    for (int x = 0; x < lblattice.halo_grid[0]; ++x) {
        for (int y = 0; y < lblattice.halo_grid[1]; ++y) {
            for (int z = 0; z < lblattice.halo_grid[2]; ++z) {
                int index = get_linear_index(x, y, z, lblattice.halo_grid);

                //outer viscosity of cell
                lbfields[index].var_visc_gamma_shear = 1. - 2. /
                                                            (6. * lbpar.viscosity * lbpar.tau /
                                                             (lbpar.agrid * lbpar.agrid) + 1.);


                //inner viscosity of cell
                //  lbfields[index].var_visc_gamma_shear = 1. - 2. /
                //         (6. * inner_fluid_visc * lbpar.tau /(lbpar.agrid * lbpar.agrid) + 1.);


                //Currently not defined variable
                // lbfields[index].var_visc_gamma_bulk = 1. - 2. /
                //         (9. * lbpar.bulk_viscosity * lbpar.tau / (lbpar.agrid * lbpar.agrid) + 1.);

            }
        }
    }
    //  print_lbnodes_variable_visc();
    // NEXT, we continue with reflagging over all cells. TODO.
}

void LBodes_variable_viscosity::update_algorithm() {





    //  print_lbnodes_variable_visc();
}

void LBodes_variable_viscosity::print_lbnodes_variable_visc() {
    //printing all LB nodes with viscosity
    for (int y = 0; y < lblattice.halo_grid[1]; ++y) {
        for (int z = 0; z < lblattice.halo_grid[2]; ++z) {
            for (int x = 0; x < lblattice.halo_grid[0]; ++x) {
                int index = get_linear_index(x, y, z, lblattice.halo_grid);
                std::cout << std::setprecision(1) << lbfields[index].flag;
            }
            std::cout << " -> Y-" << y << ", Z-" << z << std::endl;
        }
        std::cout << std::endl;
    }
}