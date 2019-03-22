//
// Created by Tibor Poštek on 2019-03-22.
//

#include "lbnodes_variable_viscosity.hpp"

void init_data_structure(){
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








void var_visc_initial_algorithm(){
    init_data_structure();
    double min_Py{DBL_MAX};
    double max_Py{DBL_MIN};
    double min_Pz{DBL_MAX};
    double max_Pz{DBL_MIN};


/** initializes the variable viscosity fields, all the fields will be constant with viscosity values given by lbfluid. */
    for (int x = 0; x < lblattice.halo_grid[0]; ++x) {
        for (int y = 0; y < lblattice.halo_grid[1]; ++y) {
            for (int z = 0; z < lblattice.halo_grid[2]; ++z) {
                int index = get_linear_index(x, y, z, lblattice.halo_grid);

                //outer viscosity of cell
                   lbfields[index].var_visc_gamma_shear = 1. - 2. /
                           (6. * lbpar.viscosity * lbpar.tau /(lbpar.agrid * lbpar.agrid) + 1.);


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

void var_visc_update_algorithm(){





  //  print_lbnodes_variable_visc();
}

void print_lbnodes_variable_visc(){
    //printing all LB nodes with viscosity
    for (int y = 0; y < lblattice.halo_grid[1]; ++y) {
        for (int z = 0; z < lblattice.halo_grid[2]; ++z) {
            for (int x = 0; x < lblattice.halo_grid[0]; ++x) {
                int index = get_linear_index(x, y, z, lblattice.halo_grid);
                std::cout << std::setprecision (1) << lbfields[index].flag;
            }
            std::cout << " -> Y-" << y << ", Z-" << z << std::endl;
        }
        std::cout << std::endl;
    }
}