//
// Created by Tibor Poštek on 2019-03-22.
//

#ifndef ESPRESSO_LBNODES_VARIABLE_VISCOSITY_HPP
#define ESPRESSO_LBNODES_VARIABLE_VISCOSITY_HPP

#include "grid_based_algorithms/lb.hpp"
#include <iomanip>
#include <iostream>
#include <float.h>

#include "particle_data.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "cells.hpp"
#include "grid.hpp"

class LBodes_variable_viscosity{
private:
    void init_data_structure();

public:
    bool looping_over_particles{false};

    bool making_initial_algorithm{false};

    void particle_from_main_loop(Particle &p);

    void var_visc_initial_algorithm();

    void var_visc_update_algorithm();

    void print_lbnodes_variable_visc();

};



#endif //ESPRESSO_LBNODES_VARIABLE_VISCOSITY_HPP
