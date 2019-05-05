//
// Created by Tibor Poštek on 2019-03-22.
//

#ifndef ESPRESSO_LBNODES_VARIABLE_VISCOSITY_HPP
#define ESPRESSO_LBNODES_VARIABLE_VISCOSITY_HPP

#include "grid_based_algorithms/lb.hpp"
#include <iomanip>
#include <iostream>
#include <float.h>
#include <math.h>

#include "particle_data.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "cells.hpp"
#include "grid.hpp"
#include "utils/math/triangle_functions.hpp"
#include <vector>
#include <grid_based_algorithms/lb.hpp>
#include "Triangle.hpp"

using Utils::get_n_triangle;

class LBodes_variable_viscosity {
private:
    void init_data_structure();

    Vector3d A_unfolded, B_unfolded, C_unfolded;
    Vector3d A_folded, B_folded, C_folded;
    Particle *p1{nullptr}, *p2{nullptr}, *p3{nullptr};
    Bonded_ia_parameters *iaparams{nullptr};
    VarViscNode ***my_grid_var_visc;

    int size_x{0};
    int size_y{0};
    int size_z{0};

    double min_Py{DBL_MAX};
    double max_Py{DBL_MIN};
    double min_Pz{DBL_MAX};
    double max_Pz{DBL_MIN};

    void reset_algorithm_parameters();

    void check_min_max_x_y(double &min_y, double &max_y, double &min_z, double &max_z, Triangle triangle);

    void findingObjectBoundary(Triangle triangle, int pY, Vector3d normal_vector, double d,
                               std::vector<Vector3d> *boundaryPoints);

    void markingObjectBoundary(std::vector<Vector3d> &boundary_points, Vector3d normal_vector);

    void markNode(int x, int y, int z, Vector3d Z_point, Flag flag);

    void markingObjectInside(int pY, int minZ, int maxZ, int minX, int maxX);

    void remarkingObjectInside(int pY, int minZ, int maxZ, int minX, int maxX);

    VarViscNode &get_node(int x, int y, int z);

public:
    bool making_initial_algorithm{false};
    bool making_update_algorithm{false};

    void particle_from_main_loop(Triangle &unfolded_triangle, Triangle &folded_triangle, int &coutOfBPoints);

    void initial_algorithm();

    void update_algorithm();

    void print_lbnodes_variable_visc();

    void marking_object_inside();

    int coutOfMarkedNodes{0};

    int count_of_readed_nodes{0};

    int count_outer{0};
    int count_boundary{0};
    int count_inner{0};
    int count_input{0};
    int count_output{0};
    int count_input_output{0};
    int count_not_defined{0};

};


#endif //ESPRESSO_LBNODES_VARIABLE_VISCOSITY_HPP
