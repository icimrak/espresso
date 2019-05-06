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
    int size_x{0};
    int size_y{0};
    int size_z{0};

    double min_Py = std::numeric_limits<double>::max();
    double max_Py = std::numeric_limits<double>::min();
    double min_Pz = std::numeric_limits<double>::max();
    double max_Pz = std::numeric_limits<double>::min();

    double inner_fluid_visc{0};

    void reset_algorithm_parameters();

    void check_min_max_x_y(double &min_y, double &max_y, double &min_z, double &max_z, Triangle triangle);

    void findingObjectBoundary(Triangle triangle, int pY, Vector3d normal_vector, double d,
                               std::vector<Vector3d> *boundaryPoints);

    void markingObjectBoundary(std::vector<Vector3d> &boundary_points, Vector3d normal_vector);

    void markingObjectBoundary_update_algorithm(std::vector<Vector3d> &boundary_points, Vector3d normal_vector);

    void markNode(int x, int y, int z, Vector3d Z_point, Flag flag);

    void markingObjectInside(int pY, int minZ, int maxZ, int minX, int maxX);

    void remarkingObjectInside(int pY, int minZ, int maxZ, int minX, int maxX);

    LB_FluidNode &get_node(int x, int y, int z);

    void initial_algorithm(Triangle &unfolded_triangle, Triangle &folded_triangle);

    void update_algorithm(Triangle &unfolded_triangle, Triangle &folded_triangle);

    void set_viscosity_to_node(bool is_inner, LB_FluidNode& node);

public:
    bool making_initial_algorithm{false};
    bool making_update_algorithm{false};

    void particle_from_main_loop(Triangle &unfolded_triangle, Triangle &folded_triangle, double inner_fluid_visc);

    void before_initial_algorithm();

    void before_update_algorithm();

    void print_lbnodes_variable_visc();

    void print_lbnodes_variable_visc(int y);

    void marking_object_inside();
};


#endif //ESPRESSO_LBNODES_VARIABLE_VISCOSITY_HPP
