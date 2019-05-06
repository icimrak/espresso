//
// Created by Tibor Poštek on 2019-03-22.
//

#include "lbnodes_variable_viscosity.hpp"

void LBodes_variable_viscosity::init_data_structure() {
    size_x = lblattice.halo_grid[0];  // size of simulation box should looks like (int) box_l[0];
    size_y = lblattice.halo_grid[1];  //(int) box_l[1];
    size_z = lblattice.halo_grid[2];  //(int) box_l[2];

    for (int x = 0; x < size_x; ++x) {
        for (int y = 0; y < size_y; ++y) {
            for (int z = 0; z < size_z; ++z) {
                LB_FluidNode &node = get_node(x, y, z);
                // initializes the variable viscosity fields, all the fields will be constant with viscosity values given by lbfluid.
                set_viscosity_to_node(false, node);
                node.varViscNode.flag = Flag::outer;
                node.varViscNode.Z_point = Vector3d{std::numeric_limits<double>::quiet_NaN(),
                                                    std::numeric_limits<double>::quiet_NaN(),
                                                    std::numeric_limits<double>::quiet_NaN()};
            }
        }
    }
}


void LBodes_variable_viscosity::particle_from_main_loop(Triangle &unfolded_triangle, Triangle &folded_triangle,
                                                        double inner_fluid_visc) {
    this->inner_fluid_visc = inner_fluid_visc;
    if (making_initial_algorithm || making_update_algorithm) {
        if (making_update_algorithm) {
            //section of update algorithm
            update_algorithm(unfolded_triangle, folded_triangle);
        } else {
            //section of initial algorithm
            initial_algorithm(unfolded_triangle, folded_triangle);
        }
    }
}

void LBodes_variable_viscosity::initial_algorithm(Triangle &unfolded_triangle, Triangle &folded_triangle) {
    const Vector3d normal_vector = get_n_triangle(unfolded_triangle.getB(), unfolded_triangle.getC(),
                                                  unfolded_triangle.getA()).normalize();
    //parameter d from equation plane
    double d = -(normal_vector[0] * unfolded_triangle.getA()[0] +
                 normal_vector[1] * unfolded_triangle.getA()[1] +
                 normal_vector[2] * unfolded_triangle.getA()[2]);
    int minY = static_cast<int>(ceil(
            std::min({unfolded_triangle.getA()[1], unfolded_triangle.getB()[1], unfolded_triangle.getC()[1]})));
    int maxY = static_cast<int>(floor(
            std::max({unfolded_triangle.getA()[1], unfolded_triangle.getB()[1], unfolded_triangle.getC()[1]})));

    std::vector<Vector3d> boundaryPoints;
    for (int pY = minY; pY <= maxY; ++pY) {
        //maybe there should be folded_triangle
        findingObjectBoundary(unfolded_triangle, pY, normal_vector, d, &boundaryPoints);
        markingObjectBoundary(boundaryPoints, normal_vector);
        boundaryPoints.clear();
    }
    //maybe it should be folded_triangle instead of unfolded_triangle
    check_min_max_x_y(reinterpret_cast<double &>(min_Py), max_Py, min_Pz, max_Pz, unfolded_triangle);
}


void LBodes_variable_viscosity::update_algorithm(Triangle &unfolded_triangle, Triangle &folded_triangle) {

}


void LBodes_variable_viscosity::marking_object_inside() {
    for (int pY = static_cast<int>(min_Py); pY < max_Py; ++pY) {
        markingObjectInside(pY, static_cast<int>(min_Pz), static_cast<int>(max_Pz), 0, size_x);
    }
    for (int pY = static_cast<int>(min_Py); pY < max_Py; ++pY) {
        remarkingObjectInside(pY, static_cast<int>(min_Pz), static_cast<int>(max_Pz), 0, size_x);
    }
}

void LBodes_variable_viscosity::markingObjectInside(int pY, int minZ, int maxZ, int minX, int maxX) {
    for (int pZ = minZ; pZ <= maxZ; ++pZ) {
        int BN{0};
        for (int pX = minX; pX < maxX; ++pX) {
            VarViscNode U = get_node(pX, pY, pZ).varViscNode;
            if (U.flag == Flag::boundary_flag && BN % 2 == 0) {
                markNode(pX, pY, pZ, Vector3d{(double) pX, (double) pY, (double) pZ}, Flag::input);
                BN++;
            } else if (U.flag != Flag::boundary_flag && BN % 2 != 0) {
                markNode(pX, pY, pZ, Vector3d{(double) pX, (double) pY, (double) pZ}, Flag::inner);
            } else if (U.flag == Flag::boundary_flag && BN % 2 != 0) {
                markNode(pX, pY, pZ, Vector3d{(double) pX, (double) pY, (double) pZ}, Flag::output);
                BN++;
            } else if (U.flag == Flag::input_output && BN % 2 == 0) {
                BN += 2;
            }
            /*TODO
             * This part, when node is not_defined is not implemented yet
             * when whole area from triangle is crossed by plane
             */
        }
    }
}

void LBodes_variable_viscosity::remarkingObjectInside(int pY, int minZ, int maxZ, int minX, int maxX) {
    for (int pZ = minZ; pZ <= maxZ; ++pZ) {
        for (int pX = minX; pX < maxX; ++pX) {
            VarViscNode U = get_node(pX, pY, pZ).varViscNode;
            if (U.flag != Flag::outer) {
                markNode(pX, pY, pZ,
                         Vector3d{std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(),
                                  std::numeric_limits<double>::quiet_NaN()}, Flag::inner);
            } else {
                markNode(pX, pY, pZ,
                         Vector3d{std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(),
                                  std::numeric_limits<double>::quiet_NaN()}, Flag::outer);
            }
        }
    }
}

void LBodes_variable_viscosity::findingObjectBoundary(Triangle triangle, int pY, Vector3d normal_vector, double d,
                                                      std::vector<Vector3d> *boundaryPoints) {
    int count_NAN{0};
    //variables, that indicate that segment is crossed or cut
    bool cuts_A = false;
    bool cuts_B = false;
    bool cuts_C = false;

    //I need to find out position P1 and P2
    std::vector<Vector3d> P_points;

    double t_AB = (pY - triangle.getA()[1]) / (triangle.getB()[1] - triangle.getA()[1]);
    if (triangle.getB()[1] - triangle.getA()[1] == 0.000000 && (t_AB == 0 || t_AB == 1)) {
        count_NAN++;
    }
    double t_BC = (pY - triangle.getB()[1]) / (triangle.getC()[1] - triangle.getB()[1]);
    if (triangle.getC()[1] - triangle.getB()[1] == 0.000000 && (t_BC == 0 || t_BC == 1)) {
        count_NAN++;
    }
    double t_AC = (pY - triangle.getA()[1]) / (triangle.getC()[1] - triangle.getA()[1]);
    if (triangle.getC()[1] - triangle.getA()[1] == 0.000000 && (t_AC == 0 || t_AC == 1)) {
        count_NAN++;
    }

    switch (count_NAN) {
        case 1:
            //When count_NAN==1 then only one edge of triangle cross the current plane pY
            if (t_AB == 0) {
                P_points.push_back(triangle.getA());
            } else if (t_AB == 1) {
                P_points.push_back(triangle.getB());
            }
            if (t_BC == 0) {
                P_points.push_back(triangle.getB());
            } else if (t_BC == 1) {
                P_points.push_back(triangle.getC());
            }
            if (t_AC == 0) {
                P_points.push_back(triangle.getA());
            } else if (t_AC == 1) {
                P_points.push_back(triangle.getC());
            }
            goto Process_Boundary_points;
        case 2:
            runtimeErrorMsg() << "Error, because 2 edges cannot lie on that plane\n";
            break;
        case 3:
            //When count_NAN==1 then whole triangle cross the current plane pY
            P_points.push_back(triangle.getA());
            P_points.push_back(triangle.getB());
            P_points.push_back(triangle.getC());
            break;
        default:
            break;
    };

    //If is any edge cut or crossed
    if (t_AB >= 0 && t_AB <= 1) {
        double pZ = triangle.getA()[2] + t_AB * (triangle.getB()[2] - triangle.getA()[2]);
        double pX = triangle.getA()[0] + t_AB * (triangle.getB()[0] - triangle.getA()[0]);
        if (t_AB == 1) {
            cuts_B = true;
        } else if (t_AB == 0) {
            cuts_A = true;
        }
        P_points.push_back(Vector3d{pX, (double) pY, pZ});
    }
    if (t_BC >= 0 && t_BC <= 1) {
        double pZ = triangle.getB()[2] + t_BC * (triangle.getC()[2] - triangle.getB()[2]);
        double pX = triangle.getB()[0] + t_BC * (triangle.getC()[0] - triangle.getB()[0]);
        if (t_BC == 1) {
            cuts_C = true;
            P_points.push_back(Vector3d{pX, (double) pY, pZ});
            goto Intersection_of_AC;
        } else if (t_BC == 0 && !cuts_B) {
            P_points.push_back(Vector3d{pX, (double) pY, pZ});
            goto Intersection_of_AC;
        }
        P_points.push_back(Vector3d{pX, (double) pY, pZ});
    }

    Intersection_of_AC:
    if (t_AC >= 0 && t_AC <= 1) {
        double pZ = triangle.getA()[2] + t_AC * (triangle.getC()[2] - triangle.getA()[2]);
        double pX = triangle.getA()[0] + t_AC * (triangle.getC()[0] - triangle.getA()[0]);
        if (t_AC == 1) {
            if (!cuts_C) {
                P_points.push_back(Vector3d{pX, (double) pY, pZ});
            }
            goto Process_Boundary_points;
        } else if (t_AC == 0) {
            if (!cuts_A) {
                P_points.push_back(Vector3d{pX, (double) pY, pZ});
            }
            goto Process_Boundary_points;
        }
        P_points.push_back(Vector3d{pX, (double) pY, pZ});
    }

    //Rounding Z-coordinate and compute X-coordinate
    Process_Boundary_points:
    double z_min, z_max;

    switch (P_points.size()) {
        case 2:
            //The intersection consists from 2 points
            if (P_points.at(0)[2] > P_points.at(1)[2]) {
                z_max = P_points.at(0)[2];
                z_min = P_points.at(1)[2];
            } else {
                z_min = P_points.at(0)[2];
                z_max = P_points.at(1)[2];
            }
            //If Z is integer
            if (floor(z_max) - ceil(z_min) >= 0) {
                int z_hlp = static_cast<int>(ceil(z_min));
                while (z_hlp <= floor(z_max)) {
                    //TODO Here should be zero division!!!!!
                    double x = -((normal_vector[1] * pY) + (normal_vector[2] * z_hlp) + d) / normal_vector[0];
                    boundaryPoints->push_back(Vector3d{x, (double) pY, (double) z_hlp});
                    z_hlp++;
                }
            }
            break;
        case 1:
            //The intersection consists from 1 point only
            z_max = P_points.at(0)[2];
            z_min = P_points.at(0)[2];
            //If Z is integer
            if (floor(z_max) - ceil(z_min) >= 0) {
                int z_hlp = static_cast<int>(ceil(z_min));
                while (z_hlp <= floor(z_max)) {
                    //TODO Here should be zero division!!!!!
                    double x = -((normal_vector[1] * pY) + (normal_vector[2] * z_hlp) + d) / normal_vector[0];
                    boundaryPoints->push_back(Vector3d{x, (double) pY, (double) z_hlp});
                    z_hlp++;
                }
            }
            break;
        case 3:
            /*TODO
             * The intersection consists from 3 points
             * This part, when node is not_defined is not implemented yet
             * when whole area from triangle is crossed by plane
             */
            runtimeErrorMsg() << "The whole area from triangle is crossed by plane. It is not implemented yet.\n";
            break;
        default:
            runtimeErrorMsg() << "For pY=" << pY << " not found the intersection, but it should! " << P_points.size();
            break;
    }
}

void LBodes_variable_viscosity::markingObjectBoundary(std::vector<Vector3d> &boundary_points, Vector3d normal_vector) {
    for (Vector3d Z_point : boundary_points) {
        int x_low = static_cast<int>(floor(Z_point[0]));
        int x_high = static_cast<int>(ceil(Z_point[0]));

        if (x_high == size_x) {
            x_high = 0;
        }

        Vector3d low_vector{x_low - Z_point[0], 0, 0};
        Vector3d high_vector{x_high - Z_point[0], 0, 0};

        double numerator_low{normal_vector[0] * low_vector[0] + normal_vector[1] * low_vector[1] +
                             normal_vector[2] * low_vector[2]};
        double numerator_high{normal_vector[0] * high_vector[0] + normal_vector[1] * high_vector[1] +
                              normal_vector[2] * high_vector[2]};

        //Without absolute value of numerator it is correct!!
        double cos_alpha_low{numerator_low /
                             sqrt(pow(normal_vector[0], 2) + pow(normal_vector[1], 2) + pow(normal_vector[2], 2)) *
                             sqrt(pow(low_vector[0], 2) + pow(low_vector[1], 2) + pow(low_vector[2], 2))};
        double cos_alpha_high{numerator_high /
                              sqrt(pow(normal_vector[0], 2) + pow(normal_vector[1], 2) +
                                   pow(normal_vector[2], 2)) *
                              sqrt(pow(high_vector[0], 2) + pow(high_vector[1], 2) + pow(high_vector[2], 2))};

        Vector3d Z_node_temp{Z_point};
        Z_node_temp[0] = x_low;

        VarViscNode actual_node = get_node((int) Z_node_temp[0], (int) Z_node_temp[1],
                                           (int) Z_node_temp[2]).varViscNode;

        if (cos_alpha_low >= 0) {
            //Z_node_temp is inside of immersed object
            if (actual_node.flag == Flag::outer &&
                (std::isnan(actual_node.Z_point[0]) || (!std::isnan(actual_node.Z_point[0]) && (
                        (Z_point[0] - Z_node_temp[0]) <
                        (actual_node.Z_point[0] - Z_node_temp[0]))))) {
                /*TODO Input_output marking
                 * ak je X rozdiel mensi, tak bude boundary, ale ak je == tak input/output
                 * to by chcelo nejako doriesit
                */
                markNode((int) Z_node_temp[0], (int) Z_node_temp[1], (int) Z_node_temp[2], Z_point,
                         Flag::boundary_flag);
            } else if (actual_node.flag == Flag::boundary_flag) {
                if (actual_node.Z_point[0] != Z_point[0]) {
                    markNode((int) Z_node_temp[0], (int) Z_node_temp[1], (int) Z_node_temp[2], Z_point,
                             Flag::input_output);
                } else {
                    markNode((int) Z_node_temp[0], (int) Z_node_temp[1], (int) Z_node_temp[2], Z_point,
                             Flag::boundary_flag);
                }
            } else {
                markNode((int) Z_node_temp[0], (int) Z_node_temp[1], (int) Z_node_temp[2], Z_point,
                         Flag::outer);
            }
        } else {
            //Z_node_temp is outside of immersed object
            markNode((int) Z_node_temp[0], (int) Z_node_temp[1], (int) Z_node_temp[2], Z_point,
                     Flag::outer);
        }
        if (x_low != x_high) {
            Z_node_temp[0] = x_high;
            actual_node = get_node((int) Z_node_temp[0], (int) Z_node_temp[1], (int) Z_node_temp[2]).varViscNode;

            if (cos_alpha_high >= 0) {
                //Z_node_temp is inside of immersed object
                if (actual_node.flag == Flag::outer &&
                    (std::isnan(actual_node.Z_point[0]) || (!std::isnan(actual_node.Z_point[0]) && (
                            (Z_node_temp[0] - Z_point[0]) <
                            (Z_node_temp[0] - actual_node.Z_point[0]))))) {
                    /*TODO Input_output marking
                     * ak je X rozdiel mensi, tak bude boundary, ale ak je == tak input/output
                     * to by chcelo nejako doriesit
                    */
                    markNode((int) Z_node_temp[0], (int) Z_node_temp[1], (int) Z_node_temp[2], Z_point,
                             Flag::boundary_flag);
                } else if (actual_node.flag == Flag::boundary_flag) {
                    if (actual_node.Z_point[0] != Z_point[0]) {
                        markNode((int) Z_node_temp[0], (int) Z_node_temp[1], (int) Z_node_temp[2], Z_point,
                                 Flag::input_output);
                    } else {
                        markNode((int) Z_node_temp[0], (int) Z_node_temp[1], (int) Z_node_temp[2], Z_point,
                                 Flag::boundary_flag);
                    }
                } else {
                    markNode((int) Z_node_temp[0], (int) Z_node_temp[1], (int) Z_node_temp[2], Z_point,
                             Flag::outer);
                }
            } else {
                //Z_node_temp is outside of immersed object
                markNode((int) Z_node_temp[0], (int) Z_node_temp[1], (int) Z_node_temp[2], Z_point,
                         Flag::outer);
            }
        }
    }
}


void LBodes_variable_viscosity::check_min_max_x_y(double &min_y, double &max_y, double &min_z, double &max_z,
                                                  Triangle triangle) {
    if (triangle.getA()[1] < min_y) {
        min_y = triangle.getA()[1];
    }
    if (triangle.getB()[1] < min_y) {
        min_y = triangle.getB()[1];
    }
    if (triangle.getC()[1] < min_y) {
        min_y = triangle.getC()[1];
    }
    if (triangle.getA()[1] > max_y) {
        max_y = triangle.getA()[1];
    }
    if (triangle.getB()[1] > max_y) {
        max_y = triangle.getB()[1];
    }
    if (triangle.getC()[1] > max_y) {
        max_y = triangle.getC()[1];
    }

    if (triangle.getA()[1] < min_z) {
        min_z = triangle.getA()[1];
    }
    if (triangle.getB()[1] < min_z) {
        min_z = triangle.getB()[1];
    }
    if (triangle.getC()[1] < min_z) {
        min_z = triangle.getC()[1];
    }
    if (triangle.getA()[1] > max_z) {
        max_z = triangle.getA()[1];
    }
    if (triangle.getB()[1] > max_z) {
        max_z = triangle.getB()[1];
    }
    if (triangle.getC()[1] > max_z) {
        max_z = triangle.getC()[1];
    }
}


LB_FluidNode &LBodes_variable_viscosity::get_node(int x, int y, int z) {
    int index = get_linear_index(x, y, z, lblattice.halo_grid);
    return lbfields[index];
}

void LBodes_variable_viscosity::markNode(int x, int y, int z, Vector3d Z_point, Flag flag) {
    LB_FluidNode &node = get_node(x, y, z);

    //I will also store Z_point, because I can see which Particle marked my given node
    node.varViscNode = VarViscNode{Z_point, flag};
    set_viscosity_to_node(flag == Flag::inner, node);
}


void LBodes_variable_viscosity::reset_algorithm_parameters() {
    min_Py = std::numeric_limits<double>::max();
    max_Py = std::numeric_limits<double>::min();
    min_Pz = std::numeric_limits<double>::max();
    max_Pz = std::numeric_limits<double>::min();
}

void LBodes_variable_viscosity::before_initial_algorithm() {
    init_data_structure();
    reset_algorithm_parameters();
}

void LBodes_variable_viscosity::before_update_algorithm() {

    //  print_lbnodes_variable_visc();
}

void LBodes_variable_viscosity::print_lbnodes_variable_visc() {
    for (int y = 0; y < size_y; ++y) {
        for (int z = 0; z < size_z; ++z) {
            for (int x = 0; x < size_x; ++x) {
                std::cout << std::setw(5) << std::setprecision(1) << get_node(x, y, z).var_visc_gamma_shear << " ";
            }
            std::cout << " -> Y-" <<std::setw(3)<< y << ", Z-" <<std::setw(3)<< z << " ...... ";
            for (int x = 0; x < size_x; ++x) {
                std::cout << get_node(x, y, z).varViscNode.flag;
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}

void LBodes_variable_viscosity::set_viscosity_to_node(bool is_inner, LB_FluidNode &node) {
    if (is_inner) {
        node.var_visc_gamma_shear = 1. - 2. / (6. * inner_fluid_visc * lbpar.tau / (lbpar.agrid * lbpar.agrid) + 1.);
    } else {
        node.var_visc_gamma_shear = 1. - 2. / (6. * lbpar.viscosity * lbpar.tau / (lbpar.agrid * lbpar.agrid) + 1.);
    }
}
