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
                lbfields[index].varViscNode.flag = Flag::outer;
            }
        }
    }
}

void LBodes_variable_viscosity::particle_from_main_loop(Particle &p) {
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
        A_unfolded = unfolded_position(*p1);
        B_unfolded = A_unfolded + get_mi_vector(p2->r.p, A_unfolded);
        C_unfolded = A_unfolded + get_mi_vector(p3->r.p, A_unfolded);
        Triangle triangle_unfolded{A_unfolded, B_unfolded, C_unfolded};

        //neviem ci folded robim spravne
        A_unfolded = folded_position(*p1);
        B_unfolded = A_unfolded + get_mi_vector(p2->r.p, A_unfolded);
        C_unfolded = A_unfolded + get_mi_vector(p3->r.p, A_unfolded);
        Triangle triangle_folded{A_unfolded, B_unfolded, C_unfolded};

        std::cout << "A " << A_unfolded[0] << "; " << A_unfolded[1] << "; " << A_unfolded[2] << std::endl;

        if (making_update_algorithm) {
            //tu zas inak zareagujem na update algorithm

        } else if (making_initial_algorithm) {
            //tu nejako zareagujem ak robim initial algorithm

            //normalovy vektor
            const Vector3d normal_vector = get_n_triangle(A_unfolded, B_unfolded, C_unfolded).normalize();
            //parameter d z rovnice roviny
            double d = -(normal_vector[0] * A_unfolded[0] + normal_vector[1] * A_unfolded[1] +
                         normal_vector[2] * A_unfolded[2]);
            int minY = ceil(std::min({A_unfolded[1], B_unfolded[1], C_unfolded[1]}));
            int maxY = floor(std::max({A_unfolded[1], B_unfolded[1], C_unfolded[1]}));

            //pre vsetky Py z toho intervalu potrebujem najst vhodne t
            std::vector<Vector3d> boundaryPoints;
            for (int pY = minY; pY <= maxY; ++pY) {
                findingObjectBoundary(triangle_folded, pY, normal_vector, d, &boundaryPoints);
                markingObjectBoundary(boundaryPoints, normal_vector);
                boundaryPoints.clear();
            }
            check_min_max_x_y(min_Py, max_Py, min_Pz, max_Pz, triangle_folded);
        }


    }
}

void LBodes_variable_viscosity::findingObjectBoundary(Triangle triangle, int pY, Vector3d normal_vector, double d,
                                                      std::vector<Vector3d> *boundaryPoints) {
    int count_NAN{0};
    //ak pretina danu usecku
    bool cuts_A = false;
    bool cuts_B = false;
    bool cuts_C = false;
    //potrebujem zistit suradnice P1 a P2
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
            //ak count_NAN==1 tak celá jedna hrana trojuholnika pretina tuto rovinu
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
            //ak count_NAN==3 tak cely trojuholnik pretina tuto rovinu
            P_points.push_back(triangle.getA());
            P_points.push_back(triangle.getB());
            P_points.push_back(triangle.getC());
            break;
        default:
            break;
    };

    //ak je priesecnik cez nejaku hranu
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
            cuts_B = true;
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
                cuts_C = true;
                P_points.push_back(Vector3d{pX, (double) pY, pZ});
            }
            goto Process_Boundary_points;
        } else if (t_AC == 0) {
            if (!cuts_A) {
                cuts_A = true;
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
            //priesecnikom su 2 body
            if (P_points.at(0)[2] > P_points.at(1)[2]) {
                z_max = P_points.at(0)[2];
                z_min = P_points.at(1)[2];
            } else {
                z_min = P_points.at(0)[2];
                z_max = P_points.at(1)[2];
            }
            //ak existuje celocislene Z
            if (floor(z_max) - ceil(z_min) >= 0) {
                int z_hlp = ceil(z_min);
                while (z_hlp <= floor(z_max)) {
                    //moze tu byt delenie nulou!!!!!
                    double x = -((normal_vector[1] * pY) + (normal_vector[2] * z_hlp) + d) / normal_vector[0];
                    boundaryPoints->push_back(Vector3d{x, (double) pY, (double) z_hlp});
                    z_hlp++;
                }
            }
            break;
        case 1:
            //priesecnikom je len jeden bod
            z_max = P_points.at(0)[2];
            z_min = P_points.at(0)[2];
            //ak existuje celocislene Z
            if (floor(z_max) - ceil(z_min) >= 0) {
                int z_hlp = ceil(z_min);
                while (z_hlp <= floor(z_max)) {
                    //moze tu byt delenie nulou!!!!!
                    double x = -((normal_vector[1] * pY) + (normal_vector[2] * z_hlp) + d) / normal_vector[0];
                    boundaryPoints->push_back(Vector3d{x, (double) pY, (double) z_hlp});
                    z_hlp++;
                }
            }
            break;
        case 3:
            //priesecnikom su 3 body....ked cely trojuholnik pretina tuto rovinu
            //Zatial to Ignorujem


            runtimeErrorMsg() << "cely trojuholnik pretina rovinu, to treba doriesit\n";


            //usporiadam si tieto body podľa suradnic, ktore prave potrebujem
            /*  sort(P_points.begin(), P_points.end(), compareByPosition_Z);
              double z_max = P_points.at(2)[2];
              double z_min = P_points.at(0)[2];

              //ak existuje celocislene Z
              if (floor(z_max) - ceil(z_min) >= 0) {
                  //prejdem celociselne body a najdem k nim X z rovnice roviny
                  int z_hlp = ceil(z_min);
                  while (z_hlp <= floor(z_max)) {
                      //skalarny vektor a moj hladany bod mi v rovnici roviny musi dat 0 => viem vypocitat bod x_folded
                      //musim zistit vsetky x_folded, ktore su v danej pY a z_hlp suradnici
                      //takze tu najdem svoje x_folded, ktore mozem odrazu aj nasekat na int pozicie
                      double x_min{};
                      double x_max{};
                      int x_hlp = ceil(x_min);
                      while (x_hlp <= floor(x_max)) {

                          x_hlp++;
                      }
                      Z_points.push_back(Vector3d{x_folded, (double) pY, (double) z_hlp});
                      z_hlp++;
                  }
              }*/
            break;
        default:
            runtimeErrorMsg() << "For pY=" << pY << " not found the intersection, but it should! " << P_points.size();
            break;
    }
}

void LBodes_variable_viscosity::markingObjectBoundary(std::vector<Vector3d> &boundary_points, Vector3d normal_vector) {
    for (size_t k = 0; k < boundary_points.size(); ++k) {
        Vector3d Z_point = boundary_points.at(k);
        int x_low = floor(Z_point[0]);
        int x_high = ceil(Z_point[0]);

        //toto by mala byt dlzka x suradnice kanala.....teda lblattice.halo_grid[0]
        if (x_high == lblattice.halo_grid[0]) {
            x_high = 0;
        }

        Vector3d low_vector{x_low - Z_point[0], 0, 0};
        Vector3d high_vector{x_high - Z_point[0], 0, 0};

        double citatel_low{normal_vector[0] * low_vector[0] + normal_vector[1] * low_vector[1] +
                           normal_vector[2] * low_vector[2]};

        double citatel_high{normal_vector[0] * high_vector[0] + normal_vector[1] * high_vector[1] +
                            normal_vector[2] * high_vector[2]};

        //len bez absolutnej hodnoty citatela to funguje!!!
        double cos_alpha_low{citatel_low /
                             sqrt(pow(normal_vector[0], 2) + pow(normal_vector[1], 2) + pow(normal_vector[2], 2)) *
                             sqrt(pow(low_vector[0], 2) + pow(low_vector[1], 2) + pow(low_vector[2], 2))};

        //len bez absolutnej hodnoty citatela to funguje!!!
        double cos_alpha_high{citatel_high /
                              sqrt(pow(normal_vector[0], 2) + pow(normal_vector[1], 2) +
                                   pow(normal_vector[2], 2)) *
                              sqrt(pow(high_vector[0], 2) + pow(high_vector[1], 2) + pow(high_vector[2], 2))};
        //nejaky if a potom mam suradnice pre y a z v Z_Points a x_folded su ako x_low alebo x_high a tieto body v meshi oflagujem

        Vector3d Z_node_temp{Z_point};
        Z_node_temp[0] = x_low;


        //Do LB nodov si budem ukladat aj Z_point bod, ktory ma preflagoval, aby som vedel rozhodnut, ze ak su
        // napr bunky strasne natesno, tak aby mi to nepreflagovali
        VarViscNode actual_node = lbfields[get_linear_index((int) Z_node_temp[0], (int) Z_node_temp[1],
                                                            (int) Z_node_temp[2], lblattice.halo_grid)].varViscNode;

        if (cos_alpha_low >= 0) {
            //  cout << "Bod " << "x_folded= " << Z_node_temp[0] << ", y= " << Z_node_temp[1] << " z= " << Z_node_temp[2]
            //       << " je vo vnutri bunky" << endl;
            if (actual_node.flag == Flag::outer &&
                (std::isnan(actual_node.Z_point[0]) || !std::isnan(actual_node.Z_point[0]) && (
                        (Z_point[0] - Z_node_temp[0]) <
                        (actual_node.Z_point[0] - Z_node_temp[0])))) {
                //ak je X rozdiel mensi, tak bude boundary, ale ak je == tak input/output
                //to by chcelo nejako doriesit
                markNode((int) Z_node_temp[0], (int) Z_node_temp[1], (int) Z_node_temp[2], Z_point, Flag::boundary_flag);
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
            //   cout << "Bod " << "x_folded= " << Z_node_temp[0] << ", y= " << Z_node_temp[1] << " z= " << Z_node_temp[2]
            //       << " je von z bunky" << endl;
            markNode((int) Z_node_temp[0], (int) Z_node_temp[1], (int) Z_node_temp[2], Z_point,
                     Flag::outer);
        }
        if (x_low != x_high) {
            Z_node_temp[0] = x_high;
            actual_node = lbfields[get_linear_index((int) Z_node_temp[0], (int) Z_node_temp[1],
                                                    (int) Z_node_temp[2], lblattice.halo_grid)].varViscNode;


            if (cos_alpha_high >= 0) {
                //    cout << "Bod " << "x_folded= " << Z_node_temp[0] << ", y= " << Z_node_temp[1] << " z= " << Z_node_temp[2]
                //        << " je vo vnutri bunky" << endl;
                if (actual_node.flag == Flag::outer && (std::isnan(actual_node.Z_point[0]) ||
                                                        !std::isnan(actual_node.Z_point[0]) && (
                                                                (Z_node_temp[0] - Z_point[0]) <
                                                                (Z_node_temp[0] -
                                                                 actual_node.Z_point[0])))) {
                    //ak je X rozdiel mensi, tak bude boundary, ale ak je == tak input/output
                    //to by chcelo nejako doriesit
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
                //    cout << "Bod " << "x_folded= " << Z_node_temp[0] << ", y= " << Z_node_temp[1] << " z= " << Z_node_temp[2]
                //         << " je von z bunky" << endl;
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


void LBodes_variable_viscosity::markNode(int x, int y, int z, Vector3d Z_point, Flag flag) {
    int index = get_linear_index(x, y, z, lblattice.halo_grid);
    //tu budem musiet nastavit aj presnu hodnotu (double)
    // lbfields[index].var_visc_gamma_shear
    lbfields[index].varViscNode = VarViscNode{Z_point, flag};
}


void LBodes_variable_viscosity::reset_algorithm_parameters() {
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
                std::cout << std::setprecision(1) << lbfields[index].varViscNode.flag;
            }
            std::cout << " -> Y-" << y << ", Z-" << z << std::endl;
        }
        std::cout << std::endl;
    }
}