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
                                 std::vector<Vector3d > *boundaryPoints) {
    int count_NAN{0};
    //ak pretina danu usecku
    bool cuts_A = false;
    bool cuts_B = false;
    bool cuts_C = false;
    //potrebujem zistit suradnice P1 a P2
    std::vector<Vector3d> P_points;

    double t_AB = (pY - triangle.getA().y) / (triangle.getB().y - triangle.getA().y);
    if (triangle.getB().y - triangle.getA().y == 0.000000 && (t_AB == 0 || t_AB == 1)) {
        count_NAN++;
    }
    double t_BC = (pY - triangle.getB().y) / (triangle.getC().y - triangle.getB().y);
    if (triangle.getC().y - triangle.getB().y == 0.000000 && (t_BC == 0 || t_BC == 1)) {
        count_NAN++;
    }
    double t_AC = (pY - triangle.getA().y) / (triangle.getC().y - triangle.getA().y);
    if (triangle.getC().y - triangle.getA().y == 0.000000 && (t_AC == 0 || t_AC == 1)) {
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
        double pZ = triangle.getA().z + t_AB * (triangle.getB().z - triangle.getA().z);
        double pX = triangle.getA().x_folded + t_AB * (triangle.getB().x_folded - triangle.getA().x_folded);
        if (t_AB == 1) {
            cuts_B = true;
        } else if (t_AB == 0) {
            cuts_A = true;
        }
        P_points.push_back(Position{pX, (double) pY, pZ});
    }
    if (t_BC >= 0 && t_BC <= 1) {
        double pZ = triangle.getB().z + t_BC * (triangle.getC().z - triangle.getB().z);
        double pX = triangle.getB().x_folded + t_BC * (triangle.getC().x_folded - triangle.getB().x_folded);
        if (t_BC == 1) {
            cuts_C = true;
            P_points.push_back(Position{pX, (double) pY, pZ});
            goto Intersection_of_AC;
        } else if (t_BC == 0 && !cuts_B) {
            cuts_B = true;
            P_points.push_back(Position{pX, (double) pY, pZ});
            goto Intersection_of_AC;
        }
        P_points.push_back(Position{pX, (double) pY, pZ});
    }

    Intersection_of_AC:
    if (t_AC >= 0 && t_AC <= 1) {
        double pZ = triangle.getA().z + t_AC * (triangle.getC().z - triangle.getA().z);
        double pX = triangle.getA().x_folded + t_AC * (triangle.getC().x_folded - triangle.getA().x_folded);
        if (t_AC == 1) {
            if (!cuts_C) {
                cuts_C = true;
                P_points.push_back(Position{pX, (double) pY, pZ});
            }
            goto Process_Boundary_points;
        } else if (t_AC == 0) {
            if (!cuts_A) {
                cuts_A = true;
                P_points.push_back(Position{pX, (double) pY, pZ});
            }
            goto Process_Boundary_points;
        }
        P_points.push_back(Position{pX, (double) pY, pZ});
    }

    //Rounding Z-coordinate and compute X-coordinate
    Process_Boundary_points:
    double z_min, z_max;

    switch (P_points.size()) {
        case 2:
            //priesecnikom su 2 body
            if (P_points.at(0).z > P_points.at(1).z) {
                z_max = P_points.at(0).z;
                z_min = P_points.at(1).z;
            } else {
                z_min = P_points.at(0).z;
                z_max = P_points.at(1).z;
            }
            //ak existuje celocislene Z
            if (floor(z_max) - ceil(z_min) >= 0) {
                int z_hlp = ceil(z_min);
                while (z_hlp <= floor(z_max)) {
                    //moze tu byt delenie nulou!!!!!
                    double x = -((normal_vector.y * pY) + (normal_vector.z * z_hlp) + d) / normal_vector.x_folded;
                    boundaryPoints->push_back(Position{x, (double) pY, (double) z_hlp});
                    z_hlp++;
                }
            }
            break;
        case 1:
            //priesecnikom je len jeden bod
            z_max = P_points.at(0).z;
            z_min = P_points.at(0).z;
            //ak existuje celocislene Z
            if (floor(z_max) - ceil(z_min) >= 0) {
                int z_hlp = ceil(z_min);
                while (z_hlp <= floor(z_max)) {
                    //moze tu byt delenie nulou!!!!!
                    double x = -((normal_vector.y * pY) + (normal_vector.z * z_hlp) + d) / normal_vector.x_folded;
                    boundaryPoints->push_back(Position{x, (double) pY, (double) z_hlp});
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
              double z_max = P_points.at(2).z;
              double z_min = P_points.at(0).z;

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
                      Z_points.push_back(Position{x_folded, (double) pY, (double) z_hlp});
                      z_hlp++;
                  }
              }*/
            break;
        default:
            runtimeErrorMsg() << "For pY=" << pY << " not found the intersection, but it should! " << P_points.size();
            break;
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
                std::cout << std::setprecision(1) << lbfields[index].flag;
            }
            std::cout << " -> Y-" << y << ", Z-" << z << std::endl;
        }
        std::cout << std::endl;
    }
}