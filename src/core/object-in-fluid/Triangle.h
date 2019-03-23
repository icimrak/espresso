//
// Created by Tibor Postek on 16/11/2018.
//

#ifndef CELL_INSIDE_ALG_TRIANGLE_H
#define CELL_INSIDE_ALG_TRIANGLE_H

#include <ostream>
#include <cmath>
#include "utils/Vector.hpp"

class Triangle {
private:
    Vector3d A, B, C;
public:
    Triangle(Vector3d A, Vector3d B, Vector3d C);

    Vector3d getA();

    Vector3d getB();

    Vector3d getC();

};


#endif //CELL_INSIDE_ALG_TRIANGLE_H
