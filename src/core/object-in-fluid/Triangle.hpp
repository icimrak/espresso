//
// Created by Tibor Postek on 16/11/2018.
//

#ifndef ESPRESSO_TRIANGLE_HPP
#define ESPRESSO_TRIANGLE_HPP

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


#endif //ESPRESSO_TRIANGLE_HPP
