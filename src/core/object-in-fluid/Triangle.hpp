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
    Utils::Vector3d A, B, C;
public:
    Triangle(Utils::Vector3d A, Utils::Vector3d B, Utils::Vector3d C);

    Utils::Vector3d getA();

    Utils::Vector3d getB();

    Utils::Vector3d getC();

};


#endif //ESPRESSO_TRIANGLE_HPP
