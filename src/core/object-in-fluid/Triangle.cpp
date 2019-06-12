//
// Created by Tibor Postek on 16/11/2018.
//

#include "Triangle.hpp"


Triangle::Triangle(Utils::Vector3d A, Utils::Vector3d B, Utils::Vector3d C)
        : A{A}, B{B}, C{C} {
}

Utils::Vector3d Triangle::getA() {
    return this->A;
}

Utils::Vector3d Triangle::getB() {
    return this->B;
}

Utils::Vector3d Triangle::getC() {
    return this->C;
}
