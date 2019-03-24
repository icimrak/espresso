//
// Created by Tibor Postek on 16/11/2018.
//

#include "Triangle.hpp"


Triangle::Triangle(Vector3d A, Vector3d B, Vector3d C)
        : A{A}, B{B}, C{C} {
}

Vector3d Triangle::getA() {
    return this->A;
}

Vector3d Triangle::getB() {
    return this->B;
}

Vector3d Triangle::getC() {
    return this->C;
}
