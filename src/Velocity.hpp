#pragma once
#include <vector>
#include "Array.hpp"

namespace point2point{

    struct Prereq{
        double ksq, E, K;
        double R1sq, R2sq, R2cb;
        double I1, I2;

        Prereq(Array2 x0, Array2 x, double delta);
    };

    Array2 biotSavart(const Prereq& p, Array2 x0, Array2 x, double S);

    Array2x2 uGradEval(const Prereq& p, Array2 x0C, Array2 xP, double SP);

    Array2 biotSavartApprox(const Array2& uC, const Array2x2& uGrad, const Array2& Dx);
}


