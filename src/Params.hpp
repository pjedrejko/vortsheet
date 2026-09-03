#pragma once
#include <cmath>

namespace Params {
    //names ~ as in paper

    //physics
    static constexpr double B = 1.0;        //buoyancy
    static constexpr double delta = 0.01;   //smoothing

    //fast velocity computation tuning
    static constexpr double epsilon = 0.05; //tree algorithm's error accepted
    static constexpr int nC = 4;            //nodes per leaf cell

    //discretization tuning
    static constexpr double kL = 0.1;       //~ dl = kL * delta
    static constexpr double lenInitial = M_PI;
    static constexpr int    nInitial   = int((lenInitial / (delta * kL))) + 1;
    static constexpr double dl         = lenInitial / (nInitial - 1);
    static constexpr double kE = 4.0/3;     //split segment when longer than kE * dl 

    //adaptive time step tuning
    static constexpr double kt1 = 0.1;
    static constexpr double kt2 = 2.0;

    //surgery tuning
    static constexpr double kN = 0.2; 
    static constexpr double kC = 0.65;
    static constexpr double kA = -0.5;//cos(120 deg)  //-0.98480775301; //cos(170 deg)

    
    static constexpr int w = 61; //width of info printed
}
