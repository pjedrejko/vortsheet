#include "Sheet.hpp"
#include "Tree.hpp"
#include "Velocity.hpp"
#include "Params.hpp"
#include "Utilities.hpp"
#include<stdexcept>

void Sheet::clear() {
    x.clear();
    u.clear();
}

void Sheet::initializeRing(double totalCirculation) {
    clear();

    logg.formPrint(1, "Initializing vortex ring {")                  .print("\n");
    logg.formPrint(   "initial circulation: %lf", totalCirculation)  .print("\n");
    logg.formPrint(   "initial length:      %lf", Params::lenInitial).print("\n");
    logg.formPrint(   "initial nodes:       %d" , Params::nInitial)  .print("\n");
    logg.formPrint(   "segment length:      %lf", Params::dl)        .print("\n");

    const double R = Params::lenInitial / M_PI;
    const double ds = Params::dl; //parametrize with initial physical length

    G.bottomAxis = 0;
    G.topAxis = Params::nInitial - 1;

    NodeIdx i0 = -1;
    NodeIdx i1 = -1;

    for (int i = 0; i < Params::nInitial; i++) {
        double s = -0.5 * M_PI + i * ds;
        double r = std::cos(s) * R;
        double z = std::sin(s) * R;
        double g = std::cos(s) * 0.5 * totalCirculation;

        if(i == G.bottomAxis || i == G.topAxis) r = 0.0; //fix on z-axis

        i1 = x.add(r, z, g);
        if(i0 != -1) G.addSegment(i0, i1, ds, true);
        i0 = i1;
    }

    G.computeWeights();
    logg.formPrint(-1, "}\n");
}


void Sheet::computeVelocities(Nodes& X, Nodes& U){
    //RK4 can pass some temporary arrays for X and U

    parallel.brodcastNodes(X);
    U.assign(X.size(), 0.0);
    
    T.setWorkingArrays(X, U);
    
    T.forwardPass ();
    T.velocity.backwardPass();

    Timer::start("computing velocities... ");
    parallel.determineLeafRange();

    #pragma omp parallel for
    for(LeafIdx l = parallel.l0; l <= parallel.l1; l++){
        T.velocity.induce(T.leaves[l], 0);
    }
    Timer::end();

    parallel.reduceVelocities(U);
}

void Sheet::produceCirculation(const Nodes& x, Nodes& dxdt){
    Timer::start("computing production... ");

    for(const Segment& s: G.seg){
        if(s.ds < 1e-14) continue;

        double dgdt = Params::B * ( x.z[s.end] - x.z[s.start] ) / s.ds;
        dxdt.g[s.start] += 0.5 * dgdt;
        dxdt.g[s.end  ] += 0.5 * dgdt;
        
        //double dGdt = Params::B * ( x.z[s.end] - x.z[s.start] );
        //dxdt.g[s.start] += 0.5 * dGdt / G.W[s.start];
        //dxdt.g[s.end  ] += 0.5 * dGdt / G.W[s.end];
    }

    Timer::end();
}

void Sheet::toFile(){
    logg.formPrint(1, "writing sheet {").print("\n");

    int stp = integrals.step;
    integrals.toFile();
    x.toFile(stp, "x");
    G.toFile(stp);
    logg.formPrint(-1, "}\n");
}

void Sheet::fromFile(int stp){
    logg.formPrint(1, "reading sheet from step %d {", stp).print("\n");

    integrals.fromFile(stp);
    x.fromFile(stp, "x");
    G.fromFile(stp);
    logg.formPrint(-1, "}\n");
}
