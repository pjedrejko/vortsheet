#pragma once
#include <vector>
#include <cmath>
#include "Nodes.hpp"
#include "Tree.hpp"
#include "Graph.hpp"

class Sheet {
public:
    Nodes x; //poistion
    Nodes u; //velocity
    std::vector<NodeWeight> W; //for ds integration

    Graph G; //Sheets structure (segments, connections etc.)
    Tree T;  //for NlogN operations (velocity, surgery)

    class RungeKutta4{
        Sheet& sheet;
        Nodes k1, k2, k3, k4;
        Nodes temp_x;

        void resize(size_t count);
        void computeRHS(Nodes& x, Nodes& dxdt);
        void euler(Nodes& out, const Nodes& x, double dt, const Nodes& u);
        
        public:
        RungeKutta4(Sheet& parent): sheet{parent}{};
        void stepEuler(double dt);
        void step(double dt);
    };
    RungeKutta4 rk4;

    class Integrals{
        const Sheet& sheet;
        public:
        int    step            = 0; 
        double time            = 0.0;
        double dt              = 0.0;
        double circulation     = 0.0;
        double axisThickness   = 0.0;
        double impulse         = 0.0;
        double volume          = 0.0;
        double parameterLength = 0.0;
        double physicalLength  = 0.0;
   
        Integrals(const Sheet& parent): sheet{parent}{};
        void compute();
        void toFile();
        void fromFile(int stp);
        void info();
        
        double computeCirculation();
        double computeAxisThickness();
        double computeImpulse();
        double computeVolume();
        double computeParameterLength();
        double computePhysicalLength();
        double computeTimeStep(double circulation, double axisThickness);
    }; 
    Integrals integrals;

    class Sync{
        public:
        Sheet& sheet;
        //int iRank =  -1;
        //int nRanks = -1;
        //int nThreads = -1;

        //leaves range (in the Tree)
        LeafIdx l0 = -1;
        LeafIdx l1 = -1;
        Sync(Sheet& sheet): sheet(sheet){};

        void info();
        void brodcastNodes(Nodes& X);
        void determineLeafRange();
        void reduceVelocities(Nodes& U);
    };
    Sync parallel;

    Sheet(): G(x, W), T(G), rk4(*this), integrals(*this), parallel(*this){};
    void clear();

    void initializeRing(double totalCirculation);

    void fromFile(int stp);
    void toFile();

    void brodcastNodes();
    void reduceVelocities();

    void computeVelocities(Nodes& X, Nodes& U);
    void produceCirculation(const Nodes& x, Nodes& dxdt);
};
