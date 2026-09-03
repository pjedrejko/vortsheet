#pragma once
#include <array>
#include <cmath>
#include "Array.hpp"

using CellIdx = int;

struct Cell{

    //nodes range
    int i0 = -1; //first node
    int i1 = -1; //last node

    //position in tree
    CellIdx parent = -1;
    std::array<CellIdx, 4> children{-1, -1, -1, -1};

    //boundaries
    Array2 rLim{{INFINITY, -INFINITY}}, zLim{{INFINITY, -INFINITY}};

    //barycenters' positions:
    Array2 xP = {0, 0}; //positive barycenter
    Array2 xM = {0, 0}; //negative barycenter

    //barycenter's strengths:
    double SP = 0.0;
    double SM = 0.0;
   
    //centroid
    Array2 xC = {0, 0};

    Cell(CellIdx parent, int i0, int i1): i0{i0}, i1{i1}, parent{parent}
    {};

    double h(){
        return std::max(rLim[1] - rLim[0], zLim[1] - zLim[0]);
    }

    double diag() const {
        return std::sqrt( sq(rLim[1] - rLim[0]) + sq(zLim[1] - zLim[0]) );
    }

    bool isLeaf() const{
        return
        ( nNodes() > 0 ) &&
        children[0] == -1 &&
        children[1] == -1 &&
        children[2] == -1 &&
        children[3] == -1;
    }

    int nNodes() const {
        if ( i0 < 0 || i1 < 0 )
            return 0;
        else
        return (i1 - i0) + 1;
    }

};

double minDistSq(const Cell& c1, const Cell& c2);
double minDistSq(const Array2& xLim, const Cell& c2);