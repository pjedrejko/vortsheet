#include "Tree.hpp"
#include "Utilities.hpp"
#include<stdexcept>
#include<iostream>

void Tree::setWorkingArrays(const Nodes& xWorking, Nodes& uWorking){
    //these doesnt have to be actuall sheets nodes!
    //can be temporary array (rk4) or segCenters (surgery)
    this->X = &xWorking;
    this->U = &uWorking;
}

void Tree::forwardPass(){
    //building the tree by sorting nodes indices in node_map
    //and assigning each cell a contiguous index range

    Timer::start("building tree     (forward pass)... ");

    //initialize identity map
    node_map.resize(X->size());
    for(NodeIdx i = 0; i < NodeIdx(node_map.size()); i++)
        node_map[i] = i;

    cells.clear();
    leaves.clear();
    addCell(-1, 0, X->size() - 1); //root

    int lvlStart = 0;
    int lvlEnd = cells.size();

    while(lvlStart < lvlEnd){
        nLvls++;

        for(CellIdx idx = lvlStart; idx < lvlEnd; idx++){

            if(isSmallEnough(idx))
                leaves.push_back(idx);
            else 
                splitCell(idx);
        }

        lvlStart = lvlEnd;
        lvlEnd = cells.size();
    }

    Timer::end();
}

NodeIdx Tree::partition(const std::vector<double>& y, double yCenter, NodeIdx i0, NodeIdx i1){
    //puts indices of the nodes with coordinate < center of the cell to the left,
    //larger coords to the right (not sorted within the sides). Returns last index of the
    //left side

    NodeIdx iLeft  = i0;
    NodeIdx iRight = i1;

    while(iLeft <= iRight){
        if( y[node_map[iLeft]] > yCenter ){
            std::swap(node_map[iLeft], node_map[iRight]);
            iRight--;
        }
        else
            iLeft++;
    }

    NodeIdx i1Left = iLeft - 1;
    return i1Left;
}

CellIdx Tree::addCell(CellIdx parent, NodeIdx i0, NodeIdx i1){

    if( i0 > i1 ) return -1;

    cells.push_back(Cell(parent, i0, i1));
    Cell& C = cells.back();

    //determining boundaries (for tighter bound on epsilon)
    for(NodeIdx i = i0; i <= i1; i++){
        C.rLim[0] = std::min( X->r[node_map[i]], C.rLim[0]);
        C.rLim[1] = std::max( X->r[node_map[i]], C.rLim[1]);

        C.zLim[0] = std::min( X->z[node_map[i]], C.zLim[0]);
        C.zLim[1] = std::max( X->z[node_map[i]], C.zLim[1]);
    }

    //computing centroid
    C.xC[0] = ( C.rLim[0] + C.rLim[1] ) * 0.5;
    C.xC[1] = ( C.zLim[0] + C.zLim[1] ) * 0.5;

    return cells.size() - 1;
}

void Tree::splitCell(CellIdx idx){
    //split nodes by coordinates, first in r (left or right from the cells center),
    //then each part split analogousely in z
    //  [parent]   =>   [01|11]  \ children
    //  [ cell ]        [00|10]  /

    const double rCenter = cells[idx].xC[0];
    const double zCenter = cells[idx].xC[1];

    NodeIdx i0_00 = cells[idx].i0;
    NodeIdx i1_11 = cells[idx].i1;

    //split in r [00, 01 | 10, 11]
    NodeIdx i1_01 = partition(X->r, rCenter, i0_00, i1_11);
    NodeIdx i0_10 = i1_01 + 1;

    //split in z [00| 01 | 10, 11]
    NodeIdx i1_00 = partition(X->z, zCenter, i0_00, i1_01);
    NodeIdx i0_01 = i1_00 + 1;

    //split in z [00| 01 | 10| 11]
    NodeIdx i1_10 = partition(X->z, zCenter, i0_10, i1_11);
    NodeIdx i0_11 = i1_10 + 1;

    NodeIdx c00 = addCell(idx, i0_00, i1_00);
    NodeIdx c01 = addCell(idx, i0_01, i1_01);
    NodeIdx c10 = addCell(idx, i0_10, i1_10);
    NodeIdx c11 = addCell(idx, i0_11, i1_11);

    //(in 2 steps cause addCell can reallocate cells)
    cells[idx].children[0] = c00;
    cells[idx].children[1] = c01;
    cells[idx].children[2] = c10;
    cells[idx].children[3] = c11;
}

bool Tree::isSmallEnough(CellIdx idx){
    return cells[idx].nNodes() <= Params::nC;
};

void Tree::toFile(int step, const std::string& prefix) const {
    std::string filename = "results/" + prefix + std::to_string(step);
    std::printf("writing %-*s\t", Params::w - 7, (filename + "...").c_str());
    
    FILE* f = std::fopen(filename.c_str(), "w");
    if (!f) 
        throw std::runtime_error("Could not open file: " + filename);
    
    std::fprintf(f, "%% rLim[0] rLim[1] zLim[0] zLim[1] isLeaf\n");
    for (const Cell& c: cells) {
        std::fprintf(f, "%.17g %.17g %.17g %.17g %d\n", c.rLim[0], c.rLim[1], c.zLim[0], c.zLim[1], c.isLeaf());
    }

    std::fclose(f);
    std::printf("done.\n");
}
