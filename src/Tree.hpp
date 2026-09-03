#pragma once
#include "Cell.hpp"
#include <vector>
#include <string>
#include <cmath>
#include "Nodes.hpp"
#include "Velocity.hpp"
#include "Params.hpp"
#include "Graph.hpp"

/*
accelerates N^2 operations to NlogN (velocity computation and surgery)
*/

using LeafIdx = int;

class Tree {
public:
    std::vector<Cell> cells;
    std::vector<NodeIdx> node_map; //for sorting nodes (->populating cells); less efficient but easier with rk4;
                                   //accessing coords in a cell always use x.r[node_map[i]] not x.r[i]!
    std::vector<CellIdx> leaves;

    const Nodes* X = nullptr; // pointer cause can be some temporary positions (rk4) or segment centers (surgery)
    Nodes*       U = nullptr; //

    Graph& G; // G.x are always the actual positions of the sheet nodes

    int nLvls = 0;

    struct Velocity{
        Tree& T;
        Velocity(Tree& T): T(T){};
        void backwardPass();
        void computeAllMoments(CellIdx idx);
        void computeMoments   (CellIdx idx);
        void biotSavartApprox (CellIdx probe, CellIdx source);
        void biotSavartExact  (CellIdx probe, CellIdx source);
        double estimateEpsilon(CellIdx probe, CellIdx source);
        void induce           (CellIdx probe, CellIdx source);
    };
    Velocity velocity;

    struct Surgery{
        Tree& T;

        Nodes& x; //the actual sheet nodes (from G); for convenience
        std::vector<NodeWeight>& w;

        std::vector<NodeIdx> replacedBy; //to inform segments about node mergers
        Nodes segCenters;                //to group segments in a tree

        Surgery(Tree& T): T(T), x(T.G.x), w(T.G.W){};

        void mergeNodes(NodeIdx kept, NodeIdx removed);
        void mergeSegments(Segment& kept, Segment& removed);
        bool mergingCriterion(const Segment& s12, const Segment& s34);

        bool nodeToBeRemoved(NodeIdx i);
        void updateSegNodes(Segment& s);
        
        bool trySegSegMerge(Segment& s12, Segment& s34);
        bool trySegCellMerge(SegIdx idx, CellIdx c);
        bool trySegTreeMerge(SegIdx s, CellIdx c);
        
        void computeSegCenters();
        
        void initializeSurgery();
        void apply();
        void finalizeSurgery();
    };
    Surgery surgery;

    Tree(Graph& G): G(G), velocity(*this), surgery(*this){};

    void setWorkingArrays(const Nodes& xWorking, Nodes& uWorking);

    //building tree methods:
    void forwardPass();
    bool isSmallEnough(CellIdx idx);
    void splitCell(CellIdx idx);
    NodeIdx partition(const std::vector<double>& y, double yCenter, NodeIdx i0, NodeIdx i1);
    CellIdx addCell(CellIdx parent, NodeIdx i0, NodeIdx i1);

    void toFile(int step, const std::string& prefix) const;
};
