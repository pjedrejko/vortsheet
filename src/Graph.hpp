#pragma once
#include<vector>
#include"Nodes.hpp"
#include"Array.hpp"

/*
stores info about Sheets structure: 
segments, their connections and node weights (computed from segments);
keeps track of two nodes on the z-axis to keep the axisymmetric sheet closed
*/

using SegIdx = int;

enum class SurgeryState{
    ToBeKept,
    ToBeRemoved,
    Untouched
};

struct Connection{
    NodeIdx to;
    SegIdx via;
    bool visited; //for DFS
    Connection(NodeIdx to, SegIdx via, bool visited = false): 
    to{to}, via{via}, visited{visited}{};
};

struct Segment {
    NodeIdx start = -1;           //index of the starting node
    NodeIdx end = -1;             //index of the ending node
    double ds = 0.0;              //parameter length
    bool interfacial = false;     //if separates buoyant region from non-buoyant
    int layers = 1;               //how many layers (>1 after surgery)

    SurgeryState state = SurgeryState::Untouched;
    Segment() = default;
    Segment(int s, int e, double ds, bool inter, int layers = 1): 
    start{s}, end{e}, ds{ds}, interfacial{inter}, layers{layers}, state{SurgeryState::Untouched} {}
    Segment& flip();
};

class Graph{
    public:
    Nodes& x;
    std::vector<NodeWeight>& W;

    NodeIdx topAxis = -1; //these nodes must stay on the z-axis
    NodeIdx bottomAxis = -1;

    std::vector<Segment> seg;
    std::vector<std::vector<Connection>> C; //for DFS

    Graph(Nodes& x, std::vector<NodeWeight>& W): x(x), W(W){};

    void determineConnectivity();
    void computeWeights();
    void refine();

    double lengthSq(SegIdx idx);
    SegIdx addSegment(NodeIdx start, NodeIdx end, double ds, bool interfacial);
    SegIdx splitSegment(SegIdx idx);

    void toFile(int step) const;
    void fromFile(int step);

    template<typename Lambda>
    void DFS(NodeIdx start, Lambda actionOnConnection) {
        std::vector<NodeIdx> stack;
        stack.push_back(start);

        while (!stack.empty()) {
            NodeIdx current = stack.back();
            stack.pop_back();

            for (auto& c : C[current]) {
                if (!c.visited) {
                    c.visited = true;

                    actionOnConnection(c);

                    stack.push_back(c.to);
                }
            }
        }
    }

};