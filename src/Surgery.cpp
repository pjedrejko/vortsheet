#include "Tree.hpp"
#include <algorithm>
#include "Utilities.hpp"

bool Tree::Surgery::nodeToBeRemoved(NodeIdx i){return replacedBy[i] != i;}

void Tree::Surgery::mergeNodes(NodeIdx kept, NodeIdx removed){
    //keep track of axial nodes
    if(removed == T.G.topAxis   ) T.G.topAxis    = kept;
    if(removed == T.G.bottomAxis) T.G.bottomAxis = kept;

    //self merging is fine
    replacedBy[removed] = kept;

    x.r[kept] = 0.5 * (x.r[kept] + x.r[removed]);
    x.z[kept] = 0.5 * (x.z[kept] + x.z[removed]); //probably weighted avg would be better

    double wTotal    =   w[kept] + w[removed];
    double circTotal = x.g[kept] * w[kept] + x.g[removed] * w[removed];

    T.G.W[kept] = wTotal;

    if(wTotal > 1e-14)
        x.g[kept] = circTotal / wTotal;
}

void Tree::Surgery::mergeSegments(Segment& kept, Segment& removed){
    mergeNodes(kept.start, removed.end);
    mergeNodes(kept.end,   removed.start);

    kept.ds     += removed.ds;
    kept.layers += removed.layers;

    if(kept.interfacial && removed.interfacial)
        kept.interfacial = false;
    //otherwise seg[kept].interfacial stays as it was
    //interfacial remove and non-interfacial kept cannot happen!

    kept.state = SurgeryState::ToBeKept;

    removed = Segment();//just for safety
    removed.state = SurgeryState::ToBeRemoved;

}

bool Tree::Surgery::mergingCriterion(const Segment& s12, const Segment& s34){
    //assumes head-to-tail merging

    if(s12.state != SurgeryState::Untouched || s34.state != SurgeryState::Untouched)
        return false; //do not involve seg in surgery more than once

    if(&s12 == &s34) 
        return false;

    NodeIdx i1 = s12.start;
    NodeIdx i2 = s12.end;
    NodeIdx i3 = s34.start;
    NodeIdx i4 = s34.end;

    Array2 S12{x.r[i2] - x.r[i1], x.z[i2] - x.z[i1]};
    Array2 S34{x.r[i4] - x.r[i3], x.z[i4] - x.z[i3]};
    Array2 r13{x.r[i3] - x.r[i1], x.z[i3] - x.z[i1]};
    Array2 r24{x.r[i4] - x.r[i2], x.z[i4] - x.z[i2]};
    Array2 r23{x.r[i3] - x.r[i2], x.z[i3] - x.z[i2]};
    Array2 r14{x.r[i4] - x.r[i1], x.z[i4] - x.z[i1]};

    double dN_sq = sq(Params::kN) * std::min(norm_sq(S12), norm_sq(S34));
    double dC_sq = sq(Params::kC) * std::max(norm_sq(S12), norm_sq(S34));

    bool distCriterion = 
    std::max({
    sq(cross(S12, r13)),
    sq(cross(S12, r24)),
    sq(cross(S34, r13)),
    sq(cross(S34, r24))}) < dN_sq
    &&
    std::max(norm_sq(r23), norm_sq(r14)) < dC_sq;

    bool angleCriterion = 
    dot(S12, S34) / std::sqrt(norm_sq(S12) * norm_sq(S34)) < Params::kA;

    return distCriterion && angleCriterion;
}

void Tree::Surgery::updateSegNodes(Segment& s){
    if(s.state != SurgeryState::ToBeRemoved){
        s.start = replacedBy[s.start];
        s.end   = replacedBy[s.end];
    }
}

bool Tree::Surgery::trySegSegMerge(Segment& s12, Segment& s34){
    //criterion assumes head-to-tail merging
    //non-interfacial segments don't have meaningful direction - check both cases
    //interfacial segment has to be kept

    updateSegNodes(s12);
    updateSegNodes(s34);

    if(s12.interfacial && s34.interfacial){
        if(mergingCriterion(s12, s34)){ 
            mergeSegments(s12, s34);
            return true;
        }
    }
    else if(s12.interfacial){
        if(mergingCriterion(s12, s34.flip())){
            mergeSegments(s12, s34);
            return true;
        }
    }
    else{ //c34.interfacial or none 
        if(mergingCriterion(s34, s12.flip())){
            mergeSegments(s34, s12);
            return true;
        }
    }
    return false;
}

void Tree::Surgery::initializeSurgery(){
    logg.formPrint(1, "initializing surgery {").print("\n");
    T.G.computeWeights();
    
    computeSegCenters();
    T.setWorkingArrays(segCenters, *(T.U)); //to compute centers based on X

    //initialize surgery records
    replacedBy.resize(x.size());
    for(NodeIdx i = 0; i < NodeIdx(replacedBy.size()); i++)
        replacedBy[i] = i;

    for(Segment& s: T.G.seg)
        s.state = SurgeryState::Untouched;

    T.forwardPass ();

    T.G.determineConnectivity(); 
    logg.formPrint(-1, "}").print("\n");
}

void Tree::Surgery::finalizeSurgery(){
    logg.formPrint(1, "finalizing surgery {").print("\n");

    int nSegmentsInit = T.G.seg.size();

    T.G.seg.erase(
        std::remove_if(T.G.seg.begin(), T.G.seg.end(), [](const Segment& s) {
            return s.state == SurgeryState::ToBeRemoved;
        }),
        T.G.seg.end()
    );
    int nSegmentsFinal = T.G.seg.size();
    logg.formPrint("segments removed: %-10d", nSegmentsInit - nSegmentsFinal).print("\n");

    //some segments were affected by surgeries on their neighbours
    //and were not updated due to not being tested directly
    for(Segment& s: T.G.seg)
        updateSegNodes(s);

    int nNodesInit    = T.G.x.size();
    //removing nodes -> reshuffle the dead to the end
    //need to inform segments about index change afterwards
    NodeIdx i = 0;
    NodeIdx last = x.size() - 1;
    while(i <= last){
        if(nodeToBeRemoved(i)){

            while(last > i && nodeToBeRemoved(last))
                last--;
           
            //keep track of nodes on the axis
            if(last == T.G.topAxis   ) T.G.topAxis    = i;
            if(last == T.G.bottomAxis) T.G.bottomAxis = i;

            x.r[i] = x.r[last];
            x.z[i] = x.z[last];
            x.g[i] = x.g[last];

            //for segments to read about the index chnage
            replacedBy[last] = i;
            replacedBy[i] = -1;

            last--;
            i++;
        }
        else 
            i++;
    }
    x.resize(last + 1);
    int nNodesFinal = T.G.x.size();
    
    logg.formPrint("nodes    removed: %-10d", nNodesInit - nNodesFinal).print("\n");

    //updating indices in the segments
    for(Segment& s: T.G.seg)
        updateSegNodes(s);

    T.G.computeWeights(); //to be sure
    logg.formPrint(-1, "}").print("\n");
}

void Tree::Surgery::computeSegCenters(){

    SegIdx n = T.G.seg.size();
    segCenters.resize(n);

    for(SegIdx i = 0; i < n; i++){
        NodeIdx s = T.G.seg[i].start;
        NodeIdx e = T.G.seg[i].end;
        segCenters.r[i] = 0.5 * ( x.r[s] + x.r[e] );
        segCenters.z[i] = 0.5 * ( x.z[s] + x.z[e] );
        segCenters.g[i] = 0.5 * ( x.g[s] + x.g[e] );
    }
}

bool Tree::Surgery::trySegCellMerge(SegIdx idx, CellIdx c){
    for(NodeIdx i = T.cells[c].i0; i <= T.cells[c].i1; i++){
        SegIdx j = T.node_map[i]; //center node associated with its seg
            if(trySegSegMerge(T.G.seg[idx], T.G.seg[j]))
                return true;
    }
    return false;
}

bool Tree::Surgery::trySegTreeMerge(SegIdx s, CellIdx c){
    if(T.G.seg[s].state != SurgeryState::Untouched) 
        return false; //was already merged, no more merging

    Array2 segCenter = {segCenters.r[s], segCenters.z[s]};

    double distSq = minDistSq(segCenter, T.cells[c]);
    constexpr double dlMax = Params::dl * Params::kE;

    if( distSq < sq(dlMax) ){
        if(T.cells[c].isLeaf()){
            trySegCellMerge(s, c);
        }
        else{
            for(CellIdx ch: T.cells[c].children)
                if(ch > -1) trySegTreeMerge(s, ch);
        }
    }
    return (T.G.seg[s].state != SurgeryState::Untouched); 
}

void Tree::Surgery::apply(){
    logg.formPrint(1, "surgery {").print("\n");

    initializeSurgery();

    Timer::start("detecting surgeries with DFS...");
    T.G.DFS(0, 
        [this](Connection c){
            trySegTreeMerge(c.via, 0);
        }
    );
    Timer::end();

    finalizeSurgery();
    
    logg.formPrint(-1, "}").print("\n");
}

