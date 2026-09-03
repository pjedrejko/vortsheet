#include"Graph.hpp"
#include"Params.hpp"
#include<stdexcept>
#include <algorithm>
#include"Utilities.hpp"

Segment& Segment::flip(){
    //useful cause non-interfacial segments dont have a well defined direction
    if(interfacial)
        throw std::runtime_error("Could not flip interfacial segment");
    else
        std::swap(start, end);

    return (*this);
}

void Graph::computeWeights(){
    logg.formPrint("updating weights... ");

    W.assign(x.size(), 0.0);

    for(const Segment& s: seg){
        W[s.start] += 0.5 * s.ds;
        W[s.end  ] += 0.5 * s.ds;
    }
    logg.print("done.\n");
}

void Graph::refine(){
    logg.formPrint(1, "refining discretization {").print("\n");

    int nodesAdded = 0;
    SegIdx i = 0;
    while(i < SegIdx(seg.size())){
        if( lengthSq(i) > sq(Params::kE * Params::dl)){
            splitSegment(i);
            nodesAdded++;
        }
        else i++;
    }
    logg.formPrint("nodes added: %-10d\n", nodesAdded);
    computeWeights();
    logg.print(-1, "}\n");
}

SegIdx Graph::addSegment(NodeIdx start, NodeIdx end, double ds, bool interfacial){
    seg.push_back(Segment(start, end, ds, interfacial));
    return seg.size() - 1;
}

SegIdx Graph::splitSegment(SegIdx idx){
    NodeIdx s = seg[idx].start;
    NodeIdx e = seg[idx].end;

    NodeIdx mid = x.add(
        0.5 * (x.r[s] + x.r[e]),
        0.5 * (x.z[s] + x.z[e]),
        0.5 * (x.g[s] + x.g[e])
    );

    //old seg as left part, new seg as right part
    SegIdx newSeg = addSegment(mid, seg[idx].end, 0.5 * seg[idx].ds, seg[idx].interfacial);
    seg[idx].end = mid;
    seg[idx].ds *= 0.5;

    return newSeg;
}

double Graph::lengthSq(SegIdx idx){
    NodeIdx s = seg[idx].start;
    NodeIdx e = seg[idx].end;

    return sq(x.r[e] - x.r[s]) + sq(x.z[e] - x.z[s]);
}


void Graph::toFile(int step) const {
    std::string filename = "results/seg" + std::to_string(step);
    
    logg.formPrint(1, "writing to %s {", filename.c_str()).print("\n");

    FILE* f = std::fopen(filename.c_str(), "w");
    if (!f)
        throw std::runtime_error("Could not open file: " + filename);

    int colWidth = 25;

    std::fprintf(f, "%% %-*s %-*s %-*s %-*s %-*s\n",
       colWidth-3, "start",
       colWidth,   "end",
       colWidth,   "ds",
       colWidth,   "layers",
       colWidth,   "interfacial");

    for (const auto& s : seg) {
        std::fprintf(f, "%-*d %-*d %-*.17e %-*d %-*d\n", 
            colWidth, s.start, 
            colWidth, s.end, 
            colWidth, s.ds, 
            colWidth, s.layers,
            1       , s.interfacial ? 1 : 0);
    }

    std::fclose(f);

    logg.formPrint("%d segments written", seg.size()).print("\n");
    logg.formPrint(-1, "}\n");
}



void Graph::fromFile(int step) {
    std::string filename = "results/seg" + std::to_string(step);

    logg.formPrint(1, "reading from %s {", filename.c_str()).print("\n");

    FILE* f = std::fopen(filename.c_str(), "r");
    if (!f)
        throw std::runtime_error("Could not open file: " + filename);

    seg.clear();

    int ch;
    while ((ch = std::fgetc(f)) != EOF) {
        if (ch == '%') {
            int dummy = std::fscanf(f, "%*[^\n]\n"); // Skip comment line
            (void) dummy;
        } else if (ch != ' ' && ch != '\t' && ch != '\n' && ch != '\r') {
            std::ungetc(ch, f); // Found data line, push back character
            break;
        }
    }

    int start_val, end_val, inter_val, layers_val;
    double ds_val;
    while (std::fscanf(f, "%d %d %lf %d %d", &start_val, &end_val, &ds_val, &layers_val, &inter_val) == 5) {
        seg.emplace_back(start_val, end_val, ds_val, inter_val != 0, layers_val);
    }

    std::fclose(f);

    logg.formPrint("%d segments read", seg.size()).print("\n");
    computeWeights();
    logg.formPrint(-1, "}\n");
}


void Graph::determineConnectivity(){
    //reset cause some of the nodes were removed, some added, some segments become non-interfacial
    logg.formPrint("determining graph connectivity... ");

    for (auto& neighbors: C) {
        neighbors.clear();
    }
    C.resize(x.size());

    for(SegIdx i = 0; i < SegIdx(seg.size()); i++){
        const Segment& s = seg[i];

        if(s.start == s.end) 
            continue;

        C[s.start].push_back(Connection{s.end, i});
        
        if(s.interfacial)
            std::swap(C[s.start][0], C[s.start].back()); //prioritize interfacial path
        else
            C[s.end].push_back(Connection{s.start, i}); //cause no meaningful direction
    }
    logg.print("done.\n");
}
