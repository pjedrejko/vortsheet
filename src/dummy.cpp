    bool Tree::trySegTreeMerge(SegIdx s, CellIdx c){
        if(G->seg[s].state != SurgeryState::Untouched) 
            return false; //was already merged, no more merging

        Array2 segCenter = {G->segCenters.r[s], G->segCenters.z[s]};

        double distSq = minDistSq(segCenter, cells[c]);
        constexpr double dlMax = Params::dl * Params::kE;

        if( distSq < sq(dlMax) ){
            if(cells[c].isLeaf()){
                trySegCellMerge(s, c);
            }
            else{
                for(CellIdx ch: cells[c].children)
                    if(ch > -1) trySegTreeMerge(s, ch);
            }
        }
        return (G->seg[s].state != SurgeryState::Untouched); 
    }


    void Tree::induceVelocity(CellIdx probe, CellIdx source){

        double eps = estimateEpsilon(probe, source);

        if( eps < Params::epsilon ){
            biotSavartApprox(probe, source);
        }
        else{
            if(cells[source].isLeaf()){
                biotSavartExact(probe, source);
            }
            else{
                for(CellIdx ch: cells[source].children)
                    if(ch > -1) induceVelocity(probe, ch);
            }
        }
    }

    template<typename AreCloseFunc, typename ShortRangeFunc, typename LongRangeFunc, typename T>
    void traverse(T probe, CellIdx source){
        double eps = estimateEpsilon(probe, source);

        if( !AreCloseFunc(probe, source) ){
            LongRangeFunc(probe, source);
        }
        else{
            if(cells[source].isLeaf()){
                ShortRangeFunc(probe, source);
            }
            else{
                for(CellIdx ch: cells[source].children)
                    if(ch > -1) traverse<AreCloseFunc, ShortRangeFunc, LongRangeFunc, T>(probe, ch);
            }
        }
    }

    void induceVelocity<
    
    
    >
