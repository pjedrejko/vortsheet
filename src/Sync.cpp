#include"Sheet.hpp"
#include"Utilities.hpp"
#include<mpi.h>
#include<omp.h>

void Sheet::Sync::info() {
    using namespace Parallel;
    MPI_Comm_rank(MPI_COMM_WORLD, &iRank);
    MPI_Comm_size(MPI_COMM_WORLD, &nRanks);
    nThreads = omp_get_max_threads();

    logg.print("\n%s parallel setup %s\n", std::string(32, '=').c_str(), std::string(32, '=').c_str());
    logg.print("  MPI ranks:   %d\n", nRanks);
    logg.print("  OMP threads: %d\n", nThreads);
    logg.print("%s\n\n", std::string(70, '=').c_str());
}

void Sheet::Sync::determineLeafRange(){
    using namespace Parallel;
    MPI_Comm_rank(MPI_COMM_WORLD, &iRank);
    MPI_Comm_size(MPI_COMM_WORLD, &nRanks);

    LeafIdx nLeaves       = LeafIdx(sheet.T.leaves.size());
    LeafIdx leavesPerRank = nLeaves / nRanks;
    l0 = iRank * leavesPerRank;
    if(iRank == nRanks - 1)
        l1 = nLeaves - 1;
    else 
        l1 = l0 + leavesPerRank - 1;
}


void Sheet::Sync::brodcastNodes(Nodes& X){
    using namespace Parallel;
    Timer::start("brodcasting nodes... ");

    int nNodes = int(X.size());
    MPI_Bcast(&nNodes, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if(iRank != 0){
              X.resize(nNodes);
        sheet.W.resize(nNodes);
    }

    MPI_Bcast(    X.r.data(), nNodes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(    X.z.data(), nNodes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(    X.g.data(), nNodes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(sheet.W.data(), nNodes, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    Timer::end();
}

void Sheet::Sync::reduceVelocities(Nodes& U){
    Timer::start("reducing velocities... ");

    int nNodes = int(U.size());
    MPI_Allreduce(MPI_IN_PLACE, U.r.data(), nNodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, U.z.data(), nNodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    Timer::end();
}