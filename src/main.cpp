#include "Sheet.hpp"
#include "Velocity.hpp"
#include <string>
#include<algorithm>
#include<mpi.h>
#include<Utilities.hpp>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    Sheet sheet;
    sheet.parallel.info();

    const double circulation = 0.2;
    sheet.initializeRing(circulation);

//    sheet.fromFile(89);

    for(int stp = 0; stp < 90; stp++){
        if(Parallel::isRoot()){
            std::printf("%s\n", std::string(80, '=').c_str());
            sheet.integrals.compute();
            sheet.integrals.info();
            sheet.toFile();
        }
        
        sheet.rk4.step(sheet.integrals.dt);
        
        if(Parallel::isRoot()){
            sheet.u.toFile(sheet.integrals.step, "u");
            sheet.G.refine();
            sheet.T.surgery.apply();
        }
    }

    logg.print("all done.\n");
    MPI_Finalize();
    return 0;
}
