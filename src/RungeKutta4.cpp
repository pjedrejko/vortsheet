#include "Sheet.hpp"
#include "Utilities.hpp"

void Sheet::RungeKutta4::resize(size_t count){
    k1.resize(count);
    k2.resize(count);
    k3.resize(count);
    k4.resize(count);
    temp_x.resize(count);
}

void Sheet::RungeKutta4::computeRHS(Nodes& x, Nodes& dxdt){
    dxdt.assign(x.size(), 0.0);
    sheet.computeVelocities (x, dxdt);
    
    if(!Parallel::isRoot()) return;
    sheet.produceCirculation(x, dxdt);
}

void Sheet::RungeKutta4::euler(Nodes& out, const Nodes& x, double dt, const Nodes& u) {
    if(!Parallel::isRoot()) return;

    const size_t N = x.size();
    for (size_t i = 0; i < N; i++) {
        out.r[i] = x.r[i] + dt * u.r[i];
        out.z[i] = x.z[i] + dt * u.z[i];
        out.g[i] = x.g[i] + dt * u.g[i];

        out.r[i] = std::max(out.r[i], 0.0); //domain guard
    }
    out.r[sheet.G.   topAxis] = 0.0; //periodicity guard
    out.r[sheet.G.bottomAxis] = 0.0;
}

void Sheet::RungeKutta4::stepEuler(double dt) {

    const size_t N = sheet.x.size();
    computeRHS(sheet.x, sheet.u);
    for(size_t i = 0; i < N; i++){
        sheet.x.r[i] += dt * sheet.u.r[i];
        sheet.x.r[i] = std::max(sheet.x.r[i], 0.0);
        sheet.x.z[i] += dt * sheet.u.z[i];
        sheet.x.g[i] += dt * sheet.u.g[i];
    }
    resize(N);
}

void Sheet::RungeKutta4::step(double dt) {
    logg.formPrint(1, "Runge-Kutta stepping {").print("\n");

    const size_t N = sheet.x.size();

    resize(N);

    logg.formPrint(1, "substep 1 {").print("\n");
    computeRHS(sheet.x, k1);
    sheet.u = k1; //store velocity for postprocessing etc.
    euler(temp_x, sheet.x, 0.5 * dt, k1);
    logg.formPrint(-1, "}\n");

    logg.formPrint(1, "substep 2 {").print("\n");
    computeRHS(temp_x, k2);
    euler(temp_x, sheet.x, 0.5 * dt, k2);
    logg.formPrint(-1, "}\n");

    logg.formPrint(1, "substep 3 {").print("\n");
    computeRHS(temp_x, k3);
    euler(temp_x, sheet.x, dt, k3);
    logg.formPrint(-1, "}\n");

    logg.formPrint(1, "substep 4 {").print("\n");
    computeRHS(temp_x, k4);

    const double dt6 = dt / 6.0;
    for (size_t i = 0; i < N; ++i) {
        sheet.x.r[i] += dt6 * (k1.r[i] + 2.0 * k2.r[i] + 2.0 * k3.r[i] + k4.r[i]);
        sheet.x.z[i] += dt6 * (k1.z[i] + 2.0 * k2.z[i] + 2.0 * k3.z[i] + k4.z[i]);
        sheet.x.g[i] += dt6 * (k1.g[i] + 2.0 * k2.g[i] + 2.0 * k3.g[i] + k4.g[i]);

        sheet.x.r[i] = std::max(sheet.x.r[i], 0.0); //domain guard
    }
    logg.formPrint(-1, "}\n");
    logg.formPrint(-1, "}\n");
}