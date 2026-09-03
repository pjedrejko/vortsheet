#include "Velocity.hpp"
#include "Params.hpp"
#include "elliptic.hpp"
#include <cmath>
#include "Tree.hpp"
#include "Utilities.hpp"

namespace point2point{
    double sq(double x){return x * x;}

    Prereq::Prereq(Array2 x0, Array2 x, double delta){
        //computing some prerequisites used later
        const auto [r0, z0] = x0;
        const auto [r,  z ] = x;

        R1sq = sq(r0 - r) + sq(z0 - z) + sq(delta);
        R2sq = sq(r0 + r) + sq(z0 - z) + sq(delta);
        R2cb = R2sq * std::sqrt(R2sq);

        ksq = 1 - (R1sq / R2sq);

        E = elliptic::E(ksq);
        K = elliptic::K(ksq);

        I1 = 4 * E / (1 - ksq);
        if(ksq > 1e-6)
            I2 = 8 * (K - E) / ksq;
        else // [0/0] -> Taylor expansion
            I2 = 2.0 * M_PI * (1.0 + (3.0 / 8.0) * ksq + (15.0 / 64.0) * ksq * ksq + (175.0 / 1024.0) * ksq * ksq * ksq);

    }

    Array2 biotSavart(const Prereq& p, Array2 x0, Array2 x, double Gr){
        const auto [r, z] = x;
        const auto [r0, z0] = x0;

        double Q = (0.25 / M_PI) * Gr / p.R2cb;
        double u0 = Q * (z0 - z) * (p.I1 - p.I2);
        double v0 = Q * (p.I1 * (r - r0) + p.I2 * r0);

        return {u0, v0};
    }

    Array2x2 uGradEval(const Prereq& p, Array2 x0C, Array2 xP, double SP){
        const auto [r0C, z0C] = x0C;
        const auto [rP, zP] = xP;

        double I3 = 4.0 * ( (2 * p.E - p.K) / (1 - p.ksq) + 2 * p.E / sq(1 - p.ksq) );
        double I4 = 1.0 * (  p.I2 + 8 * p.E / (1 - p.ksq) );
        double I5;
        if(p.ksq > 1e-6)
            I5 = 2.0 * ( (2.0 + p.ksq) * p.I2 - 8.0 * p.E ) / p.ksq;
        else // [0/0] -> Taylor expansion
            I5 = M_PI * (9.0 + p.ksq * ((15.0 / 4.0) + p.ksq * ((315.0 / 128.0) + p.ksq * (945.0 / 512.0))));

        Array2x2 A;
        A[0][0] =  (r0C - rP) * (z0C - zP) * (I4 - I3) + rP * (z0C - zP) * (I5 - I4);
        A[0][1] =  (z0C - zP) * (z0C - zP) * (I4 - I3);
        A[1][0] = -(r0C - rP) * (r0C - rP) * (I4 - I3) - rP *  r0C * I5;
        A[1][1] = -(r0C - rP) * (z0C - zP) * (I4 - I3) - rP * (z0C - zP) * I4;

        double Q = (0.25 / M_PI) * SP / p.R2cb;
        Array2x2 uGrad = {Q, Q, Q, Q};
        uGrad[0][0] *= ( A[0][0] / p.R2sq                 );
        uGrad[0][1] *= ( A[0][1] / p.R2sq + (p.I1 - p.I2) );
        uGrad[1][0] *= ( A[1][0] / p.R2sq - (p.I1 - p.I2) );
        uGrad[1][1] *= ( A[1][1] / p.R2sq                 );

        return uGrad;
    }

    Array2 biotSavartApprox(const Array2& uC, const Array2x2& uGrad, const Array2& Dx){
        return uC + matmul(uGrad, Dx);
    }
}


void Tree::Velocity::computeMoments(CellIdx idx){

    Cell& C = T.cells[idx];
    for(NodeIdx i = C.i0; i <= C.i1; i++){
        const NodeIdx j = T.node_map[i];
        const double dGamma = T.X->g[j] * T.G.W[j];

        if(dGamma > 0.0){
            C.SP    += dGamma * T.X->r[j];
            C.xP[0] += dGamma * T.X->r[j] * T.X->r[j];
            C.xP[1] += dGamma * T.X->r[j] * T.X->z[j];
        }
        else{
            C.SM    += dGamma * T.X->r[j];
            C.xM[0] += dGamma * T.X->r[j] * T.X->r[j];
            C.xM[1] += dGamma * T.X->r[j] * T.X->z[j];
        }
    }
}


void Tree::Velocity::computeAllMoments(CellIdx idx){

    Cell& C = T.cells[idx];

    if(!C.isLeaf()){
        for(CellIdx j: C.children){
            if(j < 0) continue;

            computeAllMoments(j);
            
            const Cell& Ch = T.cells[j];
            C.SP += Ch.SP;
            C.SM += Ch.SM;
            C.xP += Ch.xP;
            C.xM += Ch.xM;
        }
    }
    else
        computeMoments(idx);
}

void Tree::Velocity::backwardPass(){

    Timer::start("computing moments (backward pass)... ");

    computeAllMoments(0);

    //normalize moments to get barycenters
    for(Cell& C: T.cells){
        if(std::abs(C.SP) > 1e-14)
            C.xP /= C.SP;
        else{ 
            C.SP = 0.0;
            C.xP = {0.0, 0.0};
        }

        if(std::abs(C.SM) > 1e-14)
            C.xM /= C.SM;
        else{
            C.SM = 0.0;
            C.xM = {0.0, 0.0};
        }
    }

    Timer::end();
}


double Tree::Velocity::estimateEpsilon(CellIdx probe, CellIdx source){

    Array2 distP = T.cells[probe].xC - T.cells[source].xP;
    Array2 distM = T.cells[probe].xC - T.cells[source].xM;

    double dist = std::sqrt( std::min(dot(distP, distP), dot(distM, distM)) );
    
    if(dist < 1e-12) 
        return INFINITY;
    else
        return ( T.cells[probe].diag() + T.cells[source].diag() ) / dist;
}

void Tree::Velocity::biotSavartApprox(CellIdx probe, CellIdx source){
    namespace p2p = point2point;
    const Cell& CS = T.cells[source];
    Cell&       CP = T.cells[probe];

    //computing prerequisites (to avoid repetitions)
    p2p::Prereq pP    = p2p::Prereq(CP.xC, CS.xP, Params::delta); //positive vort. part
    p2p::Prereq pM    = p2p::Prereq(CP.xC, CS.xM, Params::delta); //negatvie vort. part

    Array2 uC{0, 0};
    uC += p2p::biotSavart(pP, CP.xC, CS.xP, CS.SP);
    uC += p2p::biotSavart(pM, CP.xC, CS.xM, CS.SM);

    Array2x2 uGradC{0, 0, 0, 0};
    uGradC += p2p::uGradEval(pP, CP.xC, CS.xP, CS.SP);
    uGradC += p2p::uGradEval(pM, CP.xC, CS.xM, CS.SM);

    for(NodeIdx i = CP.i0; i <= CP.i1; i++){
        NodeIdx j = T.node_map[i];
        Array2 x0 =  {T.X->r[j], T.X->z[j]};
        Array2 Dx = x0 - CP.xC;
        const auto [ur, uz] = p2p::biotSavartApprox(uC, uGradC, Dx);
        T.U->r[j] += ur;
        T.U->z[j] += uz;
    }
}


void Tree::Velocity::biotSavartExact(CellIdx probe, CellIdx source){
    namespace p2p = point2point;

    const Cell& CS = T.cells[source];
    Cell&       CP = T.cells[probe];

    for(NodeIdx ip = CP.i0; ip <= CP.i1; ip++){
    for(NodeIdx is = CS.i0; is <= CS.i1; is++){
        NodeIdx jp = T.node_map[ip];
        NodeIdx js = T.node_map[is];

        Array2 x0 =  {T.X->r[jp], T.X->z[jp]};
        Array2 x  =  {T.X->r[js], T.X->z[js]};
        double Gr = T.X->g[js] * T.G.W[js] * x[0];

        Array2 u = p2p::biotSavart(p2p::Prereq(x0, x, Params::delta), x0, x, Gr);
        T.U->r[jp] += u[0];
        T.U->z[jp] += u[1];
    }
    }
}

void Tree::Velocity::induce(CellIdx probe, CellIdx source){

    double eps = estimateEpsilon(probe, source);

    if( eps < Params::epsilon ){
        biotSavartApprox(probe, source);
    }
    else{
        if(T.cells[source].isLeaf()){
            biotSavartExact(probe, source);
        }
        else{
            for(CellIdx ch: T.cells[source].children)
                if(ch > -1) induce(probe, ch);
        }
    }
}



