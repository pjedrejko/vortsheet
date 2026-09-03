#include "Cell.hpp"
#include "Array.hpp"
#include <algorithm>

//double distSq1d(const Array2& axLim, const Array2& bxLim){
//    double d1 = axLim[0] - bxLim[1];
//    double d2 = bxLim[0] - axLim[1];
//
//    if( d1 <= 0.0 && d2 <= 0.0)
//        return 0.0; //overlap
//    else
//        return std::min(sq(d1), sq(d2));
//}
//
//double minDistSq(const Cell& c1, const Cell& c2){
//    return distSq1d(c1.rLim, c2.rLim) + distSq1d(c1.zLim, c2.zLim);
//}
//
//double minDistSq(const Array2& xLim, const Cell& c2){
//    return distSq1d({xLim[0], xLim[0]}, c2.rLim) 
//         + distSq1d({xLim[1], xLim[1]}, c2.zLim);
//}



double minDistSq(const Array2& pt, const Cell& c) {
    // Clamp point to cell bounding box to find closest point inside/on cell
    double closestR = std::clamp(pt[0], c.rLim[0], c.rLim[1]);
    double closestZ = std::clamp(pt[1], c.zLim[0], c.zLim[1]);

    return sq(pt[0] - closestR) + sq(pt[1] - closestZ);
}