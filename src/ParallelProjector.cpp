//=============================================================================
//  Impec2D: Implicit Pressure Explicit Concentration in two spatial dimensions
//  
//  This software solves the incompressible equations for miscible displacement 
//  for two symmetric configurations using a finite-volume solver. The under-
//  lying equations exhibit a physical instability. Small numerical errors,
//  e.g. arising from the particular choice of grid, can therefore lead to 
//  completely different solutions. Please have a look at my dissertation 
//  (goe.pdf) for more information.
//
//  Copyright:  
//
//  The software Impec2D was developed as part of my MPhil course on Scientific 
//  Computing in the Cavendish Laboratory, University of Cambridge. 
//  Copyright 2012, 2016 Peter Wirnsberger (peter.wirnsberger@gmail.com).
//
//  This file is part of Impec2D.
//
//  Impec2D is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  Impec2D is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with Impec2D.  If not, see <http://www.gnu.org/licenses/>.
//============================================================================


//============================================================================
/* ParallelProjector
 *
 * This class can map a part of the parallel solution
 * (solution for parallel configuration) onto a diagonal grid.
 * for the two symmetric configurations.
 */
//============================================================================


#include "ParallelProjector.h"
#include <iostream>
#include <stdlib.h>
#include <fstream>

using namespace std;

ParallelProjector::ParallelProjector(	double x_lb,double y_lb,	double x_ub,double y_ub, int d_N1, int d_N2, Array2d const * const parallel, int skipbounds ) {
    angle = -M_PI / 4.;

    // domain sizes of the parallel grid
    this->x_lb = x_lb;
    this->x_ub = x_ub;
    this->y_lb = y_lb;
    this->y_ub = y_ub;
    this->p = parallel;

    // resolution for the diagonal grid
    this->d_N1 = d_N1;
    this->d_N2 = d_N2;

    // resolution of parallel grid
    bounds = skipbounds;
    const int startx = (*p).base()[0]+bounds;
    const int starty = (*p).base()[1]+bounds;
    const int endx = startx +(*p).extent()[0] - 2*bounds;
    const int endy = starty + (*p).extent()[1] - 2*bounds;

    // resolution for the diagonal grid
    p_N1 = endx-startx;
    p_N2 = endy-starty;

    // dx , dy for parallel grid
    p_dx = (x_ub - x_lb)/ p_N1;
    p_dy = (y_ub - y_lb)/ p_N2;

    // dx , dy for diagonal grid
    d_dx = (x_ub - x_lb)/ (sqrt(2.) *d_N1);
    d_dy = (y_ub - y_lb)/ (sqrt(2.) * d_N2);
}

ParallelProjector::~ParallelProjector() {}

/*
 * Bilinear interpolation of the parallel solution to a point (x,y) corresponding to a cell centre of the new diagonal grid.
 */

double ParallelProjector::GetFromParallel(double x, double y ) {

    // fix negative y-coordinate using the symmetry
    y = fabs(y);

    if (x > x_ub || x < x_lb ||y > y_ub || y < y_lb) {
        cerr << "ParallelProjector" << endl;
        cerr << "\tx_lb: " << x_lb << endl;
        cerr << "\tx_ub: " << x_ub << endl;
        cerr << "\ty_lb: " << y_lb << endl;
        cerr << "\ty_ub: " << y_ub << endl;
        cerr << "\tOUT OF RANGE: (" << x << " , " << y << ")" << endl;
    }

    // four nearest neighbours for bilinear interpolation
    int x00, y00;
    int x10, y10;
    int x01, y01;
    int x11, y11;

    x00 = floor(x / (x_ub - x_lb) * p_N1);
    y00 = floor(y / (y_ub - y_lb) * p_N2);

    x10 = x00 + 1;
    y10 = y00 ;

    x01 = x00 ;
    y01 = y00 +1 ;

    x11 = x00 + 1;
    y11 = y00 + 1;

    double x_00, y_00;
    double x_11, y_11;

    x_00 = centred(x_lb, p_dx, x00);
    x_11 = centred(x_lb, p_dx, x11);
    y_00 = centred(y_lb, p_dy, y00);
    y_11 = centred(y_lb, p_dy, y11);

    return (	1. / ( (x_11-x_00) * (y_11-y_00)) * (x_11 - x) * (y_11 - y)  * (*p)(x00,y00) +
                1./ ( (x_11-x_00) * (y_11-y_00)) * ( x - x_00) * (y_11 - y) *  (*p)(x10,y10)+
                1./ ( (x_11-x_00) * (y_11-y_00)) * ( x_11 - x) * (y - y_00) *  (*p)(x01,y01) +
                1./ ( (x_11-x_00) * (y_11-y_00)) * ( x- x_00) * (y - y_00)	*  (*p)(x11, y11));
}

double ParallelProjector::centred(double lb, double dx, int idx) {
    return lb + (idx + 0.5) * dx;
}


double ParallelProjector::GetFromDiagonal(double x, double y ) {
    double p_x, p_y;
    DiagonalToParallel(x,y, p_x, p_y);
    return GetFromParallel(p_x,p_y);
}


void ParallelProjector::DiagonalToParallel(double d_x, double d_y, double & p_x, double & p_y) {
    /* Rotate coordinates <<angle>> degrees => base of parallel grid is rotated <<-angle>> degrees w.r.t the diagonal one
    *
    *  rotation matrix:	cos x	- sin x
    *  					sin x		  cos x
    */
    p_x = cos(angle) * d_x - sin(angle) * d_y;
    p_y = sin(angle) * d_x + cos(angle) * d_y;
}



void ParallelProjector::ParallelToDiagonal(double& d_x, double& d_y, double  p_x, double p_y) {
    /* Rotate coordinates <<-angle>> degrees => base of diagonal grid is rotated <<angle>> degrees w.r.t the parallel one
    *
    *  rotation matrix:	cos x	- sin x
    *  					sin x		  cos x
    */
    d_x = cos(-angle) * p_x - sin(-angle) * p_y;
    d_y = sin(-angle) * p_x + cos(-angle) * p_y;
}


void ParallelProjector::FillDiagonal(Array2d & d, int skip_bounds) {
    const int startx = d.base()[0]+skip_bounds;
    const int starty = d.base()[1]+skip_bounds;
    const int endx = startx +d.extent()[0] - 2*skip_bounds;
    const int endy = starty + d.extent()[1] - 2*skip_bounds;

    if (!( (endx-startx) == d_N1 && (endy-starty) == d_N2)) {
        cerr << "(d_N1,d_N2): (" << d_N1 << " , " << d_N2 << ")" << endl;
        cerr << "ParallelProjector:: Cannot fill array: DIMENSIONS ARE WRONG! Abort..." << endl;
        exit(1);
    }
    else {
        cout << "Diagonal array: 	(d_N1,d_N2) = (" << d_N1 << " , " << d_N2 << ")" << endl;
    }

#ifdef OMPPROJECTOR
    #pragma omp parallel for
#endif
    for (int i = startx; i < endx; i++) {
        for (int j = starty; j < endy; j++) {
            double p_x, p_y;
            double x = centred(x_lb, d_dx, i);
            double y =  centred(y_lb, d_dy, j);
            DiagonalToParallel(x,y, p_x, p_y);
            d(i, j) = GetFromParallel( p_x, p_y);
        }
    }
}

