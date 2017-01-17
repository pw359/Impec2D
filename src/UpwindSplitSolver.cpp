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


#include "UpwindSplitSolver.h"
#include "HyperbolicSplitSolver.h"
#include <iostream>
#include <omp.h>

UpwindSplitSolver::UpwindSplitSolver(SimParams *sp) : HyperbolicSplitSolver(sp) { }
UpwindSplitSolver::~UpwindSplitSolver() {}


void  UpwindSplitSolver::CalculateIntercellFlux_F(Array2d & F,const Array2d & c, const Array3d& v,	const double dt) {
#ifdef OMPHYPER
    #pragma omp parallel for
#endif
    for (int i=0; i <= sp->N1; i++) {
        for (int j=-sp->c_bounds; j < sp->N2 + sp->c_bounds; j++)
            F(i, j) = CalculateIntercellFlux(v(i, j,0), c(i-1, j), c(i, j));
    }
}

void  UpwindSplitSolver::CalculateIntercellFlux_G( Array2d & G,const Array2d & c, const Array3d& v,	const double dt) {
#ifdef OMPHYPER
    #pragma omp parallel for
#endif
    for (int i=-sp->c_bounds; i < sp->N1 +sp->c_bounds; i++) {
        for (int j=0; j <= sp->N2; j++)
            G(i, j) = CalculateIntercellFlux(v(i, j, 1),c(i, j-1), c(i, j));
    }
}

double  UpwindSplitSolver::CalculateIntercellFlux(double SPEED,double UL, double UR) {
    return (SPEED > 0) ? Flux(SPEED,UL) : Flux(SPEED, UR);
}


