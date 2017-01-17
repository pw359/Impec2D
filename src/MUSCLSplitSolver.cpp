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
/* MUSCLE solver
 *
 * See the book of Toro (Riemann Solvers and Numerical Methods for Fluid
 * Dynamics) for more information.
 */
//============================================================================


#include "MUSCLSplitSolver.h"
#include "HyperbolicSplitSolver.h"
#include <iostream>
#include "UpwindSplitSolver.h"
#include <cmath>
#include <algorithm>
#include <omp.h>

MUSCLSplitSolver::MUSCLSplitSolver(SimParams *sp) : HyperbolicSplitSolver(sp) {}
MUSCLSplitSolver::~MUSCLSplitSolver() {}

double MUSCLSplitSolver::ZetaL(double r, double w) {
    double beta_minus = 1;

    // beta_minus := beta_plus := 1 (centered scheme)
    return 2 * beta_minus * r / (1 - w + (1 + w) * r) ;
}

double MUSCLSplitSolver::ZetaR(double r, double w) {
    double beta_plus = 1;

    // beta_minus := beta_plus := 1 (centered scheme)
    return 2 * beta_plus  / (1 - w + (1 + w) * r)  ;
}

double MUSCLSplitSolver::SlopeLimiter(SLOPE_LIMITER sl, double du_m, double du_p, double w) {
    double r;
    if (du_m * du_p >0) {
        r = du_m/du_p;
        if (r > numeric_limits<double>::max())
            r = numeric_limits<double>::max();

        // Leerbee
        return std::min(    2. * r / (1. + r)    ,      ZetaR(r,w)           );
    }
    return 0;
}


double  MUSCLSplitSolver::GetSlope(double dU_prev, double dU_next) {
    double w = 0.;
    return SlopeLimiter(MINBEE, dU_prev, dU_next, w)   * 0.5 * ((1 + w) * dU_prev + (1-w) * dU_next);
}

void  MUSCLSplitSolver::Reconstruct(	double vL,double vR,double U_i, double U_ip1,double U_im1, double & UL_i, double & UR_i) {
    double dU_prev = U_i 		- U_im1;
    double dU_next = U_ip1 	- U_i;
    double d_i = GetSlope(dU_prev, dU_next);
    UR_i = U_i  + 0.5 * d_i;
    UL_i = U_i   - 0.5 * d_i;
}

/****
 * We want to calculate the flux F_im12j
 * To this end, we first need to interpolate all the states to the left and the right of each interface, respectively.
 * U_Lin(i, j,0) corresponds to the state interpolated to the interface i-12,j from the left
 * U_Lin(i, j, 1) corresponds to the state interpolated to the interface i-12,j from the right
******/

void  MUSCLSplitSolver::CalculateIntercellFlux_F(Array2d & F,const Array2d & c, const Array3d& v,	const double dt) {
    Array3d U_Lin;

    U_Lin.resize(sp->N1+3,(sp->N2+1) + 2 * sp->c_bounds,2);
    Idx3d idx3d;
    idx3d[0] = -1;
    idx3d[1] = -sp->c_bounds;
    idx3d[2] = 0;
    U_Lin.reindexSelf(idx3d);
    U_Lin = 0;

    // interpolate to right first
    double UL_i;
    double UR_i;


#ifdef OMPHYPER
    #pragma omp parallel for private(UL_i, UR_i)
#endif
    for (int i=-1; i <= sp->N1; i++) {
        for (int j=-sp->c_bounds; j < sp->N2 + sp->c_bounds; j++) {
            Reconstruct(v(i, j,0),v(i+1, j,0),  c(i, j), c(i+1, j),c(i-1, j), UL_i, UR_i);
            double df = 0.5 * dt/sp->dx * (Flux(v(i, j,0), UL_i) - Flux(v(i+1, j,0), UR_i));
            U_Lin(i+1, j,0)= 	UR_i + df;
            U_Lin(i, j, 1)= 		UL_i + df;

        }
    }
    // calculate the flux
    UpwindSplitSolver uss(sp);

#ifdef OMPHYPER
    #pragma omp parallel for
#endif
    for (int i=0; i <= sp->N1; i++) {
        for (int j=-sp->c_bounds; j < sp->N2 + sp->c_bounds; j++) {
            F(i, j) = uss.CalculateIntercellFlux(v(i, j,0), U_Lin(i, j,0), U_Lin(i, j, 1));
        }
    }
}


/****
 * We want to calculate the flux G_ij-12.
 * To this end, we first need to interpolate all the states to
 * the left and the right of each interface, respectively.
 * U_Lin(i, j,0) corresponds to the state interpolated to the interface i-12,j from the left
 * U_Lin(i, j, 1) corresponds to the state interpolated to the interface i-12,j from the right
******/

void  MUSCLSplitSolver::CalculateIntercellFlux_G( Array2d & G,const Array2d & c, const Array3d& v,	const double dt) {
    Array3d U_Lin;
    U_Lin.resize((sp->N1+1) + 2 * sp->c_bounds,sp->N2+3, 2);
    Idx3d idx3d;
    idx3d[0] =-sp->c_bounds;
    idx3d[1] = -1;
    idx3d[2] = 0;
    U_Lin.reindexSelf(idx3d);
    U_Lin = 0;

    // interpolate to right first
    double UL_i;
    double UR_i;

#ifdef OMPHYPER
    #pragma omp parallel for private (UL_i, UR_i)
#endif
    for (int i=-sp->c_bounds; i < sp->N1 + sp->c_bounds; i++) {
        for (int j=-1; j <= sp->N2; j++) {
            Reconstruct(v(i, j, 1),v(i, j+1, 1),  c(i, j), c(i, j+1),c(i, j-1), UL_i, UR_i);
            double df = 0.5 * dt/sp->dy * (Flux(v(i, j, 1), UL_i) - Flux(v(i, j+1, 1), UR_i));
            U_Lin(i, j+1,0)= 	UR_i + df;
            U_Lin(i, j, 1)= 		UL_i + df;
        }
    }
    // calculate the flux
    UpwindSplitSolver uss(sp);


#ifdef OMPHYPER
    #pragma omp parallel for
#endif
    for (int i=-sp->c_bounds; i < sp->N1 +sp->c_bounds; i++) {
        for (int j=0; j <= sp->N2; j++) {
            G(i, j) = uss.CalculateIntercellFlux(v(i, j, 1), U_Lin(i, j,0), U_Lin(i, j, 1));
        }
    }
}
