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
/* HyperbolicSplitSolver
 *
 * This class represents an additional layer of abstraction for the
 * HyperbolicSolver and encapsulates functionalities shared by all
 * split-solvers.
*/
//============================================================================


#include "HyperbolicSplitSolver.h"
#include <omp.h>

HyperbolicSplitSolver::HyperbolicSplitSolver(SimParams *sp) : HyperbolicSolver(sp) {}
HyperbolicSplitSolver::~HyperbolicSplitSolver() {}


bool HyperbolicSplitSolver::AdvectSplit(Array2d & c, const Array3d& v, double const * const q, const double dt) {
    bool ret = true;

    // copy boundary values from interior of the domain
    FillBoundaries(c);

    // determine which sweep comes first, depending on the history and the settings
    switch (sp->dim_splitting) {
    case ALTERNATING:
        if (prev_first_sweep == Y_SWEEP) {
            ret = ret && Advect_X(c, v,  dt);
            ret = ret && Advect_Y(c, v,  dt);
            prev_first_sweep = X_SWEEP;
        }
        else {
            ret = ret &&Advect_Y(c, v,  dt);
            ret = ret &&Advect_X(c,v,dt);
            prev_first_sweep = Y_SWEEP;
        }
        break;
    case X_FIRST:
        ret = ret &&Advect_X(c, v, dt);
        ret = ret &&Advect_Y(c, v, dt);
        prev_first_sweep = X_SWEEP;
        break;
    case Y_FIRST:
        ret = ret &&Advect_Y(c, v,  dt);
        ret = ret &&Advect_X(c,v,  dt);
        prev_first_sweep = Y_SWEEP;
        break;
    case STRANG:
        ret = ret &&Advect_X(c, v,  0.5*dt);
        ret = ret &&Advect_Y(c,v,  dt);
        ret = ret &&Advect_X(c,v,  0.5*dt);
        prev_first_sweep = X_SWEEP;
        break;
    case SYMMETRIC_PARALLEL:
        // concentration
        Array2d ctemp =c.copy();
        ret = ret && Advect_X(c, v, dt);
        ret = ret && Advect_Y(c, v,  dt);
        ret = ret && Advect_Y(ctemp, v,  dt);
        ret = ret && Advect_X(ctemp, v, dt);

#ifdef OMPHYPER
        #pragma omp parallel for
#endif
        for (int i=-sp->c_bounds; i < sp->N1 +sp->c_bounds; i++) {
            for (int j=-sp->c_bounds; j < sp->N2 +sp->c_bounds; j++)
                c(i, j) = 0.5 * (c(i, j) +ctemp(i, j));
        }
        prev_first_sweep = X_SWEEP;
        break;
    }
    return ret;
}


/*
 * Carry out x-sweep.
 */

bool HyperbolicSplitSolver::Advect_X(Array2d & c, const Array3d& v, const double dt) {
    Array2d F;

    F.resize(sp->N1 +1,sp->N2 + sp->c_bounds * 2);
    Idx2d idx;
    idx[0] = 0;
    idx[1] = -sp->c_bounds;
    F.reindexSelf(idx);
    F = 0;

    CalculateIntercellFlux_F(F, c, v,dt);

#ifdef OMPHYPER
    #pragma omp parallel for
#endif
    for (int i=0; i < sp->N1; i++) {
        for (int j=-sp->c_bounds; j < sp->N2 + sp->c_bounds; j++)
            c(i, j) += dt / sp->dx * (F(i, j) -  F(i+1, j) );
    }
    return true;
}

/*
 * Carry out y-sweep.
 */

bool HyperbolicSplitSolver::Advect_Y(Array2d & c, const Array3d& v, const double dt) {
    Array2d G;
    G.resize(sp->N1 + sp->c_bounds * 2,sp->N2 + 1);
    Idx2d idx;
    idx[0] = -sp->c_bounds;
    idx[1] = 0;
    G.reindexSelf(idx);
    G = 0;

    CalculateIntercellFlux_G(G, c, v,dt);
    // perform y-sweep

#ifdef OMPHYPER
    #pragma omp parallel for
#endif
    for (int i=-sp->c_bounds; i < sp->N1+sp->c_bounds; i++) {
        for (int j=0; j < sp->N2; j++)
            c(i, j) += dt / sp->dy * (G(i, j) -  G(i, j+1) );
    }
    return true;
}


bool HyperbolicSplitSolver::Solve(Array2d & c, const Array3d& v, double const * const q, const double dt) {
    bool ret = true;

    // 2nd order splitting

#ifdef SPLIT2ND

    // operational splitting: update by 1/2 dt
    // i) apply source term with operator splitting
    ret = ret && SourceSplit(c,q, 0.5*dt);

    // ii) apply diffusion term with operator splitting
    FillBoundaries(c);
    if (sp->diffusion)
        ret = ret && DiffusionSplit(c,v, 0.5 *dt);

#endif

    /// operational splitting: advect by dt
    // iii)  advection step
    ret = ret && AdvectSplit(c,v,q,dt);


#ifdef SPLIT2ND
    double dt_split = 0.5 * dt;
#else
    double dt_split = dt;
#endif

    // iv) apply source term with operator splitting
    ret = ret && SourceSplit(c,q,dt_split);
    // v) apply diffusion term with operator splitting
    FillBoundaries(c);
    if (sp->diffusion)
        ret = ret && DiffusionSplit(c,v, dt_split);
    return ret;
}


/*
 * Call these functions if you want to check, whether you have any influx from the domain boundaries.
 */

bool HyperbolicSplitSolver::BoundaryFluxCheckF(const Array2d & F) {
    double sum_influx = 0.;
    for (int j=0; j < sp->N2; j++) {

        // left
        sum_influx += fabs(F(0,j));
        if (fabs(F(0,j)) > 1.e-4) {
            cerr << "F[0][" << j << "]:" << F(0,j) << endl;
        }

        // right
        sum_influx += fabs((double)(F(sp->N1, j)));
        if (fabs((double)(F(sp->N1, j))) > 1.e-4) {
            cerr << "F[" << sp->N1 << "][" << j << "]:" << F(sp->N1, j) << endl;
        }
    }
    if (sum_influx > 1.e-8)
        cerr << "F:: INFLUX FROM BOUNDARY:  " << sum_influx << endl;
    return true;
}


bool HyperbolicSplitSolver::BoundaryFluxCheckG(const Array2d & G) {
    double sum_influx = 0.;
    for (int i=0; i < sp->N1; i++) {

        // bottom
        sum_influx += fabs(G(i,0));
        if (fabs(G(i,0)) > 1.e-4) {
            cerr << "G[" << i << "][0]:" << G(i,0) << endl;
        }

        // top
        sum_influx += fabs((double)(G(i, sp->N2)));
        if (fabs((double)(G(i, sp->N2))) > 1.e-4) {
            cerr << "G[" << i << "][" << sp->N2 <<"]:" << G(i, sp->N2) << endl;
        }
    }
    if (sum_influx > 1.e-8)
        cerr << "G::  INFLUX FROM BOUNDARY"<< endl;
    return true;
}
