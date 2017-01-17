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
/* HyperbolicSolver
 *
 * This class represents a layer of abstraction and encapsulates
 * functionalities required by all hyperbolic solvers.
 */
//============================================================================

#include "HyperbolicSolver.h"
#include "Util.h"
#include <libconfig.h++>
#include <iostream>
#include <cstdlib>
#include "TRIntegrator.h"
#include <omp.h>
using namespace std;
using namespace libconfig;

/*
 * Fill boundary cells according to the specified type
 */

void HyperbolicSolver::FillBoundaries(Array2d &c) {
    if (sp->is_test)
        FillGhostCellsPeriodicScalar(c,sp->c_bounds);
    else
        FillGhostCellsNeumannScalar(c, sp->c_bounds);
}

HyperbolicSolver::HyperbolicSolver(SimParams * sp) {
    this->sp = sp;
    mg = NULL;
    Initialize();
}

HyperbolicSolver::~HyperbolicSolver() {
    delete ti;
}

/*
 * Some schemes need the multigrid solver
 * NOTE: Initially this was not intended. However, I retained the structure for compatibility reasons.
 */

void HyperbolicSolver::AssignMGSolver(MGSolver *mg) {
    this->mg = mg;
}

void HyperbolicSolver::Initialize() {
    total_mass_outflux = 0;
    prev_first_sweep = Y_SWEEP;

    // use trapezoidal rule for time integration
    ti = new TRIntegrator(sp);
}

/*
 * Calcualte the maximum stable time step.
 */
double HyperbolicSolver::MaxStableTimestep(const Array2d & c,const  Array3d& v) {
    if (sp->is_test)
        return MaxStableAdvectionTimestep(c,v);
    else
        return MaxStableTimestepCombined(c,v);
}

double HyperbolicSolver::MaxStableAdvectionTimestep(const Array2d & c, const Array3d& v) {

    double dt;
    double dt_max = numeric_limits<double>::max();
    double v_centred, u_centred;

    for (int i=0; i<sp->N1; i++) {
        for (int j=0; j<sp->N2; j++) {
            u_centred = fabs((double)(0.5*(v(i, j,0) + v(i+1, j,0) )));
            v_centred = fabs((double)(0.5*(v(i, j, 1) + v(i, j+1, 1))));
            if (u_centred < numeric_limits<double>::min())
                u_centred = numeric_limits<double>::min();
            if (v_centred < numeric_limits<double>::min())
                v_centred = numeric_limits<double>::min();
            dt = min(sp->dx/u_centred,sp-> dy/v_centred);
            if (dt < dt_max)
                dt_max = dt;
        }
    }
    dt_max *= sp->cfl;
    return dt_max;
}

double HyperbolicSolver::MaxStableTimestepCombined(const Array2d & c, const Array3d& v) {

    // compute maximum stable timestep for advection
    blitz::Array<double,1> maxvel_threads;

    #pragma omp parallel
    {
        double v_centred, u_centred;
        int id = omp_get_thread_num();
        if (id == 0) {
            maxvel_threads.resize(omp_get_num_threads());
            maxvel_threads = 0;
        }
        #pragma omp barrier
        #pragma omp for
        for (int i=0; i<sp->N1; i++) {
            for (int j=0; j<sp->N2; j++) {
                u_centred = fabs((double)(0.5*(v(i, j,0) + v(i+1, j,0) )));
                v_centred = fabs((double)(0.5*(v(i, j, 1) + v(i, j+1, 1))));
                if (u_centred < numeric_limits<double>::min())
                    u_centred = numeric_limits<double>::min();
                if (v_centred < numeric_limits<double>::min())
                    v_centred = numeric_limits<double>::min();
                if (maxvel_threads(id) < fabs(u_centred) )
                    maxvel_threads(id) =  fabs(u_centred);
                if (maxvel_threads(id) < fabs(v_centred) )
                    maxvel_threads(id) =  fabs(v_centred);
            }
        }
    }
    double max_vel = max(maxvel_threads);
    double dt_advec_max = min(sp->dx/max_vel,sp-> dy/max_vel);
    dt_advec_max *= sp->cfl;

    double dt_diff_max = numeric_limits<double>::max();
    if (sp->diffusion) {

        // compute maximum stable timestep for diffusion
        blitz::Array<double,1> dt_diff_threads;
        #pragma omp parallel
        {
            double al = sp->alpha_l;
            double at = sp->alpha_t;
            double dx = sp->dx;
            int id = omp_get_thread_num();
            if (id == 0) {
                dt_diff_threads.resize(omp_get_num_threads());
                dt_diff_threads = 1.e10;
            }
            #pragma omp barrier
            #pragma omp for
            for (int i = 0; i < sp->N1; i++) {
                for (int j = 0; j <sp->N2; j++) {
                    double vij = 0.5 * (v(i,j,1) + v(i,j+1,1));
                    double uij = 0.5 * (v(i,j,0) + v(i+1,j,0));
                    double vij2 = vij * vij;
                    double uij2 = uij *uij;
                    double uij2pvij2  = vij2 + uij2;
                    double absu = sqrt(uij2pvij2);
                    double d11 = dx /absu * (al * uij2 + at* vij2);
                    double d22 = dx/ absu * (al * vij2 + at* uij2);
                    double dt = 1./2 *  dx * dx /  ( d11 + d22);
                    if (dt < dt_diff_threads(id))
                        dt_diff_threads(id) = dt;
                }
            }
        }
        dt_diff_max = min(dt_diff_threads) * sp->dt_diff_limit;
    }

    // if the parabolic scheme is Crank-Nicolson, then there is no timestep restriction
    if (sp->ptype == CRANKNICOLSON)
        return dt_advec_max;
    return min(dt_diff_max, dt_advec_max);

}


bool HyperbolicSolver::SourceSplit(Array2d & c, double const * const q, const double dt) {
    ti->Integrate(c,q,dt);
    this->total_mass_outflux += fabs(ti->GetTotalMassOutflux());
    return true;
}


double HyperbolicSolver::Flux (double vel, double conc) {
    return vel*conc/sp->porosity;
}


bool HyperbolicSolver::DiffusionSplit( Array2d &c, const Array3d &v, const double dt) {
    return this->DiffusionSplitCN(c,v,dt);
}

bool HyperbolicSolver::DiffusionSplitCN( Array2d &c, const Array3d &v, const double dt) {
    if (mg == NULL) {
        cout << "HyperbolicSolver:: Error> A Multigrid-Solver is needed for Crank-Nicolson!" << endl;
        exit(1);
    }
    return mg->Solve_Parabolic_cc_to_c(c, sp->c_bounds, v, sp->v_bounds, NULL, sp->ls_bounds, dt);
}
