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
/* TRIntegrator.cpp
 *
 * This class integrates c' = q using the trapezoidal rule.
 */
//============================================================================

#include "TRIntegrator.h"
#include "TimeIntegrator.h"
#include <omp.h>


TRIntegrator::TRIntegrator(SimParams *sp) : TimeIntegrator(sp) {}
TRIntegrator::~TRIntegrator() {}

int TRIntegrator::Integrate(Array2d & c, double const * const q, const double dt) {
    total_mass_outflux = 0.;
    double tmo = 0;
#ifdef OMPHYPER
    #pragma omp parallel for reduction(+:tmo)
#endif
    for (int i=0; i < sp->N1; i++) {
        for (int j=0; j<sp->N2; j++) {
            double qij =  q[i*sp->N2 +j];
            double massdiff;
            if (qij < 0) {
                massdiff = c(i,j);
                c(i,j) = c(i,j) * (1 + 0.5 * dt * qij)/ (1 - 0.5* dt * qij);
                massdiff = c(i,j) - massdiff;
                tmo += massdiff;
            }
            else if (qij > 0)
                c(i,j) += dt * sp->inj_fluid_conc * qij;
        }
    }
    total_mass_outflux = tmo;
    return 0;
}
