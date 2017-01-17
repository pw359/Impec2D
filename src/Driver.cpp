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
/* Driver.cpp
 *
 * The core functionality is encapsulated in Impec2DSolver.Solve().
 * Parameters are read from the provided file or otherwise from default.cfg.
 */
//============================================================================


#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <libconfig.h++>
#include "Impec2DSolver.h"
#include <string>
#include <omp.h>

using namespace std;
using namespace libconfig;

#include "SimParams.h"
#include "mpi.h"


int main(int argc, char ** argv)
{
    // Initialize MPI //
    MPI_Init(&argc, &argv);

    string configfile;
    if (argc == 2)
        configfile = argv[1];
    else
        configfile = "p_injection_diag.cfg";

    int myid;
    int num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    // CFD Solver 
    Impec2DSolver solver(configfile, myid, num_procs);

    if (solver.Run())
    {
        if (myid == 0)
            cout << "SUCCESSFUL" << endl;
    }
    else
        cerr << "ERROR: process " << myid << " did not terminate successfully." << endl;

    MPI_Finalize();
    return 0;
}


