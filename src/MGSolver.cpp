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


/*==============================================================================

 This class uses the HYPRE multigrid solver package (from LLNL), that allows
 for solving very large systems (>1000x1000) in parallel using MPI.
 You can download the documentation here:
 https://computation.llnl.gov/casc/hypre/software.html.

//===========================================================================*/

#include "MGSolver.h"
#include "_hypre_utilities.h"
#include "HYPRE_struct_ls.h"
#include <fstream>
#include <cmath>
#include "mpi.h"
using namespace std;


MGSolver::MGSolver(SimParams *sp) : EllipticSolver(sp) {
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    // specify resolution and grid spacing
    Np = sqrt(num_procs);
    if (Np * Np != num_procs) {
        cout
                << "MGSolver::ConfigError: Processor grid has to be a perfect square! "
                << endl;
        exit(1);
    }

    // check whether the resolution is compatible with the processor grid
    if (sp->N1%Np != 0 || sp->N2%Np!=0) {
        cout << "MGSolver:: InitializationError> The processor grid is not compatible with the resolution (N1,N2)." << endl;
    }

    // dimension of the patch owned by this processor
    nx 	= sp->N1/Np;
    ny 	= sp->N2/Np;

    // specifiy domain for each processor in the grid
    pj = myid / Np;
    pi = myid - pj * Np;

    // extents of the patch
    ilower[0] = pi * nx;
    ilower[1] = pj * ny;

    iupper[0] = ilower[0] + nx - 1;
    iupper[1] = ilower[1] + ny - 1;

    // solver specific settings
    solver_id 	= 0;
    n_pre 		= 1;
    n_post 		= 1;
    max_iterations = 10000;

    // define stencil
    nentries = 9;
    stencil_indices = new int[nentries];
    nu = sp->nu;

    // n x n set up the grid
    SetupGrid();

    // set up the stencil
    SetupStencil();

    // viscosity of the reservoir fluid
    mu_res = sp->mu_res;

    // pre calculate this for spead-up
    mob1over4 = pow(sp->mu_res / sp->mu_inj, 0.25);

    // specify the cells for which you want Dirichlet boundary conditions
    if (sp->stype == PARALLEL || sp->stype == DIAGONAL) {
        N_ell_dir = 1;
        ell_dir  = (int**)malloc(N_ell_dir * sizeof(int*));
        for (int i=0; i<N_ell_dir; i++) {
            ell_dir[i] = (int*) malloc(2*sizeof(int));
            ell_dir[i][0] = 1;
            ell_dir[i][1] = 1;
        }
        N_para_dir = 0;
    }
    else if (sp->stype == HORIZONTAL) {

        // left border
        N_ell_dir = sp->N2;

        // N_ell_dir = 1;
        ell_dir  = (int**)malloc(N_ell_dir * sizeof(int*));
        for (int i=0; i<N_ell_dir; i++) {
            ell_dir[i] = (int*) malloc(2*sizeof(int));
            ell_dir[i][0] = 1;
            ell_dir[i][1] = i;
        }

        // right border
        N_para_dir = sp->N2;
        para_dir = (int**)malloc(N_para_dir *sizeof(int*));
        for (int i=0; i<N_para_dir; i++) {
            para_dir[i] = (int*) malloc(2*sizeof(int));
            para_dir[i][0] = sp->N1-2;
            para_dir[i][1] = i;
        }
    }
    else 
        cerr << "ERROR: Unknown source type!" << endl;
}

bool MGSolver::SolveLinearMultigrid(Array2d & c_c, int c_bounds, Array2d &p_c, int p_bounds, double const * const q_v) {
    bool stat =Solve_Elliptic_cc_to_c(c_c,c_bounds,p_c,p_bounds, q_v);
    this->last_its = this->num_iterations;
    return stat;
}


void MGSolver::SetupGrid() {
    HYPRE_StructGridCreate(MPI_COMM_WORLD, 2, &grid);
    HYPRE_StructGridSetExtents(grid, ilower, iupper);
    HYPRE_StructGridAssemble(grid);
}

void MGSolver::SetupStencil() {

    // define the stencil
    HYPRE_StructStencilCreate(2, nentries, &stencil);

    // relative index for stencil
    int offsets[9][2] = { { 0, 0 }, { -1, 0 }, { -1, -1 }, { 0, -1 }, { 1, -1 },
        { 1, 0 }, { 1, 1 }, { 0, 1 }, { -1, 1 }
    };
    for (int entry = 0; entry < nentries; entry++)
        HYPRE_StructStencilSetElement(stencil, entry, offsets[entry]);
}

void MGSolver::InitializeComponents() {

    // create empty matrix for stencil
    HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, stencil, &A);

    // coefficients are now ready to be set
    HYPRE_StructMatrixInitialize(A);

    HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &b);
    HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &x);

    // Indicate that the vector coefficients are ready to be set //
    HYPRE_StructVectorInitialize(b);
    HYPRE_StructVectorInitialize(x);
}

void MGSolver::DestroyComponents() {
    HYPRE_StructMatrixDestroy(A);
    HYPRE_StructVectorDestroy(b);
    HYPRE_StructVectorDestroy(x);
}

MGSolver::~MGSolver() {
    delete[] stencil_indices;
    if (N_ell_dir > 0) {
        for (int i=0; i<N_ell_dir; i++)
            free(ell_dir[i]);
        free(ell_dir);
    }
    if (N_para_dir >0) {
        for (int i=0; i<N_para_dir; i++)
            free(para_dir[i]);
        free(para_dir);
    }
    // Free memory //
    HYPRE_StructGridDestroy(grid);
    HYPRE_StructStencilDestroy(stencil);
}


/*******************************************************************************
 This method solves an elliptic equation of the form

    div ( K(x,y) /mu(c)  grad (p)) = -q

 where K is a 2x2 Tensor and mu an arbitrary function
 depending on an external quantity c.

 NOTE: ".. _cc_to_c" means centre to centre
*******************************************************************************/


bool MGSolver::Solve_Elliptic_cc_to_c(const Array2d & c_c, int c_bounds, Array2d &p_c, int p_bounds, double const * const q_c) {
    if (c_bounds <=0) {
        cout << "MGSolver::Error> At least one layer of ghost cells is required for the concentration aray!" << endl;
        return false;
    }

    // matrix and vectors coefficients are now ready to be set
    InitializeComponents();

    // set all the matrix entries
    SetMatrixCoefficients_Elliptic(c_c);

    // initialize b according to q
    SetRHS_Elliptic(q_c);

    // solve systems with selected solver
    SolveSystem(sp->eps);

    // map solution to desired array structure
    MapHYPREVectorToVertex(p_c, p_bounds);

    // free memory for A,x,b
    DestroyComponents();

    return true;
}


/************************************************************
 *
 * Set up the matrix for zero gradient Neumann boundaries.
 * If f(c) := 1 we have the following stencils for the elements
 * interior                      corners: (eq 0,0)
 *
 *		NW	 	N		NE
 *
 *		W		C		E
 *
 *		SW		S		SE
 *
 * For some cells (dir_i, dir_j) we use Dirichlet
 * boundaries, so that the problem is well-posed.
 *
 ************************************************************/

void MGSolver::SetMatrixCoefficients_Elliptic(const Array2d & c) {
    int bc_ilower[2];
    int bc_iupper[2];
    int idx[2];
    double vals[9];

    // set stencil for the interior everywhere and correct afterwards
    for (int j = 0; j < nentries; j++)
        stencil_indices[j] = j;

    // set the matrix coefficients in the interior
    double W, NW, N, NE, E, SE, S, SW, C;

    for (int k = ilower[0] ; k <= iupper[0]; k++) {
        for (int l = ilower[1]; l <= iupper[1]; l++) {
            GetWeightedMatrixCoefficientsForCell(nu, k, l, W, NW, N, NE, E,
                                                 SE, S, SW, C, c);

            vals[0] = C;
            vals[1] = W;
            vals[2] = SW;
            vals[ 3] = S;
            vals[ 4] = SE;
            vals[ 5] = E;
            vals[ 6] = NE;
            vals[ 7] = N;
            vals[ 8] = NW;

            idx[0] = k;
            idx[1] = l;
            HYPRE_StructMatrixSetValues(A, idx, nentries, stencil_indices,
                                        vals);
        }
    }


    // Processors at y = lb
    if (pj == 0) {
        bc_ilower[0] = pi * nx;
        bc_ilower[1] = pj * ny;
        bc_iupper[0] = bc_ilower[0] + nx - 1;
        bc_iupper[1] = bc_ilower[1];

        for (int k = bc_ilower[0]; k <= bc_iupper[0]; k++) {
            for (int l = bc_ilower[1]; l <= bc_iupper[1]; l++) {
                GetWeightedMatrixCoefficientsForCell(nu, k, l, W, NW, N, NE, E,
                                                     SE, S, SW, C, c);

                vals[0] = C+S;
                vals[1] = W+SW;
                vals[2] = 0;
                vals[3] = 0;
                vals[4] = 0;
                vals[5] = E+SE;
                vals[6] = NE;
                vals[7] = N;
                vals[8] = NW;


                idx[0] = k;
                idx[1] = l;
                HYPRE_StructMatrixSetValues(A, idx, nentries, stencil_indices,
                                            vals);
            }
        }
    }

    // Processors at y = ub  //
    if (pj == Np - 1) {
        bc_ilower[0] = pi * nx;
        bc_ilower[1] = pj * ny + ny - 1;
        bc_iupper[0] = bc_ilower[0] + nx - 1;
        bc_iupper[1] = bc_ilower[1];

        for (int k = bc_ilower[0]; k <= bc_iupper[0]; k++) {
            for (int l = bc_ilower[1]; l <= bc_iupper[1]; l++) {
                GetWeightedMatrixCoefficientsForCell(nu, k, l, W, NW, N, NE, E,
                                                     SE, S, SW, C, c);

                vals[0] = C+N;
                vals[1] = W+NW;
                vals[2] = SW;
                vals[3] = S;
                vals[4] = SE;
                vals[5] = E+NE;
                vals[6] = 0;
                vals[7] = 0;
                vals[8] = 0;
                idx[0] = k;
                idx[1] = l;
                HYPRE_StructMatrixSetValues(A, idx, nentries, stencil_indices,
                                            vals);
            }
        }
    }

    // Processors at x = lb 
    if (pi == 0) {
        bc_ilower[0] = pi * nx;
        bc_ilower[1] = pj * ny;

        bc_iupper[0] = bc_ilower[0];
        bc_iupper[1] = bc_ilower[1] + ny - 1;

        for (int k = bc_ilower[0]; k <= bc_iupper[0]; k++) {
            for (int l = bc_ilower[1]; l <= bc_iupper[1]; l++) {
                GetWeightedMatrixCoefficientsForCell(nu, k, l, W, NW, N, NE, E,
                                                     SE, S, SW, C, c);

                vals[0] = C + W;
                vals[1] = 0;
                vals[2] = 0;
                vals[3] = S+SW;
                vals[4] = SE;
                vals[5] = E;
                vals[6] = NE;
                vals[7] = N+NW;
                vals[8] = 0;

                idx[0] = k;
                idx[1] = l;
                HYPRE_StructMatrixSetValues(A, idx, nentries, stencil_indices,
                                            vals);
            }
        }
    }

    // Processors at x = ub
    if (pi == Np - 1) {
        bc_ilower[0] = pi * nx + nx - 1;
        bc_ilower[1] = pj * ny;
        bc_iupper[0] = bc_ilower[0];
        bc_iupper[1] = bc_ilower[1] + ny - 1;

        for (int k = bc_ilower[0]; k <= bc_iupper[0]; k++) {
            for (int l = bc_ilower[1]; l <= bc_iupper[1]; l++) {
                GetWeightedMatrixCoefficientsForCell(nu, k, l, W, NW, N, NE, E,
                                                     SE, S, SW, C, c);
                vals[0] = C+E;
                vals[1] = W;
                vals[2] = SW;
                vals[3] = S+SE;
                vals[4] = 0;
                vals[5] = 0;
                vals[6] = 0;
                vals[7] = N+NE;
                vals[8] = NW;

                idx[0] = k;
                idx[1] = l;
                HYPRE_StructMatrixSetValues(A, idx, nentries, stencil_indices,
                                            vals);
            }
        }
    }

    // Connections for the corner point (0,0)
    if (pi == 0 && pj == 0) {
        int k = 0;
        int l = 0;

        GetWeightedMatrixCoefficientsForCell(nu, k, l, W, NW, N, NE, E, SE, S,
                                             SW, C, c);
        vals[0] = C+W+S+SW;
        vals[1] = 0;
        vals[2] = 0;
        vals[3] = 0;
        vals[4] = 0;
        vals[5] =  E + SE;
        vals[6] = NE;
        vals[7] = N + NW;
        vals[8] = 0;

        idx[0] = k;
        idx[1] = l;
        HYPRE_StructMatrixSetValues(A, idx, nentries, stencil_indices, vals);
    }

    // Connections for the corner point (0,N*n-1)
    if (pi == 0 && pj == Np - 1) {
        int k = 0;
        int l = pj * ny + ny - 1;

        GetWeightedMatrixCoefficientsForCell(nu, k, l, W, NW, N, NE, E, SE, S,
                                             SW, C, c);

        vals[0] = C + W + N + NW;
        vals[1] = 0;
        vals[2] = 0;
        vals[3] = S + SW;
        vals[4] = SE;
        vals[5] = E + NE;
        vals[6] = 0;
        vals[7] = 0;
        vals[8] = 0;

        idx[0] = k;
        idx[1] = l;
        HYPRE_StructMatrixSetValues(A, idx, nentries, stencil_indices, vals);
    }

    // Connections for the corner point (N*n-1,N*n-1)
    if (pi == Np - 1 && pj == Np - 1) {
        int k = pi * nx + nx - 1;
        int l = pj * ny + ny - 1;

        GetWeightedMatrixCoefficientsForCell(nu, k, l, W, NW, N, NE, E, SE, S,
                                             SW, C, c);


        vals[0] = C + N + NE + E;
        vals[1] = W + NW ;
        vals[2] = SW;
        vals[3] = S + SE;
        vals[4] = 0;
        vals[5] = 0;
        vals[6] = 0;
        vals[7] = 0;
        vals[8] = 0;

        idx[0] = k;
        idx[1] = l;
        HYPRE_StructMatrixSetValues(A, idx, nentries, stencil_indices, vals);
    }

    // Connections for the corner point (N*n-1,0)
    if (pi == Np - 1 && pj == 0) {
        int k = pi * nx + nx - 1;
        int l = pj * ny;

        GetWeightedMatrixCoefficientsForCell(nu, k, l, W, NW, N, NE, E, SE, S,
                                             SW, C, c);

        vals[0] = C + S + E + SE;
        vals[1] = W + SW;
        vals[2] = 0;
        vals[3] = 0;
        vals[4] = 0;
        vals[5] = 0;
        vals[6] = 0;
        vals[7] = N + NE;
        vals[8] = NW;

        idx[0] = k;
        idx[1] = l;
        HYPRE_StructMatrixSetValues(A, idx, nentries, stencil_indices, vals);

    }

    // Set Dirichlet boundaries for specific cells
    for (int m=0; m<N_ell_dir; m++) {
        int dir_i = ell_dir[m][0];
        int dir_j = ell_dir[m][1];
        if (ilower[0] <= dir_i && ilower[1] <= dir_j && iupper[0] >= dir_i
                && iupper[1] >= dir_j) {

            // if it is a test case, then
            // choose coefficient for cell (dir_i, dir_j) s.t.
            // the pressure matches the pressure for the exact solution
            double kappa;
            if (sp->is_test) {
                double qij;
                double x = sp->x_lb + (dir_i+0.5) * sp->dx;
                double y = sp->y_lb + (dir_j+0.5) * sp->dy;
                switch (sp->testno) {
                case 1:
                    qij = qTest1(x,y);
                    break;
                case 2:
                    qij = qTest2(x,y);
                    break;
                case 3:
                    qij = qTest3(x,y);
                    break;
                case 4:
                    qij = qTest4(x,y);
                    break;
                }
                double pij = pTest(x,y);
                kappa = -dA * qij/pij;
            }
            // otherwise just set it to an arbitrary value
            else
                kappa = 1.;

            vals[0] =  kappa/ dA;
            vals[1] = 0;
            vals[2] = 0;
            vals[3] = 0;
            vals[4] = 0;
            vals[5] = 0;
            vals[6] = 0;
            vals[7] = 0;
            vals[8] = 0;

            idx[0] = dir_i;
            idx[1] = dir_j;
            HYPRE_StructMatrixSetValues(A, idx, nentries, stencil_indices, vals);
        }
    }

    HYPRE_StructMatrixAssemble(A);
}


void MGSolver::SetRHS_Elliptic(double const * const q) {
    int nvalues = nx * ny;
    double *values = (double*) calloc(nvalues, sizeof(double));

    HYPRE_StructVectorSetBoxValues(x, ilower, iupper, values);
    int idx[2];
    for (int k = ilower[0] ; k <= iupper[0]; k++) {
        for (int l = ilower[1]; l <= iupper[1]; l++) {
            idx[0] = k;
            idx[1] = l;
            HYPRE_StructVectorSetValues(b, idx, -q[k*sp->N2 + l]);
        }
    }
    free(values);
    // this is a collective call finalizing the vector assembly.
    HYPRE_StructVectorAssemble(b);
    HYPRE_StructVectorAssemble(x);
}



void MGSolver::SolveSystem(double tol) {

    if (solver_id == 0) {
        HYPRE_StructPCGCreate(MPI_COMM_WORLD, &solver);
        HYPRE_StructPCGSetMaxIter(solver, this->max_iterations);
        HYPRE_StructPCGSetTol(solver, tol);
        HYPRE_StructPCGSetTwoNorm(solver, 1);
        HYPRE_StructPCGSetRelChange(solver, 0);
        HYPRE_StructPCGSetPrintLevel(solver, 0); 
        HYPRE_StructPCGSetLogging(solver, 0);

        // Use symmetric SMG as preconditioner //
        HYPRE_StructSMGCreate(MPI_COMM_WORLD, &precond);
        HYPRE_StructSMGSetMemoryUse(precond, 0);
        HYPRE_StructSMGSetMaxIter(precond, 1);
        HYPRE_StructSMGSetTol(precond, 0.0);
        HYPRE_StructSMGSetZeroGuess(precond);
        HYPRE_StructSMGSetNumPreRelax(precond, 1);
        HYPRE_StructSMGSetNumPostRelax(precond, 1);

        // Set the preconditioner and solve //
        HYPRE_StructPCGSetPrecond(solver, HYPRE_StructSMGSolve,
                                  HYPRE_StructSMGSetup, precond);
        HYPRE_StructPCGSetup(solver, A, b, x);
        HYPRE_StructPCGSolve(solver, A, b, x);

        // Get some info on the run //
        HYPRE_StructPCGGetNumIterations(solver, &num_iterations);
        HYPRE_StructPCGGetFinalRelativeResidualNorm(solver, &final_res_norm);

        // Clean up //
        HYPRE_StructPCGDestroy(solver);

        // CLEANUP Pre conditioner
        HYPRE_StructSMGDestroy(precond);
    }

    if (solver_id == 1) {
        HYPRE_StructSMGCreate(MPI_COMM_WORLD, &solver);
        HYPRE_StructSMGSetMemoryUse(solver, 0);
        HYPRE_StructSMGSetMaxIter(solver, this->max_iterations);
        HYPRE_StructSMGSetTol(solver, tol);
        HYPRE_StructSMGSetRelChange(solver, 0);
        HYPRE_StructSMGSetNumPreRelax(solver, n_pre);
        HYPRE_StructSMGSetNumPostRelax(solver, n_post);
 
       // Logging must be on to get iterations and residual norm info below //
        HYPRE_StructSMGSetLogging(solver, 1);

        // Setup and solve //
        HYPRE_StructSMGSetup(solver, A, b, x);
        HYPRE_StructSMGSolve(solver, A, b, x);

        // Get some info on the run //
        HYPRE_StructSMGGetNumIterations(solver, &num_iterations);
        HYPRE_StructSMGGetFinalRelativeResidualNorm(solver, &final_res_norm);

        // Clean up //
        HYPRE_StructSMGDestroy(solver);
    }

    if (solver_id == 2) {
        HYPRE_StructJacobiCreate(MPI_COMM_WORLD, &solver);

        HYPRE_StructJacobiSetTol(solver, tol);
        HYPRE_StructJacobiSetMaxIter(solver, this->max_iterations);
        HYPRE_StructJacobiSetNonZeroGuess(solver);
        HYPRE_StructJacobiSetup(solver, A, b, x);
        HYPRE_StructJacobiSolve(solver, A, b, x);
        HYPRE_StructJacobiGetNumIterations(solver, &num_iterations);
        HYPRE_StructJacobiGetFinalRelativeResidualNorm(solver, &final_res_norm);
        HYPRE_StructJacobiDestroy(solver);
    }
}


void MGSolver::MapHYPREVectorToVertex(Array2d &p, int p_bounds) {

    int size_i, size_tot;
    size_i = nx*ny;
    size_tot = size_i * num_procs;
    double * buffer = new double[size_i];
    double val[1];
    int idx[2];
    for (int k=ilower[0]; k<= iupper[0]; k++) {
        for (int l=ilower[1]; l<= iupper[1]; l++) {
            idx[0] = k;
            idx[1] = l;
            HYPRE_StructVectorGetValues(x, idx, val);
            buffer[(k-ilower[0])*ny+(l-ilower[1])] = val[0];
        }
    }

    int tilower[2];
    int tpj, tpi;

    double * collect_buffer = new double[size_tot];
    MPI_Allgather(buffer, size_i, MPI_DOUBLE, collect_buffer, size_i, MPI_DOUBLE, MPI_COMM_WORLD);
    for (int pid = 0; pid<num_procs; pid++) {
        // specifiy domain for each processor in the grid
        tpj = pid / Np;
        tpi = pid - tpj * Np;
        tilower[0] = tpi * nx;
        tilower[1] = tpj * ny;

        for (int k=0; k<nx; k++) {
            for (int l=0; l<ny; l++)
                p(tilower[0]+k,tilower[1] +l) = collect_buffer[pid*size_i + k*ny+l];
        }
    }
    delete[] collect_buffer;
    delete [] buffer;
}



void MGSolver::GetWeightedMatrixCoefficientsForCell(double nu, int i, int j,
        double&W, double &NW, double & N, double & NE, double & E, double & SE,
        double &S, double &SW, double &C, const Array2d &c) {

    double C_p, W_p, SW_p, S_p, SE_p, E_p, NE_p, N_p, NW_p;
    double C_d, W_d, SW_d, S_d, SE_d, E_d, NE_d, N_d, NW_d;

    // coefficients for parallel stencil
    C_p 			= 0;
    W_p 			= 0;
    SW_p 	  	    = 0;
    S_p 			= 0;
    SE_p 		    = 0;
    E_p 			= 0;
    NE_p 		    = 0;
    N_p 			= 0;
    NW_p 		    = 0;

    // coefficients for diagonal stencil
    C_d 			= 0;
    W_d 			= 0;
    SW_d 	   	    = 0;
    S_d 			= 0;
    SE_d 		    = 0;
    E_d 			= 0;
    NE_d 		    = 0;
    N_d 			= 0;
    NW_d 		    = 0;

    // get the coefficients for the parallel stencil (five-point)
    GetParallelCoefficientsForCell(i, j, W_p, S_p, C_p, N_p, E_p, c);

    // get the coefficients for the diagonal stencil (nine-point)
    GetDiagonalCoefficientsForCell(i, j, W_d, NW_d, N_d, NE_d, E_d, SE_d, S_d,
                                   SW_d, C_d, c);


    // weight the coefficients according to the weighting factor nu
    C 	= nu * C_p 		+ (1 - nu) * C_d;
    W 	= nu * W_p 		+ (1 - nu) * W_d;
    SW 	= nu * SW_p 	+ (1 - nu) * SW_d;
    S 	= nu * S_p 		+ (1 - nu) * S_d;
    SE 	= nu * SE_p		+ (1 - nu) * SE_d;
    E 	= nu * E_p 		+ (1 - nu) * E_d;
    NE 	= nu * NE_p 	+ (1 - nu) * NE_d;
    N 	= nu * N_p 		+ (1 - nu) * N_d;
    NW 	= nu * NW_p 	+ (1 - nu) * NW_d;
}


void MGSolver::GetParallelCoefficientsForCell(int i, int j, double&W, double& S,
        double&C, double& N, double& E, const Array2d &c) {

    double K11_im12j 	= KoverMu_11_im12j((double) c(i-1,j), (double) c(i, j), i,	j);
    double K11_ip12j 	= KoverMu_11_im12j((double) c(i  ,j), (double) c(i + 1, j),i + 1, j);
    double K22_ijm12 	= KoverMu_22_ijm12((double) c(i  ,j - 1), (double) c(i, j), i,	j);
    double K22_ijp12 	= KoverMu_22_ijm12((double) c(i, j), (double) c(i, j + 1), i,j + 1);

    W 	  = 1. / (dx*dx) 	* K11_im12j;
    E 		= 1. / (dx*dx) 	* K11_ip12j;
    S 		= 1. / (dy*dy) 	* K22_ijm12;
    N 		= 1. / (dy*dy) 	* K22_ijp12;
    C 		= -1./ (dy*dy)	* (K22_ijm12 + K22_ijp12) - 1. / (dx*dx) * (K11_im12j + K11_ip12j);
}



void MGSolver::GetDiagonalCoefficientsForCell(int i, int j, double&W,
        double &NW, double & N, double & NE, double & E, double & SE, double &S,
        double &SW, double &C, const Array2d &c) {
    double M11_ip12j	= 0;
    double M11_im12j	= 0;
    double M22_ijp12 	= 0;
    double M22_ijm12	= 0;
    double M12_ip12j	= 0;
    double M21_ijp12	= 0;
    double M12_im12j	= 0;
    double M21_ijm12	= 0;

    GetTransformedMobilityTensor(i,j, M11_ip12j,	M11_im12j, M22_ijp12, M22_ijm12,
                                 M12_ip12j, M21_ijp12, M12_im12j, M21_ijm12,  c);

    NE 	= M11_ip12j 	/ (dxi * dxi);
    SW 	= M11_im12j 	/ (dxi * dxi);
    NW 	= M22_ijp12 	/ (deta * deta);
    SE 	= M22_ijm12 	/ (deta * deta);

    C 		= 		-1. / (dxi * dxi) 		* (M11_ip12j + M11_im12j)
                    -1. / (deta * deta)	  * (M22_ijp12 + M22_ijm12);

    N 		=  1. / (deta * dxi) * (M12_ip12j + M21_ijp12);
    W 	    = -1. / (deta * dxi) * (M21_ijp12 + M12_im12j);
    E 		= -1. / (deta * dxi) * (M12_ip12j + M21_ijm12);
    S 		=  1. / (deta * dxi) * (M12_im12j + M21_ijm12);
}


void MGSolver::GetTransformedMobilityTensor(int i,int j, double& M11_ip12j,double&	M11_im12j,double& M22_ijp12,double& M22_ijm12,
        double& M12_ip12j,double& M21_ijp12,double& M12_im12j,double& M21_ijm12,const Array2d &c) {

    const double K11_ip12jp12 = KoverMu_11_ip12jp12(i, j, c(i, j), c(i + 1, j),
                                c(i, j + 1), c(i + 1, j + 1));
    const double K11_im12jm12 = KoverMu_11_ip12jp12(i - 1, j - 1,
                                c(i - 1, j - 1), c(i, j - 1), c(i - 1, j), c(i, j));
    const double K11_im12jp12 = KoverMu_11_ip12jp12(i - 1, j, c(i - 1, j),
                                c(i, j), c(i - 1, j + 1), c(i, j + 1));
    const double K11_ip12jm12 = KoverMu_11_ip12jp12(i, j - 1, c(i, j - 1),
                                c(i + 1, j - 1), c(i, j), c(i + 1, j));

    const double K22_ip12jp12 = KoverMu_22_ip12jp12(i, j, c(i, j), c(i + 1, j),
                                c(i, j + 1), c(i + 1, j + 1));
    const double K22_im12jm12 = KoverMu_22_ip12jp12(i - 1, j - 1,
                                c(i - 1, j - 1), c(i, j - 1), c(i - 1, j), c(i, j));
    const double K22_im12jp12 = KoverMu_22_ip12jp12(i - 1, j, c(i - 1, j),
                                c(i, j), c(i - 1, j + 1), c(i, j + 1));
    const double K22_ip12jm12 = KoverMu_22_ip12jp12(i, j - 1, c(i, j - 1),
                                c(i + 1, j - 1), c(i, j), c(i + 1, j));


    M11_ip12j 	= secamt2 * (  K11_ip12jp12 *  cosa2      +  K22_ip12jp12    *     sina2 );
    M11_im12j 	= secamt2 * (  K11_im12jm12 *  cosa2      +  K22_im12jm12    *     sina2 );

    M22_ijp12 	= secamt2 * (  K22_im12jp12 *  cost2      +  K11_im12jp12    *     sint2 );
    M22_ijm12 	= secamt2 * (  K22_ip12jm12 *  cost2      +  K11_ip12jm12    *     sint2 );

    M12_ip12j 	= secamt2 * (K22_ip12jp12   *  costsina   -  K11_ip12jp12    *     cosasint);
    M21_ijp12 	= secamt2 * (K22_im12jp12   *  costsina   -  K11_im12jp12    *     cosasint);

    M12_im12j 	= secamt2 * (K22_im12jm12   *  costsina   -  K11_im12jm12    *     cosasint);
    M21_ijm12 	= secamt2 * (K22_ip12jm12   *  costsina   -  K11_ip12jm12    *     cosasint);

}




void MGSolver::SetupDiffusionTensorF(double dx, double dy, const Array3d &v, int i, int j,	double &d11_im12j, double &d12_im12j) {

    double u_im12j		= v(i, j,0);
    double v_im12j		= 1./4 * (v(i-1, j+1, 1) + v(i, j+1, 1) + v(i-1, j, 1) + v(i, j, 1));
    double u2_im12j = u_im12j * u_im12j;
    double v2_im12j = v_im12j * v_im12j;

    double phi = sp->porosity;
    double maxdx = max(dx,dy);
    double alpha_l_weight;
    double alpha_t_weight;

    /*****************************************
     *			WEIGHTS FOR TRANSVERSE AND
     *			LONGITUDINAL DIFFUSION
     *****************************************/
    alpha_l_weight 	= sp->alpha_l;
    alpha_t_weight 	= sp->alpha_t;

    double alpha_l		  =	alpha_l_weight * maxdx;
    double alpha_t		  =	alpha_t_weight * maxdx;

    double fac_l_im12j 	=	phi * alpha_l/ (sqrt(u2_im12j + v2_im12j));
    double fac_t_im12j 	=	phi * alpha_t/ (sqrt(u2_im12j + v2_im12j));

    // longintudinal contribution
    double d11_l_im12j 	= fac_l_im12j	*  	u2_im12j;
    double d12_l_im12j	= fac_l_im12j	* 	(u_im12j * v_im12j);

    //	transverse contribution
    double d11_t_im12j 	= fac_t_im12j	*  v2_im12j;
    double d12_t_im12j	= fac_t_im12j	*	(- u_im12j * v_im12j);

    d11_im12j   = d11_l_im12j 	+ 	d11_t_im12j;
    d12_im12j	  = d12_l_im12j 	+	  d12_t_im12j;
}




void MGSolver::SetupDiffusionTensorG(double dx, double dy, const Array3d &v, int i, int j,	double &d22_ijm12, double &d21_ijm12) {

    double u_ijm12 		=	1./4 * (v(i+1, j-1,0)   + v(i+1, j,0) + v(i, j-1,0) + v(i, j,0));
    double v_ijm12		= v(i, j, 1);
    double u2_ijm12 = u_ijm12 * u_ijm12;
    double v2_ijm12 = v_ijm12 * v_ijm12;

    double phi = sp->porosity;
    double maxdx = max(dx,dy);
    double alpha_l_weight;
    double alpha_t_weight;

    /*****************************************
     *			WEIGHTS FOR TRANSVERSE AND
     *			LONGITUDINAL DIFFUSION
     *****************************************/
    alpha_l_weight 	= sp->alpha_l;
    alpha_t_weight 	= sp->alpha_t;
    double alpha_l	  	=	alpha_l_weight * maxdx;
    double alpha_t			=	alpha_t_weight * maxdx;

    double fac_l_ijm12 	=	phi * alpha_l/ (sqrt(u2_ijm12 + v2_ijm12));
    double fac_t_ijm12 	=	phi * alpha_t/ (sqrt(u2_ijm12 + v2_ijm12));

    // longintudinal contribution
    double d21_l_ijm12		=	fac_l_ijm12 * 	(u_ijm12 * v_ijm12);
    double d22_l_ijm12		=	fac_l_ijm12	* 	v2_ijm12;
    //	transverse contribution
    double d21_t_ijm12	=	fac_t_ijm12 * 	(- u_ijm12 * v_ijm12);
    double d22_t_ijm12	=	fac_t_ijm12	* 	u2_ijm12;

    d21_ijm12	=	d21_l_ijm12	+	d21_t_ijm12;
    d22_ijm12	=	d22_l_ijm12	+	d22_t_ijm12;
}



void MGSolver::GetCoefficients_Parabolic(int i, int j, double&W, double &NW, double & N, double & NE, double & E,
        double & SE, double &S, double &SW, double &C,const Array3d &v, const double dt, Array2d const * const ls) {

    double d11_im12j, d12_im12j;
    SetupDiffusionTensorF(dx,dy, v,i,j,d11_im12j,d12_im12j);

    double d22_ijm12, d21_ijm12;
    SetupDiffusionTensorG(dx,dy, v,i,j,d22_ijm12, d21_ijm12);

    double d21_ijp12, d22_ijp12;
    SetupDiffusionTensorG(dx,dy, v,i,j+1,d22_ijp12, d21_ijp12);

    double d12_ip12j, d11_ip12j;
    SetupDiffusionTensorF(dx,dy, v,i+1,j,d11_ip12j,d12_ip12j);


    NW  =  dt/(8.*dx*dy)    * (d12_im12j + d21_ijp12);
    NE  = -dt/(8.*dx*dy)    * (d12_ip12j + d21_ijp12);
    N   = -dt/8. * (4./(dy*dy) * d22_ijp12 + 1./(dx*dy) * (-d12_im12j + d12_ip12j));
    SW  = -dt/(8.*dx*dy)       * (d12_im12j + d21_ijm12);
    SE  =  dt/(8.*dx*dy)    * (d12_ip12j + d21_ijm12);
    S   = -dt/8. * (4./(dy*dy) * d22_ijm12  +  1./(dx*dy) * (-d12_ip12j + d12_im12j));
    W   = -dt/8. * (4./(dx*dx) * d11_im12j  +  1./(dx*dy) * (-d21_ijp12 + d21_ijm12));
    C   = 1. +  dt/(2.*dx*dx) * (d11_ip12j + d11_im12j) + dt/(2.*dy*dy)* (d22_ijp12 + d22_ijm12);
    E   = -dt/8. * (  4./(dx*dx) * d11_ip12j + 1./(dx*dy) * (-d21_ijm12 + d21_ijp12)  );
}



bool MGSolver::Solve_Parabolic_cc_to_c(Array2d & c_c, int c_bounds, const Array3d &v, int v_bounds,
                                       Array2d const * const ls, int ls_bounds, const double dt) {

    // Allocate matrix and vectors
    InitializeComponents();

    // set all the matrix entries for the parabolic problem
    SetMatrixCoefficients_Parabolic(v,ls, dt);

    // initialize b according to q
    SetRHS_Parabolic(c_c, v,dt,ls);

    // solve systems with selected solver
    SolveSystem(sp->para_eps);

    MapHYPREVectorToVertex(c_c, c_bounds);

    // free memory for A,x,b
    DestroyComponents();

    // save number of iterations for output
    this->last_parabolic_its = this->num_iterations;

    return true;
}


/************************************************************
 *
 * Set up the matrix for zero gradient Neumann boundaries.
 * for the parabolic problem.
 * If f(c) := 1 we have the following stencils for the elements
 * 	interior                      corners: (eq 0,0)
 *
 *		NW	 	N	        NE
 *
 *		W	   	C		E
 *
 *		SW		S		SE
 *
 * For one cell in the grid (dir_i, dir_j) we use 0 dirichlet
 * boundaries, so that the problem is well-posed.
 *
 ************************************************************/

void MGSolver::SetMatrixCoefficients_Parabolic(const Array3d & v, Array2d const * const ls, const double dt) {
    int bc_ilower[2];
    int bc_iupper[2];
    int idx[2];
    double vals[9];

    // set stencil for the interior everywhere and correct afterwards
    for (int j = 0; j < nentries; j++)
        stencil_indices[j] = j;

    double W, NW, N, NE, E, SE, S, SW, C;

    for (int k = ilower[0] ; k <= iupper[0]; k++) {
        for (int l = ilower[1]; l <= iupper[1]; l++) {
            GetCoefficients_Parabolic(k, l, W, NW, N, NE, E, SE,S, SW, C, v,dt, ls);

            vals[0] = C;
            vals[1] = W;
            vals[2] = SW;
            vals[ 3] = S;
            vals[ 4] = SE;
            vals[ 5] = E;
            vals[ 6] = NE;
            vals[ 7] = N;
            vals[ 8] = NW;

            idx[0] = k;
            idx[1] = l;
            HYPRE_StructMatrixSetValues(A, idx, nentries, stencil_indices,
                                        vals);
        }
    }


    // Processors at y = lb
    if (pj == 0) {
        // The stencil at the boundary nodes is
        //		NW+SW		 N+S 		NE + SE
        //
        //			W			      C				E
        //
        //			0				  0			 	0

        bc_ilower[0] = pi * nx;
        bc_ilower[1] = pj * ny;
        bc_iupper[0] = bc_ilower[0] + nx - 1;
        bc_iupper[1] = bc_ilower[1];

        for (int k = bc_ilower[0]; k <= bc_iupper[0]; k++) {
            for (int l = bc_ilower[1]; l <= bc_iupper[1]; l++) {
                GetCoefficients_Parabolic(k, l, W, NW, N, NE, E, SE,S, SW, C, v,dt,ls);

                vals[0] = C+S;
                vals[1] = W+SW;
                vals[2] = 0;
                vals[3] = 0;
                vals[4] = 0;
                vals[5] = E+SE;
                vals[6] = NE;
                vals[7] = N;
                vals[8] = NW;

                idx[0] = k;
                idx[1] = l;
                HYPRE_StructMatrixSetValues(A, idx, nentries, stencil_indices,vals);
            }
        }
    }

    // Processors at y = ub  //
    if (pj == Np - 1) {
        // The stencil at the boundary nodes is
        //			0				  0			 	0
        //
        //			W			      C				E
        //
        //		NW+SW		 N+S 		NE + SE

        bc_ilower[0] = pi * nx;
        bc_ilower[1] = pj * ny + ny - 1;
        bc_iupper[0] = bc_ilower[0] + nx - 1;
        bc_iupper[1] = bc_ilower[1];

        for (int k = bc_ilower[0]; k <= bc_iupper[0]; k++) {
            for (int l = bc_ilower[1]; l <= bc_iupper[1]; l++) {
                GetCoefficients_Parabolic(k, l, W, NW, N, NE, E, SE,S, SW, C, v,dt,ls);

                vals[0] = C+N;
                vals[1] = W+NW;
                vals[2] = SW;
                vals[3] = S;
                vals[4] = SE;
                vals[5] = E+NE;
                vals[6] = 0;
                vals[7] = 0;
                vals[8] = 0;

                idx[0] = k;
                idx[1] = l;
                HYPRE_StructMatrixSetValues(A, idx, nentries, stencil_indices,
                                            vals);
            }
        }
    }

    // Processors at x = lb 
    if (pi == 0) {

        // The stencil at the boundary nodes is
        //			0				  N			 	NW+NE
        //
        //			0			      C				W+E
        //
        //		 	0 	 			  S 				SE+SW

        bc_ilower[0] = pi * nx;
        bc_ilower[1] = pj * ny;

        bc_iupper[0] = bc_ilower[0];
        bc_iupper[1] = bc_ilower[1] + ny - 1;

        for (int k = bc_ilower[0]; k <= bc_iupper[0]; k++) {
            for (int l = bc_ilower[1]; l <= bc_iupper[1]; l++) {
                GetCoefficients_Parabolic(k, l, W, NW, N, NE, E, SE,S, SW, C, v,dt,ls);

                vals[0] = C + W;
                vals[1] = 0;
                vals[2] = 0;
                vals[3] = S+SW;
                vals[4] = SE;
                vals[5] = E;
                vals[6] = NE;
                vals[7] = N+NW;
                vals[8] = 0;

                idx[0] = k;
                idx[1] = l;
                HYPRE_StructMatrixSetValues(A, idx, nentries, stencil_indices,
                                            vals);
            }
        }
    }

    // Processors at x = ub
    if (pi == Np - 1) {
        // The stencil at the boundary nodes is
        //			NW+NE	  N			 	0
        //
        //			W+E		      C				0
        //
        //		 	SE+SW		  S 				0

        bc_ilower[0] = pi * nx + nx - 1;
        bc_ilower[1] = pj * ny;
        bc_iupper[0] = bc_ilower[0];
        bc_iupper[1] = bc_ilower[1] + ny - 1;

        for (int k = bc_ilower[0]; k <= bc_iupper[0]; k++) {
            for (int l = bc_ilower[1]; l <= bc_iupper[1]; l++) {
                GetCoefficients_Parabolic(k, l, W, NW, N, NE, E, SE,S, SW, C, v,dt,ls);

                vals[0] = C+E;
                vals[1] = W;
                vals[2] = SW;
                vals[3] = S+SE;
                vals[4] = 0;
                vals[5] = 0;
                vals[6] = 0;
                vals[7] = N+NE;
                vals[8] = NW;
                idx[0] = k;
                idx[1] = l;
                HYPRE_StructMatrixSetValues(A, idx, nentries, stencil_indices,
                                            vals);
            }
        }
    }

    // Connections for the corner point (0,0)
    if (pi == 0 && pj == 0) {

        // The stencil at the boundary nodes is
        //		0			S + N	 	NE + SE + SW + NW
        //
        //		0		       C		    W + E
        //
        //		0		  	   0			0

        int k = 0;
        int l = 0;

        GetCoefficients_Parabolic(k, l, W, NW, N, NE, E, SE,S, SW, C, v,dt,ls);

        vals[0] = C+W+S+SW;
        vals[1] = 0;
        vals[2] = 0;
        vals[3] = 0;
        vals[4] = 0;
        vals[5] =  E + SE;
        vals[6] = NE;
        vals[7] = N + NW;
        vals[8] = 0;

        idx[0] = k;
        idx[1] = l;
        HYPRE_StructMatrixSetValues(A, idx, nentries, stencil_indices, vals);
    }

    // Connections for the corner point (0,N*n-1)
    if (pi == 0 && pj == Np - 1) {
        // The stencil at the boundary nodes is
        //		0			 	0		 		0
        //
        //		0		        C			W + E
        //
        //		0		  	   N + S	NW + SW + NE + SE

        int k = 0;
        int l = pj * ny + ny - 1;

        GetCoefficients_Parabolic(k, l, W, NW, N, NE, E, SE,S, SW, C, v,dt,ls);

        vals[0] = C + W + N + NW;
        vals[1] = 0;
        vals[2] = 0;
        vals[3] = S + SW;
        vals[4] = SE;
        vals[5] = E + NE;
        vals[6] = 0;
        vals[7] = 0;
        vals[8] = 0;

        idx[0] = k;
        idx[1] = l;
        HYPRE_StructMatrixSetValues(A, idx, nentries, stencil_indices, vals);
    }

    // Connections for the corner point (N*n-1,N*n-1)
    if (pi == Np - 1 && pj == Np - 1) {
        // The stencil at the boundary nodes is
        //		0			 							0		 		0
        //
        //		W+E	       							 C			 	0
        //
        //		NW + SW + NE + SE	  N + S			0

        int k = pi * nx + nx - 1;
        int l = pj * ny + ny - 1;

        GetCoefficients_Parabolic(k, l, W, NW, N, NE, E, SE,S, SW, C, v,dt,ls);


        vals[0] = C + N + NE + E;
        vals[1] = W + NW ;
        vals[2] = SW;
        vals[3] = S + SE;
        vals[4] = 0;
        vals[5] = 0;
        vals[6] = 0;
        vals[7] = 0;
        vals[8] = 0;

        idx[0] = k;
        idx[1] = l;
        HYPRE_StructMatrixSetValues(A, idx, nentries, stencil_indices, vals);
    }

    // Connections for the corner point (N*n-1,0)
    if (pi == Np - 1 && pj == 0) {
        // The stencil at the boundary nodes is
        //		NW + SW + NE + SE		N + S		 		0
        //
        //		W+E	       							 C					 	0
        //
        //		  	0									 0						0

        int k = pi * nx + nx - 1;
        int l = pj * ny;

        GetCoefficients_Parabolic(k, l, W, NW, N, NE, E, SE,S, SW, C, v,dt,ls);
        vals[0] = C + S + E + SE;
        vals[1] = W + SW;
        vals[2] = 0;
        vals[3] = 0;
        vals[4] = 0;
        vals[5] = 0;
        vals[6] = 0;
        vals[7] = N + NE;
        vals[8] = NW;
        idx[0] = k;
        idx[1] = l;
        HYPRE_StructMatrixSetValues(A, idx, nentries, stencil_indices, vals);
    }





    for (int m=0; m<N_para_dir; m++) {
        int dir_i = para_dir[m][0];
        int dir_j = para_dir[m][1];

        // Set Dirichlet boundaries for  cell (dir_i,dir_j)
        if (ilower[0] <= dir_i && ilower[1] <= dir_j && iupper[0] >= dir_i
                && iupper[1] >= dir_j) {

            vals[0] = 1./ dA;
            vals[1] = 0;
            vals[2] = 0;
            vals[3] = 0;
            vals[4] = 0;
            vals[5] = 0;
            vals[6] = 0;
            vals[7] = 0;
            vals[8] = 0;
            idx[0] = dir_i;
            idx[1] = dir_j;

            HYPRE_StructMatrixSetValues(A, idx, nentries, stencil_indices, vals);
        }
    }


    HYPRE_StructMatrixAssemble(A);
}


void MGSolver::SetRHS_Parabolic(const Array2d &c, const Array3d &v, const double dt, Array2d const * const ls) {
    int nvalues = nx * ny;
    double * rhs = new double[nvalues];

    for (int i=ilower[0]; i  <= iupper[0];  i++) {
        for (int j=ilower[1]; j  <= iupper[1]; j++) {

            // Calculate the elements of the diffusion tensor at the vertices
            double d11_im12j, d12_im12j;
            SetupDiffusionTensorF(dx,dy, v,i,j,d11_im12j,d12_im12j);

            double d22_ijm12, d21_ijm12;
            SetupDiffusionTensorG(dx,dy, v,i,j,d22_ijm12, d21_ijm12);

            double d21_ijp12, d22_ijp12;
            SetupDiffusionTensorG(dx,dy, v,i,j+1,d22_ijp12, d21_ijp12);

            double d12_ip12j, d11_ip12j;
            SetupDiffusionTensorF(dx,dy, v,i+1,j,d11_ip12j,d12_ip12j);

            // set rhs values
            rhs[(i-ilower[0])*ny+(j-ilower[1])] =  c(i,j) + dt/8. *(
                    4./(dx*dx) * ( d11_ip12j *  (c(i+1,j) - c(i,j))
                                   -d11_im12j *  (c(i,j)   - c(i-1,j)))
                    + 1./(dx*dy) * ( d12_ip12j *  (c(i,j+1) + c(i+1,j+1) - c(i,j-1)   - c(i+1,j-1))
                                     -d12_im12j *  (c(i,j+1) + c(i-1,j+1) - c(i-1,j-1) - c(i,j-1)))
                    + 1./(dx*dy) * ( d21_ijp12 *  (c(i+1,j) + c(i+1,j+1) - c(i-1,j+1) - c(i-1,j))
                                     -d21_ijm12 *  (c(i+1,j) + c(i+1,j-1) - c(i-1,j)   - c(i-1,j-1)))
                    + 4./(dy*dy) * ( d22_ijp12 *  (c(i,j+1) - c(i,j))
                                     -d22_ijm12 *  (c(i,j) - c(i,j-1)))
                                                   );
        }
    }



    // set RHS to zero for those cells with Dirichlet boundaries!
    for (int m=0; m<N_para_dir; m++) {
        int i=para_dir[m][0];
        int j=para_dir[m][1];
        if (i >= ilower[0] && i <= iupper[0] && j >= ilower[1] && j <= iupper[1]) {
            rhs[(i-ilower[0])*ny+(j-ilower[1])] = 0;
        }
    }

    AssignRHS_Parabolic(rhs);
    delete [] rhs;
}



void MGSolver::AssignRHS_Parabolic(double const * const rhs) {
    int nvalues = nx * ny;
    double *values = (double*) calloc(nvalues, sizeof(double));

    HYPRE_StructVectorSetBoxValues(x, ilower, iupper, values);
    int idx[2];
    for (int k = ilower[0] ; k <= iupper[0]; k++) {
        for (int l = ilower[1]; l <= iupper[1]; l++) {
            idx[0] = k;
            idx[1] = l;
            HYPRE_StructVectorSetValues(b, idx, rhs[(k-ilower[0])*ny + (l-ilower[1])]);
        }
    }
    free(values);
    // this is a collective call finalizing the vector assembly.
    HYPRE_StructVectorAssemble(b);
    HYPRE_StructVectorAssemble(x);
}

