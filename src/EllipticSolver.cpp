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
/* EllipticSolver
 *
 * This class provides an additional layer of abstraction for the
 * elliptic solver. The system, once distcretized, can in theory be solved
 * using any suitable solver (e.g. multigrid, cg, cuda-cg).
 */
//============================================================================

#include "EllipticSolver.h"
#include <iostream>
#include <libconfig.h++>
#include <cmath>
#include "Printer.h"
#include <fstream>
#include "GenSimDefs.h"
#include "Util.h"
#include <omp.h>
#include "mpi.h"
#include "MGSolver.h"


using namespace std;


EllipticSolver::EllipticSolver(SimParams * sp) {
    this->sp = sp;

    // domain
    x_ub  = sp->x_ub;
    x_lb  = sp->x_lb;
    y_ub  = sp->y_ub;
    y_lb  = sp->y_lb;

    dx  = (x_ub - x_lb) / (sp->N1);
    dy  = (y_ub - y_lb) / (sp->N2);
    dA  = dx * dy;

    // rotation angle of base vector e_2
    alpha   = atan(dx/dy);
    // and e_1
    theta   = atan(dy/dx);

    // grid spacing in rotated coordinate system
    dxi     = dx/cos(theta);
    deta    = dy/cos(alpha);

    // precalculate trigonometric functions for better performance
    secamt    = 1./cos(alpha-theta);
    secamt2   = secamt * secamt;
    cosa2     = cos(alpha) * cos(alpha);
    sina2     = sin(alpha) * sin(alpha);
    costsina  = cos(theta) * sin(alpha);
    cosasint  = cos(alpha) * sin(theta);
    cost2     = cos(theta) * cos(theta);
    sint2     = sin(theta) * sin(theta);
    cosa      = cos(alpha);
    sint      = sin(theta);
    cost      = cos(theta);
    sina      = sin(alpha);

    // default weighting factor between the two stencils
    SetNu(sp->nu);

    last_parabolic_its = 0;
    Initialize();
}

int EllipticSolver::GetLastIts() {
    return last_its;
}

int EllipticSolver::GetLastParabolicIts() {
    return last_parabolic_its;
}


double EllipticSolver::GetAvgIts() {
    if (calls > 0)
        return ((double)total_its)/calls;
    else
        return 0;
}


EllipticSolver::~EllipticSolver() {}



void EllipticSolver::Initialize() {

    // precalculate this factor for better performance
    mob1over4 = pow (sp->mu_res/sp->mu_inj , 0.25);
    last_its = 0;
    total_its = 0;
    calls = 0;
}


bool EllipticSolver::Solve(Array2d & p, Array2d &c, double * q) {
    bool ret = true;
    calls++;

    // setup and solve system
    ret = ret && SolveLinearMultigrid(c,sp->c_bounds, p,sp->p_bounds, q);

    total_its += last_its;
    return ret;
}



/******************************************************************************
  We want to solve d²/dx² p + d²/dy² p =  -q * my_res (1)  satisfying that the
  integral of q over Omega vanishes.

  Here, we require p(x,y) = cos(x) * cox(y), satisfying the boundary conditions
  and calcualte q, v and K such that p is a solution to eqn. (1).
******************************************************************************/


bool EllipticSolver::TestNeumann(GnuPlotter *gp) {

    sp->btype = NEUMANN;
    bool stat = true;

    Array2d p_ex;
    Array2d p;
    Array2d c;
    Array3d v;
    Array3d v_ex;

    double * q;
    sp->N1 = sp->init_x_cells;
    sp->N2 = sp->init_y_cells;

    double *normarr = new double[sp->levels];
    double *xarr = new double[sp->levels];

    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    ofstream norm_stream;
    bool fstat = true;
    if (myid == 0) {
        norm_stream.open(sp->norm_file.c_str());
        fstat = norm_stream.is_open();
        norm_stream.precision(30);
    }
    if (fstat) {
        for (int k = 0; k < sp->levels; k++) {
            if (myid == 0)
                cout << "-> (N1,N2):\t(" << sp->N1 << " , " << sp->N2 << ")" ;

            sp->dx = (sp->x_ub - sp->x_lb)/sp->N1;
            sp->dy = (sp->y_ub - sp->y_lb)/sp->N2;
            sp->M = sp->N1*sp->N2;

            EllipticSolver * mg_test = new MGSolver(sp);
            Initialize();

            // pressure
            p.resize(sp->N1 + 2 *sp->p_bounds,sp->N2 + 2 * sp->p_bounds);
            Idx2d idx;
            idx[0] = -sp->p_bounds;
            idx[1] = -sp->p_bounds;
            p.reindexSelf(idx);
            p_ex.resize(sp->N1 + 2 * sp->p_bounds,sp->N2 + 2 * sp->p_bounds);
            p_ex.reindexSelf(idx);

            // concentration
            c.resize(sp->N1 + 2 * sp->p_bounds,sp->N2 + 2 * sp->p_bounds);
            c.reindexSelf(idx);

            //velocity
            v.resize((sp->N1+1) + 2 * sp->v_bounds,(sp->N2+1) + 2 * sp->v_bounds,2);
            Idx3d idx3d;
            idx3d[0] = -sp->v_bounds;
            idx3d[1] = -sp->v_bounds;
            idx3d[2] = 0;
            v.reindexSelf(idx3d);

            v_ex.resize((sp->N1+1) + 2 * sp->v_bounds,(sp->N2+1) + 2 * sp->v_bounds,2);
            v_ex.reindexSelf(idx3d);

            p = 0;
            p_ex = 0;
            v = 0;
            v_ex = 0;
            c = 0;

            // source term
            q = new double[sp->M];

            // initialize Arrays and calculate exact solutions p_ex and v_ex
            InitializeEllipticNeumannTest(c,p,p_ex, q, v_ex);

            // fill ghost cells
            FillGhostCellsNeumannScalar(c, sp->c_bounds);

            // setup and solve syste
            stat = stat && mg_test->Solve(p,c,q);
            if (myid == 0)
                cout << "\tITS:\t" << mg_test->last_its << endl;

            // calculate the velocity field
            mg_test->CalcVelocities(p,v,c);

            // calculate pressure norms
            double L2_p_err, p_sum=0;
            for (int i=0; i<sp->N1; i++) {
                for (int j=0; j<sp->N2; j++)
                    p_sum += fabs(p(i, j) - p_ex(i, j));
            }

            p_sum *= sp->dx * sp->dy;
            L2_p_err = p_sum;

            // calculate velocity norms
            double L2_u_err, u_sum=0;
            double L2_v_err, v_sum=0;

            for (int i=0; i<=sp->N1; i++) {
                for (int j=0; j<sp->N2; j++)
                    u_sum += fabs((v_ex(i, j,0)  - v(i, j,0)));
            }

            for (int i=0; i<sp->N1; i++) {
                for (int j=0; j<=sp->N2; j++)
                    v_sum += fabs((v_ex(i, j, 1) - v(i, j, 1)));
            }

            u_sum *=sp->dx * sp->dy;
            v_sum *=sp->dx * sp->dy;
            L2_u_err = u_sum;
            L2_v_err = v_sum;

            // write results to file for first level
            if (myid == 0 &&  k==0) {
                Printer pr;
                pr.PrintDToPlotFile(   sp->num_result_file, p,		          sp->x_lb, 	sp->y_lb, sp->dx,     sp->dy,   sp->p_bounds);
                pr.PrintDToPlotFile(   sp->ex_result_file,  p_ex, 	        sp->x_lb,	  sp->y_lb, sp->dx,     sp->dy,   sp->p_bounds);
                pr.PrintVelToPlotFile( sp->num_u_file,	    sp->num_v_file, v, 			    sp->x_lb, sp->y_lb,   sp->dx,   sp->dy,sp->v_bounds);
                pr.PrintVelToPlotFile( sp->ex_u_file,		    sp->ex_v_file, 	v_ex, 	    sp->x_lb,sp-> y_lb,   sp->dx,   sp->dy,sp->v_bounds);
                cout << "\t(Writing files " << sp->num_result_file << ", " << sp->ex_result_file << ", "  << sp->num_u_file << ", "  << sp->ex_u_file << ")"<< endl;
            }

            if (myid == 0) {
                norm_stream << log(sp->dx) << "\t" << log(L2_p_err) << "\t" << log(L2_u_err) << "\t" << log(L2_v_err) << endl;
                normarr[k] = log(L2_p_err);
                xarr[k] = log(sp->dx);
            }
            delete [] q;
            delete mg_test;

            // double the resolution
            sp->N1 *= 2;
            sp->N2 *= 2;
        }
    }
    else {
        cerr << "Cannot write norm-file! Convergence test aborted..." << endl;
        stat = false;
    }

    if (myid == 0) {
        norm_stream.close();

        // perform fit
        if (sp->levels > 1) {
            if (sp->gnuplot) {
                gp->PlotConvergence( sp->norm_file);
            }
        }
    }
    delete []normarr;
    delete []xarr;
    return stat;
}


/******************************************************************************
 Initialization for test cases 1-4.
******************************************************************************/


double EllipticSolver::qTest1(double x, double y) {
    return 2./sp->mu_res * cos(x) * cos(y);
}

double EllipticSolver::qTest2(double x, double y) {
    double c1 = -0.3;
    double c2 = -0.3;
    double a1 = exp(c1*x);
    double a2 = exp(c2*y);
    double a1_x = c1* exp(c1*x);
    double a2_y = c2* exp(c2*y);
    return 1./sp->mu_res * ((a1 + a2) * cos(x) * cos(y)  + a1_x * sin(x) * cos(y) + a2_y * cos(x) * sin(y));
}

double EllipticSolver::qTest3(double x, double y) {
    double t = 2.;
    double a1 = 2. + sin(t*x) * sin(y);
    double a2 = 2. + sin(t*x) * sin(y);
    double a1_x = t*cos(t*x) * sin(y);
    double a2_y = cos(y) * sin(t*x);
    return 1./sp->mu_res * ((a1 + a2) * cos(x) * cos(y)  + a1_x * sin(x) * cos(y) + a2_y * cos(x) * sin(y));
}

double EllipticSolver::qTest4(double x, double y) {
    double t = 2.;
    double a1 = 2. + sin(t*x) * sin(y);
    double a2 = 2. + sin(x) * sin(t*y);
    double a1_x = t*cos(t*x) * sin(y);
    double a2_y = t* cos(t*y) * sin(x);
    return 1./sp->mu_res * ((a1 + a2) * cos(x) * cos(y)  + a1_x * sin(x) * cos(y) + a2_y * cos(x) * sin(y));
}

double EllipticSolver::pTest(double x, double y) {
    return cos(x) * cos(y);
}



void EllipticSolver::InitializeEllipticNeumannTest(Array2d &c ,Array2d &p, Array2d &p_ex, double *q, Array3d & v_ex) {
    double x,y;
    double c1,c2,a1,a2,a1_x,a2_y,t;

    // initialize exact solution for the velocity
    double x_u, y_u;
    for (int i=0; i<=sp->N1; i++) {
        for (int j=0; j<sp->N2; j++) {
            x_u =sp->x_lb + i * sp->dx;
            y_u = sp->y_lb + (j+0.5) * sp->dy;
            v_ex(i, j,0) = K_11(x_u,y_u)/sp->mu_res * sin(x_u) * cos(y_u);
        }
    }

    // initialize exact solution for the velocity
    double x_v, y_v;
    for (int i=0; i<sp->N1; i++) {
        for (int j=0; j<=sp->N2; j++) {
            x_v = sp->x_lb + (i+0.5)  *sp->dx;
            y_v = sp->y_lb + j * sp->dy;
            v_ex(i, j, 1) = K_22(x_v,y_v)/sp->mu_res * cos(x_v) * sin(y_v);
        }
    }


    // initialize exact solution for pressure and initialize q
    for (int i=0; i<sp->N1; i++) {
        for (int j=0; j<sp->N2; j++) {
            x = sp->x_lb + (i+0.5)  * sp->dx;
            y = sp->y_lb + (j+0.5) * sp->dy;
            p_ex(i,j) = pTest(x,y);
            c(i, j) = 0;
            p(i, j) = 0;

            switch (sp->testno) {
            case 1:
                q[i * sp->N2 + j]  = qTest1(x,y);
                break;
            case 2:
                q[i * sp->N2 + j]  = qTest2(x,y);
                break;
            case 3:
                q[i * sp->N2 + j]  = qTest3(x,y);
                break;
            case 4:
                q[i * sp->N2 + j]  = qTest4(x,y);
                break;
            default:
                cerr << "EllipticSolver::TestNeumann::Initialization> Unknown test no.!" << endl;
                exit(1);
            }
        }
    }
}



/******************************************************************************
 Diagonal components of K
******************************************************************************/

double EllipticSolver::K_11(double x, double y) {
    double ret;
    double t = 2.;
    if (sp->is_test) {
        switch (sp->testno) {
        case 1:
            ret = 1.;
            break;
        case 2:
            ret = exp(-0.3*x);
            break;
        case 3:
            ret = 2. + sin(t*x) * sin(y);
            break;
        case 4:
            ret = 2. + sin(t*x) * sin(y);
            break;
        default:
            cerr << "Unknown test no!" << endl;
            exit(1);
        }
    }
    else
        ret = 1.;
    return ret;
}

double EllipticSolver::K_22(double x, double y) {
    double ret;
    double t = 2.;
    if (sp->is_test) {
        switch (sp->testno) {
        case 1:
            ret = 1.;
            break;
        case 2:
            ret = exp(-0.3*y);
            break;
        case 3:
            ret = 2. + sin(t*x) * sin(y);
            break;
        case 4:
            ret = 2. + sin(x) * sin(t*y);
            break;
        default:
            cerr << "Unknown test no!" << endl;
            exit(1);
        }
    }
    else
        ret = 1.;
    return ret;
}



/******************************************************************************
 Calculate the mobility tensor (M=K/mu) at the cell interfaces.
******************************************************************************/

double EllipticSolver::KoverMu_11_im12j(double c_im1j,double c_ij, int i, int j) {
    double x=  sp->x_lb + (i+0.5) * sp->dx;
    double y = sp->y_lb + (j+0.5) * sp->dy;
    double K11_im1j 	= 	K_11(x - sp->dx,  y);
    double K11_ij 			= 	K_11(x ,  y);
    double mu_ij =		mu(c_ij);
    double mu_im1j = 	mu(c_im1j);

    double ret;
    switch (sp->avgproc) {
    case HARMONIC:
        ret =  	2./ ( mu_im1j / K11_im1j   + mu_ij / K11_ij );
        break;
    case MEANCONC:
        ret = 	1./2 * ( K11_im1j / mu_im1j  + K11_ij/mu_ij);
        break;
    default:
        cerr << "KoverMu_11_im12j: Unknown averaging procedure!" << endl;
        exit(1);
        break;
    }
    return ret;
}

double EllipticSolver::KoverMu_22_ijm12(double c_ijm1,double c_ij, int i, int j) {
    double x=  sp->x_lb + (i+0.5) * sp->dx;
    double y = sp->y_lb + (j+0.5) * sp->dy;
    double K22_ijm1 	= 	K_22(x,  y - sp->dy);
    double K22_ij 			= 	K_22(x ,  y);
    double mu_ij = 		mu(c_ij);
    double mu_ijm1 =  mu(c_ijm1);

    double ret;
    switch (sp->avgproc) {
    case HARMONIC:
        ret = 2./ ( mu_ijm1 / K22_ijm1 + mu_ij / K22_ij);
        break;
    case MEANCONC:
        ret = 1./2 * (K22_ijm1 / mu_ijm1 + K22_ij/mu_ij);
        break;
    default:
        cerr << "KoverMu_22_ijm12: Unknown averaging procedure!" << endl;
        exit(1);
        break;
    }
    return ret;
}

void EllipticSolver::SetNu(double nu) {
    if (nu >= 0. && nu <= 1.) {
        this->nu = nu;
    }
    else {
        cerr << "EllipticSolver::SetNu> Invalid value for nu: " << nu << endl;
        exit(1);
    }
}


/******************************************************************************
 Calculate the parallel-stencil contribution to the velocity.
******************************************************************************/

double EllipticSolver::GetParallelU1_im12j(const Array2d &c, const Array2d &p, int i, int j)  {
    return -KoverMu_11_im12j((double)c(i-1, j), (double)c(i, j), i, j)  * (p(i, j)-p(i-1, j))/dx;
}

double EllipticSolver::GetParallelU2_ijm12(const Array2d &c, const Array2d &p, int i, int j)  {
    return -KoverMu_22_ijm12  ((double)c(i, j-1), (double)c(i, j), i, j)  *  (p(i, j)-p(i, j-1))/dy;
}

/******************************************************************************
 Calculate relevant components of M for the diagonal contribution to the
 velocity.
******************************************************************************/

void EllipticSolver::M_im12j(double & M11_im12j, double & M12_im12j, double &M21_im12j, double &M22_im12j, const Array2d&c, int i, int j) {
    const double K11_im12jm12			=	KoverMu_11_ip12jp12(i-1, j-1, c(i-1, j-1), c(i, j-1), c(i-1, j), c(i, j));
    const double K22_im12jm12			= KoverMu_22_ip12jp12(i-1, j-1, c(i-1, j-1), c(i, j-1), c(i-1, j), c(i, j));

    M11_im12j   = secamt2 * ( K11_im12jm12 *  cosa2      +  K22_im12jm12    *     sina2 );
    M12_im12j   = secamt2 * ( K22_im12jm12 *  costsina   -  K11_im12jm12    *     cosasint);
    M21_im12j   = secamt2 * ( K22_im12jm12 *  costsina   -  K11_im12jm12    *     cosasint);
    M22_im12j   = secamt2 * ( K11_im12jm12 *  sint2      +  K22_im12jm12    *     cost2 );
}


void EllipticSolver::GetDiagonalW_im12j(double &w1_im12j, double &w2_im12j, const Array2d &c, const Array2d &p, const int i,const int j) {
    double d1dp_im12j = 1./dxi 		* (p(i, j) - p(i-1, j-1));
    double d2dp_im12j = 1./deta 	* (p(i-1, j) - p(i, j-1));
    double M11_im12j, M12_im12j,M21_im12j,M22_im12j;

    M_im12j( M11_im12j, M12_im12j,M21_im12j,M22_im12j, c, i,j);

    w1_im12j = - d1dp_im12j * M11_im12j - d2dp_im12j * M12_im12j;
    w2_im12j = - d1dp_im12j * M21_im12j - d2dp_im12j * M22_im12j;
}

/******************************************************************************
 Transform (x,y) to (xi,eta)  <=>  multiply by Rinv.
******************************************************************************/

void EllipticSolver::XYToXiEta(const double v1, const double v2, double &w1, double &w2) {
    w1  = cosa * secamt * v1  + sina * secamt * v2;
    w2  = -secamt * sint * v1 + cost * secamt * v2;
}

/******************************************************************************
 Transform (xi,eta) to (x,y)  <=> multiply by R.
******************************************************************************/

void EllipticSolver::XYFromXiEta(double &v1, double& v2, const double w1, const double w2) {
    v1 = cost * w1 - sina * w2;
    v2 = sint * w1 + cosa * w2;
}


/******************************************************************************
 Calculate velocity in (xi,eta) coordinates, transform it to (x,y) and finally
 interpolate.
******************************************************************************/

double EllipticSolver::GetDiagonalU1_im12j(const Array2d &c, const Array2d &p, int i, int j)  {

    // calculate velocity w in (xi,eta) coordinates at im12j and ijp12
    // im12j   	(xi,eta)
    double w1_im12j, w2_im12j;
    GetDiagonalW_im12j(w1_im12j, w2_im12j, c,p,i,j);

    // ijp12		(xi,eta)
    double w1_ijp12, w2_ijp12;
    GetDiagonalW_im12j(w1_ijp12, w2_ijp12, c,p,i,j+1);

    // transform velocities w back to (x,y)-coordiantes by multiplication with R.
    double u1_im12jm12, u2_im12jm12;
    XYFromXiEta(u1_im12jm12, u2_im12jm12, w1_im12j, w2_im12j);

    double u1_im12jp12, u2_im12jp12;
    XYFromXiEta(u1_im12jp12, u2_im12jp12, w1_ijp12, w2_ijp12);

    // interpolate linearly
    return 0.5 * (u1_im12jp12 + u1_im12jm12);

}


double EllipticSolver::GetDiagonalU2_ijm12(const Array2d &c, const Array2d &p, int i, int j)  {
    // calculate velocity w in (xi,eta) coordinates at im12j and ijm12
    // im12j   	(xi,eta)
    double w1_im12j, w2_im12j;
    GetDiagonalW_im12j(w1_im12j, w2_im12j, c,p,i,j);

    // ijm12		(xi,eta)
    double w1_ijm12, w2_ijm12;
    GetDiagonalW_im12j(w1_ijm12, w2_ijm12, c,p,i+1,j);

    // transform velocities w back to (x,y)-coordiantes by multiplication with R.
    double u1_im12jm12, u2_im12jm12;
    XYFromXiEta(u1_im12jm12, u2_im12jm12, w1_im12j, w2_im12j);

    double u1_ip12jm12, u2_ip12jm12;
    XYFromXiEta(u1_ip12jm12, u2_ip12jm12, w1_ijm12, w2_ijm12);

    // interpolate linearly
    return 0.5 * (u2_im12jm12 + u2_ip12jm12);

}




/******************************************************************************
 This method fills the velocity array v using the pressure array p
 NOTE: 		u_00 and v_00 don't coincide in space!
 	       	u_00 is at (0,0.5) with regards to the origin
         	v_00 is at (0.5, 0) with regards to the origin
         	p_00 is at (0.5,0.5) w.r.t.t.o
 (u,v) = -K/mu * grad(p)

 The total velocity consists of a parallel and a diagonal contribution!
******************************************************************************/

void EllipticSolver::CalcVelocities(Array2d & p, Array3d &v, Array2d& c) {

    // Fill the pressure ghost cells first
    FillGhostCellsNeumannScalar(p, sp->p_bounds);

#ifdef OMPELLIPTIC
    #pragma omp parallel for
#endif
    for (int i=0; i<=sp->N1; i++) {
        for (int j=0; j<sp->N2; j++)
            v(i, j,0) = nu *  GetParallelU1_im12j(c,p,i,j) +    	(1.-nu) * GetDiagonalU1_im12j(c,p,i,j)  ;
    }

#ifdef OMPELLIPTIC
    #pragma omp parallel for
#endif
    for (int i=0; i<sp->N1; i++) {
        for (int j=0; j<=sp->N2; j++)
            v(i, j, 1) = nu *GetParallelU2_ijm12(c, p, i, j) +	 	(1.-nu) *  GetDiagonalU2_ijm12(c,p,i,j)  ;
    }

    // Fill the velocity boundary cells
    FillGhostCellsNeumannVelocity(v, sp->v_bounds);
}



/******************************************************************************
 Calculate mu(c).  NOTE: Don't use pow! This method is invoked many times and
 pow is very expensive. This gives better performance.
******************************************************************************/

double EllipticSolver::mu(double c) {
    double x = (1. - c + mob1over4 * c);
    return 1. / (x * x * x * x) * sp->mu_res;
}

double EllipticSolver::KoverMu_11_ip12jp12(int i, int j, double c_ij, double c_ip1j, double c_ijp1, double c_ip1jp1) {
    double x=  x_lb + (i+0.5) * dx;
    double y = y_lb + (j+0.5) * dy;

    double K11_ip1j      =   K_11(x + dx,  y);
    double K11_ij        =   K_11(x ,  y);
    double K11_ijp1      =   K_11(x ,  y + dy);
    double K11_ip1jp1    =   K_11(x +dx ,  y + dy);
    double mu_ip1j       = mu(c_ip1j);
    double mu_ij         =   mu( c_ij);
    double mu_ijp1       =   mu(c_ijp1);
    double mu_ip1jp1     = mu( c_ip1jp1);

    double ret;
    switch (sp->avgproc) {
    case HARMONIC:
        ret =  4./ ( mu_ip1j / K11_ip1j   + mu_ij / K11_ij + mu_ijp1/ K11_ijp1 + mu_ip1jp1/K11_ip1jp1);
        break;
    case MEANCONC:
        ret = 1./4 * (K11_ip1j/mu_ip1j + K11_ij/mu_ij + K11_ijp1/mu_ijp1 + K11_ip1jp1/mu_ip1jp1);
        break;
    default:
        cerr << "KoverMu_11_ip12jp12: Unknown averaging procedure!" << endl;
        exit(1);
        break;
    }
    return ret;
}


double EllipticSolver::KoverMu_22_ip12jp12(int i, int j, double c_ij, double c_ip1j, double c_ijp1, double c_ip1jp1) {
    double x=  x_lb + (i+0.5) * dx;
    double y = y_lb + (j+0.5) * dy;

    double K22_ip1j         =   K_22(x + dx,  y);
    double K22_ij           =   K_22(x ,  y);
    double K22_ijp1         =   K_22(x ,  y + dy);
    double K22_ip1jp1       =   K_22(x +dx ,  y + dy);
    double mu_ip1j        = mu(c_ip1j);
    double mu_ij          =   mu(c_ij);
    double mu_ijp1        =   mu(c_ijp1);
    double mu_ip1jp1      = mu( c_ip1jp1);
    double ret;

    switch (sp->avgproc) {
    case HARMONIC:
        ret =  4./ ( mu_ip1j / K22_ip1j   + mu_ij / K22_ij + mu_ijp1/ K22_ijp1 + mu_ip1jp1/K22_ip1jp1);
        break;
    case MEANCONC:
        ret = 1./4 * (K22_ip1j/mu_ip1j + K22_ij/mu_ij + K22_ijp1/mu_ijp1 + K22_ip1jp1/mu_ip1jp1);
        break;
    default:
        cerr << "KoverMu_22_ip12jp12: Unknown averaging procedure!" << endl;
        exit(1);
        break;
    }
    return ret;
}

