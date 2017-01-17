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
/* Util.cpp
 *
 * This file provides a variety of routines invoked from different places.
 */
//============================================================================


#include "Util.h"
#include "SimParams.h"
#include <cmath>
#include <sys/stat.h>
#include <fstream>
#include "GenSimDefs.h"

extern "C" void   	dgesv_	(int *n, 	int *nrhs, double *a, int *lda, int *ipiv,	double *b, int *ldb, int *info );


/*
 * Initialize top-hat profile (for advection test).
 */
double Tophat(SimParams const * const sp, double x, double y) {
    double width = 0.3 * (sp->x_ub - sp->x_lb);
    double height = 0.3 * (sp->y_ub - sp->y_lb);

    double xc = (sp->x_ub - sp->x_lb) / 2.;
    double yc = (sp->y_ub - sp->y_lb) / 2.;
    if ( fabs(x-xc) <= width/2. && fabs(y-yc) <= height/2.)
        return 1.;
    return 0.;
}

/*
 * For initialization with a gaussian bell (for advection test).
 * Return corresponding value for the point (x,y).
 */
double Gaussbell(SimParams const * const sp, double x, double y) {
    double xc = (sp->x_ub - sp->x_lb) / 2.;
    double yc = (sp->y_ub - sp->y_lb) / 2.;

    // Integral of Exp[-64 x*x] over [x-dx/2,x+dx/2]
    double intx = 1./16 * sqrt(M_PI) * ( -  erf(8. * xc - 4. * sp->dx - 8. * x)        + erf( 8. * xc + 4. * sp->dx -8. * x));
    double inty = 1./16 * sqrt(M_PI) * ( -  erf(8. * yc - 4. * sp->dy - 8. * y)        + erf( 8. * yc + 4. * sp->dy -8. * y));
    return 1./(sp->dx * sp->dy) * intx * inty;
}

/*
 * For initialization with a diagonal sine wave (for advection test).
 * Return corresponding value for the point (x,y).
 */
double DiagonalSinewave(SimParams const * const sp,double x, double y) {
    double xc = (sp->x_ub - sp->x_lb) / 2.;
    double yc = (sp->y_ub - sp->y_lb) / 2.;
    return 1. +0.2 * sin(4.*M_PI * ((x-xc) + (y-yc)));
}

/*
 * For initialization with a sine wave (for advection test).
 * Return corresponding value for the point (x,y).
 */
double Sinewave(SimParams const * const sp,double x, double y) {
    return  1. + 1./( sp->dx * sp->dy * 4. * M_PI * M_PI) *
            (cos(2. * M_PI * (x+sp->dx/2.)) - cos(2. * M_PI * (x-sp->dx/2.))) *
            (cos(2. * M_PI * (y+sp->dy/2.)) - cos(2. * M_PI * (y-sp->dy/2.)));
}

/*
 * For initialization with a horizontal sine wave (for advection test).
 * Return corresponding value for the point (x,y).
 */
double HorizontalSine(SimParams const * const sp,double x, double y) {
    return 1. +  1./(sp->dx) * (-1./(2. * M_PI)) * ( cos( 2 * M_PI * (x + sp->dx / 2.)) - cos( 2 * M_PI * (x - sp->dx / 2.)));
}

/*
 * For initialization with a vertical sine wave (for advection test).
 * Return corresponding value for the point (x,y).
 */
double VerticalSine(SimParams const * const sp,double x, double y) {
    return 1. +  1./(sp->dy) * (-1./(2. * M_PI)) * ( cos( 2 * M_PI * (y + sp->dy / 2.)) - cos( 2 * M_PI * (y - sp->dy / 2.)));
}

/*
 * Calculate L2-norm for general exact solution specified by the function fptr.
 */
double CalcAdvectionL2Norm(SimParams const * const sp,const Array2d& c, double t, double (*fptr)(SimParams const*const,double, double)) {
    double l2sum = 0.;
    for (int i=0; 		i<sp->N1; i++) {
        for (int j=0; 	j<sp->N2; j++) {
            double x = sp->x_lb + (i + 0.5) * sp->dx;
            double y = sp->y_lb + (j + 0.5) * sp->dy;
            double ex = fptr(sp,x,y);
            l2sum += (ex - c(i, j)) * (ex - c(i, j));
        }
    }
    return  sqrt(sp->dx * sp->dy *l2sum);
}

/*
 * Calculate L2-norm for general exact solution specified by the function fptr.
 */
double CalcAdvectionL1Norm(SimParams const * const sp,const Array2d& c, double t, double (*fptr)(SimParams const*const,double, double)) {
    double l1sum = 0.;
    for (int i=0; 		i<sp->N1; i++) {
        for (int j=0; 	j<sp->N2; j++) {
            double x = sp->x_lb + (i + 0.5) * sp->dx;
            double y = sp->y_lb + (j + 0.5) * sp->dy;
            double ex = fptr(sp,x,y);
            l1sum += fabs((double)(ex - c(i, j)));
        }
    }
    l1sum *= sp->dx * sp->dy;
    return l1sum;
}

/*
 *	Create folders for output files.
 */

void SetupDirectories(SimParams const * const sp) {
    struct stat st;
    if (sp->log_to_file) {
        // Create logging-directory if necessary
        if(!stat(sp->log_dir.c_str(),&st) == 0) {
            umask(0);
            if (mkdir(sp->log_dir.c_str(), 0755) == 0)
                cout << "Directory " << sp->log_dir << " was created..." << endl;
            else {
                cerr << "Error: could not create " << sp->log_dir << endl;
                exit(1);
            }
        }
    }

    // Create results-directory if necessary
    if(!stat(sp->results_dir.c_str(),&st) == 0) {
        umask(0);
        if (mkdir(sp->results_dir.c_str(), 0755) == 0)
            cout << "Directory " << sp->results_dir << " was created..." << endl;
        else {
            cerr << "Error: could not create " << sp->results_dir << endl;
            exit(1);
        }
    }
}



/*
 * Fill ghost cells to impose shifted Neumann boundary conditions.
 */

void FillGhostCellsShiftedNeumannScalar(SimParams const * const sp,Array2d &arr) {

    for (int i=0; i<sp->N1; i++) {

        // bottom
        arr(i, -1) = arr(i, 1);

        // top
        arr(i, sp->N2) = arr(i, sp->N2-2);
    }
    for (int j=0; j<sp->N2; j++) {

        // left
        arr(-1, j) = arr(1, j);

        // right
        arr(sp->N1, j) = arr(sp->N1-2, j);
    }
    arr(-1, -1) = arr(1, 1);
    arr(-1, sp->N2) = arr(1, sp->N2-2);
    arr(sp->N1, sp->N2) = arr(sp->N1-2, sp->N2-2);
    arr(sp->N1, -1) = arr(sp->N1-2, 1);
}


/*
 * Fill ghost cells to impose periodic boundary conditions.
 */

void FillGhostCellsPeriodicScalar(Array2d &arr, int n) {

    // start and end index for both dimensions
    const int x0 = arr.base()[0]+n;
    const int y0 = arr.base()[1]+n;
    const int xN = x0 +arr.extent()[0] - 1 -2*n;
    const int yN = y0 + arr.extent()[1] - 1 - 2*n;

    // TOP & BOTTOM  (excluding corners)
    for (int i=0; i<=xN; i++) {
        for (int j = 0; j < n; j++) {
            arr(i, yN+j+1) = arr(i, j);
            arr(i, -(j+1)) = arr(i, yN-j);
        }
    }

    // LEFT & RIGHT  (excluding corners)
    for (int i = 0; i < n; i++) {
        for (int j=0; j<=yN; j++) {
            arr(-(i+1), j) = arr(xN-i, j);
            arr(xN+i+1, j) = arr(i, j);
        }
    }

    // CORNERS
    for (int i = 0; i < n; i++) {
        for (int j=0; j<n; j++) {

            // BOTTOM LEFT
            arr(-(i+1), -(j+1)) = arr(xN-i, yN-j);

            // BOTTOM RIGHT
            arr(xN+i+1, -(j+1)) = arr(i, yN-j);

            // TOP LEFT
            arr(-(i+1), yN+j+1) = arr(xN-i, j);

            // TOP RIGHT
            arr(xN+i+1, yN+j+1) = arr(i, j);
        }
    }
}

/*
 * Check for Neumann boundary type.
 */
void CheckNeumannScalar(SimParams const * const sp,const Array2d &arr) {
    double checksum = 0.;
    for (int i=-sp->c_bounds; i < sp-> N1 + sp->c_bounds ; i++) {
        checksum += fabs( arr(i,0) - arr(i, -1));
        checksum += fabs( arr(i, sp->N2) - arr(i, sp->N2-1));
    }

    for (int j=-sp->c_bounds; j < sp-> N2 + sp->c_bounds ; j++) {
        checksum += fabs( arr(0,j) - arr(-1, j));
        checksum += fabs( arr(sp->N1, j) - arr(sp->N1-1, j));
    }

    if (checksum != 0) {
        cerr << "BOUNDARIES VALUES ARE NOT CORRECT!:  " << checksum  << endl;
        exit(1);
    }
}

/*
 * Fill ghost layers for Neumann boundaries.
 */

void FillGhostCellsNeumannScalar(Array2d &arr, int n) {

    // start and end index for both dimensions
    const int x0 = arr.base()[0]+n;
    const int y0 = arr.base()[1]+n;
    const int xN = x0 +arr.extent()[0] - 1 -2*n;
    const int yN = y0 + arr.extent()[1] - 1 - 2*n;

    // TOP & BOTTOM  (excluding corners)
    for (int i=0; i<=xN; i++) {
        for (int j = 0; j < n; j++) {
            arr(i, yN+j+1) = arr(i, yN-j);
            arr(i, -(j+1)) = arr(i, j);
        }
    }

    // LEFT & RIGHT  (excluding corners)
    for (int i = 0; i < n; i++) {
        for (int j=0; j<=yN; j++) {
            arr(-(i+1), j) = arr(i, j);
            arr(xN+i+1, j) = arr(xN - i, j);
        }
    }

    // CORNERS
    for (int i = 0; i < n; i++) {
        for (int j=0; j<n; j++) {

            // BOTTOM LEFT
            arr(-(i+1), -(j+1)) = arr(i, j);

            // BOTTOM RIGHT
            arr(xN+i+1, -(j+1)) = arr(xN-i, j);

            // TOP LEFT
            arr(-(i+1), yN+j+1) = arr(i, yN-j);

            // TOP RIGHT
            arr(xN+i+1, yN+j+1) = arr(xN-i, yN-j);
        }
    }
}


/*
 * Fill velocity ghost cells for Neumann boundaries.
 */

void FillGhostCellsNeumannVelocity(Array3d &v, int n) {

    // Assuming n boundary cells in both dimensions
    const int x0 = v.base()[0]+n;
    const int y0 = v.base()[1]+n;
    const int xN = x0 +v.extent()[0] - 1 -2*n;
    const int yN = y0 + v.extent()[1] - 1 - 2*n;

    // deal with u-components first: v[:][:][0] 

    // bottom and top
    for (int i=0; i <= xN; i++) {
        for (int j = 0; j < n; j++) {

            // bottom
            v(i, -(j+1),0) 	= v(i, j,0);

            // top
            v(i, yN+j,0) 	= v(i, yN-1-j,0);
        }
    }

    // left top corner and left bottom corner
    for (int i=0; i < n; i++) {
        for (int j = 0; j < n; j++) {

            // top
            v(-(i+1), yN+j,0) 	=	-v(i+1, yN+j,0);

            // bottom
            v(-(i+1), -(j+1),0) 	= -v(i+1, -(j+1),0);
        }
    }

    // right top corner and right bottom corner
    for (int i=0; i < n; i++) {
        for (int j = 0; j < n; j++) {

            // top
            v(xN+i+1, yN+j,0) 		= -v(xN-1-i, yN+j,0);

            // bottom
            v(xN+i+1, -(j+1),0) 	= -v(xN-i-1, -(j+1),0);
        }
    }

    // left and right
    for (int i=0; i < n; i++) {
        for (int j = 0; j < yN; j++) {

            // left
            v(-(i+1), j,0) 		= -v(i+1, j,0);

            // right
            v(xN+i+1, j,0) 	= -v(xN-i-1, j,0);
        }
    }

    // set v-components: v[:][:][1] 

    // left and right
    for (int i=0; i < n; i++) {
        for (int j = 0; j <= yN; j++) {

            // left
            v(-(i+1), j, 1) 		= v(i, j, 1);

            // right
            v(xN+i, j, 1) 		= v(xN-i-1, j, 1);
        }
    }

    // bottom and top
    for (int i=0; i < xN; i++) {
        for (int j = 0; j < n; j++) {

            // bottom
            v(i, -(j+1), 1) 		= -v(i, j+1, 1);

            // top
            v(i, yN+j+1, 1) 	= -v(i, yN-1-j, 1);
        }
    }

    // left top corner and left bottom corner
    for (int i=0; i < n; i++) {
        for (int j = 0; j < n; j++) {

            // top
            v(-(i+1), yN+1+j, 1) 	= v(i, yN+1+j, 1);

            // bottom
            v(-(i+1), -(j+1), 1) 		= v(i, -(j+1), 1);
        }
    }

    // right top corner and right bottom corner
    for (int i=0; i < n; i++) {
        for (int j = 0; j < n; j++) {

            // top
            v(xN+i, yN+1+j, 1) 		=	v(xN-1-i, yN+j+1, 1);

            // bottom
            v(xN+i, -(j+1), 1) 		= v(xN-1-i, -(j+1), 1);
        }
    }
}

/*
 * Calculate the exact volum fraction for initialization with a circular well.
 */

double ExactVolFracCircle(double x1, double y1, double r, double dx, double dy) {
    if (x1 * x1 + y1*y1 < r*r) {
        double y2 = y1+dy;
        double x2 = x1+dx;
        if (  x2 * x2 + y2 * y2> r*r) {
            double xs = sqrt(r*r - y1 *y1);
            double ys = sqrt(r*r - x1 *x1);

            double area = 0;
            double xub, xlb, yub;
            if ( ys <= y2 && xs <= x2 ) {
                xlb = x1;
                xub = xs;
                yub = ys;
            }
            else if (ys < y2 && xs >x2) {
                yub = ys;
                xlb = x1;
                xub = x2;
            }
            else if (ys > y2 && xs  <x2) {
                yub = y2;
                xlb = sqrt(r*r - yub * yub);
                xub = xs;
                area += (xlb - x1) * dy;
            }
            else {
                yub = y2;
                xlb = sqrt(r*r - yub * yub);
                xub = x2;
                area += dy * (xlb - x1);
            }
            double t1 = asin(xlb/r);
            double t2 = asin(xub/r);
            area += r*r/2 * (  sin(2*t2)/2. + t2 - 		 sin(2*t1)/2. - t1 );
            area -= y1 * (xub-xlb);
            return area / (dx * dy);
        }
        else
            return 1.;
    }
    else
        return 0;
}


/*
 * Initialize with circular well.
 */

void InitCircleExact(SimParams *sp, Array2d & arr, double r, double val) {
    arr = 0;
    for (int i = 0; i < sp->N1; i++) {
        for (int j = 0; j < sp->N2; j++) {
            arr(i,j) = val *ExactVolFracCircle( sp->x_lb + i * sp->dx,sp->y_lb + j * sp->dy,r,sp->dx, sp->dy) ;
        }
    }
}


/*
 * Initialize a quarter of a circle for the diagonal solution.
 */

void InitQuarterCircleDiagonal(SimParams *sp, double * s, double r, double val) {
    Array2d qI(sp->N1 + 2 * sp->c_bounds,sp->N2 + 2 * sp->c_bounds);
    Idx2d idx;
    idx[0] = -sp->c_bounds;
    idx[1] = -sp->c_bounds;
    qI.reindexSelf(idx);
    InitCircleExact(sp, qI, r, val);

    // DIAG
    Array2d q = qI.copy();
    qI.reverseSelf(0);
    qI.reverseSelf(1);
    q = q - qI;

    for (int i = 0; i < sp->N1; i++) {
        for (int j = 0; j < sp->N2; j++)
            s[i*sp->N2 + j] = 	q(i,j);
    }
}

/*
 * Initialize a quarter of a circle for the parallel solution.
 */

void InitQuarterCircleParallel(SimParams *sp, double * s, double r, double val) {
    Array2d qI(sp->N1 + 2 * sp->c_bounds,sp->N2 + 2 * sp->c_bounds);
    Idx2d idx;
    idx[0] = -sp->c_bounds;
    idx[1] = -sp->c_bounds;
    qI.reindexSelf(idx);

    InitCircleExact(sp, qI, r, val);

    // DIAG
    Array2d q = qI.copy();
    qI.reverseSelf(0);
    qI.reverseSelf(1);
    q = q + qI;

    // PARALLEL
    qI.reverseSelf(0);
    q = q - qI;

    qI.reverseSelf(0);
    qI.reverseSelf(1);
    q = q - qI;

    for (int i = 0; i < sp->N1; i++) {
        for (int j = 0; j < sp->N2; j++)
            s[i*sp->N2 + j] = 	q(i,j);
    }
}

/*
 * Use point wells for parallel grid.
 */

void InitOneCellParallel(SimParams *sp, double * q,  double val) {
    for (int i = 0; i< sp->M; i++)		q[i] = 0;
    // injection wells
    q[0]       =	val;						// corresponds to p_00 				(bottom-right corner, origin)
    q[sp->M-1] = 	val;						// corresponds to p_N1-1,N2-1 	(top-right corner)
    // production wells
    q[sp->N2-1]          = -val;				// corresponds to p_0N2-1 			(top-left corner)
    q[(sp->N1-1)*sp->N2] = -val;	            // corresponds to p_N1-1,0 			(bottom-right corner)
}


/*
 * Use point wells for diagonal grid.
 */

void InitOneCellDiagonal	(SimParams *sp, double * q,  double val) {
    for (int i = 0; i< sp->M; i++)		q[i] = 0;
    // injection well
    q[0] = val;									// corresponds to p_00 				(bottom-right corner, origin)
    // production well
    q[sp->M-1] = -val;						    // corresponds to p_N1-1,N2-1 	(top-right corner)
}


/*
 * Calculate the average velocity.
 */

double CalculateAverageVelocity(const Array3d & v, SimParams*sp) {
    double u_avg = 0;
    for (int i=0;  i< sp->N1; i++) {
        for (int j=0;  j< sp->N2; j++) {
            double u1 = 0.5 * (v(i,j,0) + v(i+1,j,0));
            double u2 = 0.5 * (v(i,j,1) + v(i,j+1,1));
            double uij_abs = sqrt(u1*u1 + u2*u2);
            u_avg += uij_abs;
        }
    }
    return u_avg / (sp->N1 * sp->N2);
}


void InitHorizontalFront(SimParams *sp, double *q	, double t) {

    int qlayer_width = round(sp->well_width * (sp->x_ub - sp->x_lb) / sp->dx);
    cout << "layer width: " << qlayer_width << endl;

    if (qlayer_width < 1) {
        cerr << "Injection layer cannot be represented on this grid: (qlayer_width < 1)" << endl;
        exit(1);
    }
    double val = sp->I/(sp->N2 * sp->dx * sp->dy * qlayer_width);
    for (int i = 0; i< sp->M; i++)
        q[i] = 0;

    if (sp->im == PULSING) {
        double tau = 1./sp->q_freq;
        double x = fmod(t, tau);
        if (x > tau/2.) {
            val *= sp->q_high_multiple;
        }
    }

    else if (sp->im == LINEAR) {
        double vl = val;
        double vh = val * sp->q_high_multiple;
        val = vl + (t-sp->tstart) * (vh-vl) / (sp->tfinal-sp->tstart);
        cout << "q( " << t << ") = " << val << endl;
    }

    for  (int i=0; i<sp->N2 * qlayer_width; i++) {
        q[i] = val;
        q[sp->M-1-i] = -val;
    }
}
