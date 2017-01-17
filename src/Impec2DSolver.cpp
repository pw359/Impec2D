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
/* Impec2DSolver
 *
 * This class is the central element of Impec2D.
 * Specific tasks, such as
 * i)   solving the elliptic equation
 * ii)  solving the hyperbolic equation
 * iii) writing files
 * etc, are invoked here.
 */
//============================================================================

#include "Impec2DSolver.h"
#include <libconfig.h++>
#include "MGSolver.h"
#include "MUSCLSplitSolver.h"
#include "ParallelProjector.h"
#include "GenSimDefs.h"
#include <sys/time.h>
#include <algorithm>
#include <fstream>
#include "Util.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include "stdlib.h"
#include "Printer.h"
#include <string>
#include <sstream>


using namespace std;
using namespace libconfig;

Impec2DSolver::Impec2DSolver(const string& configfile, int myid, int num_procs) {

    this->myid = myid;
    this->num_procs = num_procs;

    q = NULL;
    q_vertex = NULL;
    es = NULL;
    hs = NULL;
    gp = NULL;

    // all processes read parameters from configuration file, one at at time  and store them in an object.
    double buffer[1];
    int from = myid -1;
    int to = myid + 1;
    int tag = 123;
    MPI_Status status;

    if (num_procs > 1) {
        if (myid > 0 && myid < num_procs-1) {
            MPI_Recv(buffer, 1, MPI_DOUBLE, from, tag, MPI_COMM_WORLD, &status);
            sp = new SimParams(configfile);
            MPI_Send(buffer, 1, MPI_DOUBLE, to, tag, MPI_COMM_WORLD);
        }
        else if (myid == 0) {
            sp = new SimParams(configfile);
            MPI_Send(buffer, 1, MPI_DOUBLE, to, tag, MPI_COMM_WORLD);
        }
        else  {
            MPI_Recv(buffer, 1, MPI_DOUBLE, from, tag, MPI_COMM_WORLD, &status);
            sp = new SimParams(configfile);
        }
    }
    else
        sp = new SimParams(configfile);

    if (sp->gnuplot )
        gp = new GnuPlotter(sp->levelset, sp->log_dir);
    else
        gp = NULL;

    Initialize();
}

Impec2DSolver::~Impec2DSolver() {
    if (q)			delete [] q;
    if (q_vertex)	delete [] q_vertex;
    if (es)	        delete es;
    if (hs)	        delete hs;
    if (gp)	        delete gp;

    delete sp;
}



void Impec2DSolver::Initialize() {

    // create output directories, if necessary
    // NOTE: directory can only be created at one level; e.g. if cwd is "." then ./test will work but ./test/diagonal1 not
    SetupDirectories(sp);

    es = new MGSolver(sp);
    hs = new MUSCLSplitSolver(sp);
    hs->AssignMGSolver((MGSolver*)es);

    // resize the arrays (p, v and c) according to their boundaries cells specified in cfg
    ResizeArrays();

    // initialize concentration array
    InitializeConcentration();

    // initialize source term q
    InitializeSource();
}



/*
 *  Allocate the arrays taking into account the boundary cell layer specified in cfg
 */

void Impec2DSolver::ResizeArrays() {

    // pressure
    p.resize(sp->N1 + 2 *sp->p_bounds,sp->N2 + 2 * sp->p_bounds);
    Idx2d idx;
    idx[0] = -sp->p_bounds;
    idx[1] = -sp->p_bounds;
    p.reindexSelf(idx);

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

    // Initialize to zero
    p = 0;
    v = 0;
    c = 0;

    // source term
    q = new double[sp->M];

    // vertex based version
    q_vertex = new double[(sp->N1+1)* (sp->N2+1)];
}



/*
 *  Purpose:
 *  Initialize the concentration array according to the configuration file
 */

void Impec2DSolver::InitializeConcentration() {

    // NOTE: assume ordering by columns (dimension check omitted)
    if (sp->init_file.length() >0) {
        string line;
        double x, y, conc;
        int ctr = 0;
        ifstream myfile (sp->init_file.c_str());
        if (myfile.is_open())   {
            while ( myfile.good() )   {
                getline (myfile,line);
                if (line.length()>0) {
                    stringstream ss;
                    ss << line;
                    ss>>x >> y >> conc;
                    int i = ctr / sp->N2;
                    int j = ctr % sp->N2;
                    c(i, j) = conc;
                    ctr++;
                }
            }
            myfile.close();
            if (ctr != sp->N1 * sp->N2) {
                cerr << "Concentration Input FILE:: wrong dimensions!" << endl;
                exit(1);
            }
        }
        else
            cout << "Unable to open file" << sp->init_file << "; Initializing to zero." << endl;
    }
    else {

#ifdef OMPIMPEC
        #pragma omp parallel for
#endif
        for (int i = -sp->c_bounds; i <= sp->N1 +  sp->c_bounds-1; i++) {
            for (int j = -sp->c_bounds; j <= sp->N2 +  sp->c_bounds-1; j++) {
                c(i, j) = 0.;
            }

        }
    }
}


/*
 *  Purpose:
 *  Initialize the source term q according to the specification in the cfg class
 */

void Impec2DSolver::InitializeSource() {

    double Inj_Rate = sp->I;
    double r = sp->well_rc *(sp->x_ub - sp->x_lb);
    double circval = Inj_Rate * 4. / ( M_PI * r*r);
    double qval = Inj_Rate/(sp->dx * sp->dy);

    if (sp->stype == PARALLEL) {
        switch (sp->well_profile) {
        case CIRCULAR:
            InitQuarterCircleParallel(sp, q,r,circval);
            break;
        case GAUSSIAN:
            break;
        case ONECELL:
            InitOneCellParallel(sp, q, qval);
            break;
        }
    }
    else if (sp->stype == DIAGONAL) {
        switch (sp->well_profile) {
        case CIRCULAR:
            InitQuarterCircleDiagonal(sp, q,r,circval);
            break;
        case GAUSSIAN:
            break;
        case ONECELL:
            InitOneCellDiagonal	(sp, q, qval);
            break;
        }
    }

    else if (sp->stype == HORIZONTAL) {
        InitHorizontalFront(sp, q,  0);
    }
    else {
        cerr << "SOURCE-TERM: unknown source term!" << endl;
        exit(1);
    }
}


/*
 *  Purpose:
 *  This method advances the array c in time from tstart to tfinal.
 *  This involves solving an elliptic pressure equation first, and subsequently
 *  updating the convection equation. Both parts are dealt with by subclasses
 *  of the abstract classes EllipticSolver and HyperbolicSolver.
 *
 */


bool Impec2DSolver::Solve() {

    bool state = true;;
    struct timeval start, end;
    long mtime, seconds, useconds;
    gettimeofday(&start, NULL);

    double t =sp->tstart;
    unsigned int tidx;
    double dt;
    int counter = 0;


    if (myid == 0)
        cout << "Final time: " << sp->tfinal << endl;

    // read the desired output times from configuration file
    if (  find (sp->times.begin(), sp->times.end(), sp->tfinal) == sp->times.end())
        sp->times.push_back(sp->tfinal);
    sort (sp->times.begin(), sp->times.begin()+sp->times.size());
    tidx = 0;


    // set format for console output
    cout.setf(ios_base::fixed, ios_base::floatfield);
    cout.precision(9);

    // store diagnostic information
    vector<DiagState> diagstates;
    DiagState ds;


    double L = sp->x_ub -sp->x_lb;
    if (sp->stype == DIAGONAL)
        L *= sqrt(2.);

    if (sp->diffusion && myid == 0)
        cout <<"->Peclet numbers: (L/alpha_x); (Pe_l, Pe_t)=(" << L/(sp->dx *sp->alpha_l) << " ,  " << L/(sp->dx *sp->alpha_t) << ")" << endl;
    else if (myid == 0)
        cout <<"->Peclet numbers: (L/alpha_x); (Pe_l, Pe_t)=(inf, inf)" << endl;


    double maxc = 0;
    double minc = 10e10;
    double avg_fillperdt = 0;

    while (state && t < sp->tfinal) {
        counter++;

        // fill boundary cells for the elliptic step
        hs->FillBoundaries(c);

        /*****************************************************************************
         * 		SOLVE THE GENERAL PRESSURE EQUATION
         *****************************************************************************/
        state = state && es->Solve(p, c,q);

        if (!state) {
            cerr << "ERROR::EllipticSolver: [" << counter << "]\t" << "t:\t" << t << endl;
            break;
        }

        /*****************************************************************************
         * 		CALCULATE THE NEW VELOCITIES
         *****************************************************************************/
        es->CalcVelocities(p,v,c);

        if (sp->diffusion && counter == 2 && myid == 0) {
            double u_avg = CalculateAverageVelocity(v,sp);
            cout << "Dl x <|u|>: " << sp->alpha_l * sp->dx  * u_avg << endl;
            cout << "Dt x <|u|>: " << sp->alpha_t * sp->dx  * u_avg << endl;
            cout <<"->Peclet numbers: (2Q/<|u|>); (Pe_l, Pe_t)=(" << 2.*sp->I/(u_avg * sp->alpha_l * sp->dx) << " ,  " << 2.*sp->I/(u_avg * sp->alpha_t * sp->dx) << ")" << endl;
        }

        /*****************************************************************************
         * 		CALCULATE THE MAXIMUM STABLE TIMESTEP
         *****************************************************************************/
        dt = hs->MaxStableTimestep(c,v);

        /*****************************************************************************
         * 		MATCH OUTPUT TIMES
         *****************************************************************************/
        if (tidx < sp->times.size() &&  t + dt >= sp->times[tidx]) {
            dt = sp->times[tidx] - t;
            tidx++;
        }

        /*****************************************************************************
         * 		SOLVE THE HYPERBOLIC EQUATION (Finite Volume)
         *****************************************************************************/
        state = state && hs->Solve(c,v,q,dt);

        if (!state) {
            cerr << "ERROR::HyperbolicSolver: [" << counter << "]\t" << "t:\t"<< "\tdt:\t" << dt << endl;
            break;
        }

        // Update time
        t += dt;
        ds.ctr 				= counter;
        ds.t				= t;
        ds.dt				= dt;
        ds.its				=	es->GetLastIts();
        ds.fillperdt		= dt * sp->I / (sp->dx * sp->dy);
        avg_fillperdt		+= ds.fillperdt;


        if (sp->check_mass_conservation)
            ds.mass_err		= CalculateMassError(t);
        else
            ds.mass_err		= -1;

        if (sp->check_noflow_boundary)
            ds.vel_err 		= CalculateVelocityBoundaryError() ;
        else
            ds.vel_err		= -1;

        if (sp->check_divergence)
            ds.diff_err 	= CalculateDivergenceError();
        else
            ds.diff_err		=	-1;

        diagstates.push_back(ds);

        if (myid == 0)
            OutputState(ds);

        // log results to file
        if (sp->log_to_file && myid == 0) {
            if (IsLoggingTime(t))
                LogResultsMatch( t, true);

            // log_skip_steps = 0 encodes "log only steps given in the list"
            if ( IsLoggingStep(counter) || (sp->log_skip_steps != 0 &&   (counter == 1 || counter % sp->log_skip_steps == 0)))
                LogResultsMatch( counter, false);
        }

        double cmax = max(c) ;
        double cmin = min(c);
        if (cmax > maxc)
            maxc= cmax;
        if (cmin < minc)
            minc = cmin;

        if (sp->stype == HORIZONTAL && sp->im != CONSTANT) {
            InitHorizontalFront(sp, q,  t);
        }

        /**********************************************/
    }
    avg_fillperdt /= counter;
    if (myid == 0) {
        cout << "overall maximum: " << maxc << endl;
        cout << "overall minimum: " << minc << endl;
        cout << "average filled cells per dt: " << avg_fillperdt << endl;

        cout.unsetf(ios_base::fixed);

        // timing
        gettimeofday(&end, NULL);
        seconds = end.tv_sec - start.tv_sec;
        useconds = end.tv_usec - start.tv_usec;
        mtime = ((seconds) * 1000 + useconds / 1000.0) + 0.5;
        cout << "---------------------------------------------------------" << endl;
        if (counter > 0)
            cout << "Average no. of Iterations: " << es->GetAvgIts() << endl;
        cout << "Elapsed time [ms]: " << mtime << endl;

        LogResultsMatch(-1, false);
        WriteDiagnostics(diagstates);
    }

    return state;
}


void Impec2DSolver::WriteDiagnostics(const vector<DiagState> & diagstates) {
    string fname = sp->results_dir + "/diagnostics.dat";
    ofstream file(fname.c_str());
    if (file.is_open()) {
        string SEP = "    ";
        DiagState ds;
        for (size_t i = 0; i < diagstates.size(); i++) {
            ds = diagstates[i];
            file << ds.ctr		<< SEP;
            file <<	ds.its		<< SEP;
            file <<	ds.t			<< SEP;
            file <<	ds.dt			<< SEP;
            file << ds.fillperdt << SEP;
            if (sp->check_mass_conservation)
                file <<	ds.mass_err		<< SEP;
            if (sp->check_noflow_boundary)
                file <<	ds.vel_err			<< SEP;
            if (sp->check_divergence)
                file <<	ds.diff_err		<< SEP;
            file << endl;
        }
        file.close();
    }
    else
        cerr << "Unable to write to diagnostics-file: " <<fname<< endl;
}




void Impec2DSolver::OutputState(const DiagState & ds) {
    stringstream ss;
    int prec = 6;
    ss.setf(ios_base::fixed, ios_base::floatfield);
    ss.precision(prec);
    ss << "[" << ds.ctr << "]\t" << "ITS_e: " << ds.its;
    ss << "\tITS_p: " << es->GetLastParabolicIts();
    ss << "\t t: " << ds.t << "\tdt: " << ds.dt << "\tcells/dt: " << ds.fillperdt;

    if (sp->check_mass_conservation) {
        ss.precision(prec);
        ss  << "\tE_m: ";
        ss << ds.mass_err*100;
        ss.precision(prec);
    }
    if (sp->check_noflow_boundary) {
        ss.precision(prec);
        ss  << "\tE_v: ";
        ss <<	ds.vel_err;
        ss.precision(prec);
    }
    if (sp->check_divergence) {
        ss.precision(prec);
        ss  << "\tE_div: ";
        ss <<	ds.diff_err;
        ss.precision(prec);
    }
    cout << ss.str() << endl;
}


double Impec2DSolver::CalculateMassError(double t) {
    double mass = GetTotalMass();
    double num_cells_per_vol = 1./(sp->dx * sp->dy);
    double mass_influx = (t-sp->tstart) * sp->I  * num_cells_per_vol     ;

    if (sp->stype == PARALLEL)
        mass_influx *= 2.;

    double mass_outflux = hs->GetTotalMassOutflux();
    return fabs((mass + mass_outflux - mass_influx)/mass_influx);
}

double Impec2DSolver::CalculateDivergenceError() {
    double err = 0.;
    double sum_source = 0.;

#ifdef OMPIMPEC
    #pragma omp parallel for reduction(+:err,sum_source)
#endif
    for (int i=0; i < sp->N1; i++) {
        for (int j=0; j < sp->N2; j++) {
            double div = (v(i, j+1, 1) - v(i, j, 1))/sp->dy +  (v(i+1, j,0) - v(i, j,0))/sp->dx;
            err += fabs(div - q[i * sp->N2 +j]);
            sum_source += fabs(q[i * sp->N2 +j]);
        }
    }
    return err * sp->dx *  sp->dx;
}




double Impec2DSolver::CalculateVelocityBoundaryError() {
    double err_sum = 0.;
    double retval = -1.;
#ifdef OMPIMPEC
    #pragma omp parallel for reduction(+:err_sum)
#endif
    for (int i=0; i < sp->N1; i++) {
        // bottom
        err_sum += fabs((double)v(i,0,1));

        // top
        err_sum += fabs((double)v(i, sp->N2, 1));
    }

#ifdef OMPIMPEC
    #pragma omp parallel for reduction(+:err_sum)
#endif
    for (int j=0; j < sp->N2; j++) {
        //left
        err_sum += fabs((double)v(0,j,0));

        //right
        err_sum += fabs((double)v(sp->N1, j,0));
    }
    retval = err_sum;
    return retval;
}





bool  Impec2DSolver::IsLoggingStep(int counter) {

    // return true if counter is in the list <<coutners>>
    return (find (sp->counters.begin(), sp->counters.end(), counter) != sp->counters.end());
}

bool  Impec2DSolver::IsLoggingTime(double t) {

    // return true if t is in the list <<times>>
    return (find (sp->times.begin(),sp-> times.end(), t) != sp->times.end());
}



double Impec2DSolver::GetTotalMass() {
    double sum_mass = 0;


#ifdef OMPIMPEC
    #pragma omp parallel for reduction(+:sum_mass)
#endif
    for (int i=0; i<sp->N1; i++) {
        for (int j=0; j<sp->N2; j++)
            sum_mass += c(i, j);
    }
    return sum_mass;
}

void Impec2DSolver::LogResultsMatch( double postfix, bool logTime) {
    string dir, spostfix;
    //postfix < 0 marks the end of the simulation
    if (postfix<0) {
        dir = sp->results_dir;
        spostfix = ".dat";
    }
    else {
        dir = sp->log_dir;
        stringstream ss;
        ss << "_";
        if (logTime)
            ss << "t_";
        ss << postfix ;
        spostfix 		=	ss.str() + ".dat";
    }

    string cfile 		= dir + "/" + "concentration" 	+ spostfix;
    string pfile		= dir + "/" + "pressure"		+ spostfix;
    string ufile 		= dir + "/velocity_u"			+ spostfix;
    string vfile 		= dir + "/velocity_v"			+ spostfix;
    Printer pr;


    bool log_p  = (postfix < 0);
    bool log_v = (postfix < 0)	;
    bool log_c = true;

    if (sp->stype == PARALLEL ) {
        int d_N1 , d_N2;
        double d_dx, d_dy ;
        Array2d diag_c, diag_ls, diag_p;

        // concentration
        if (log_c)
            Project(diag_c, &c, sp->c_bounds, d_N1 ,d_N2, d_dx, d_dy);

        if (log_c)
            pr.PrintDToPlotFile(cfile, diag_c,sp->x_lb, sp->y_lb, d_dx, d_dy, sp->c_bounds);

        // pressure
        if (log_p) {
            Project(diag_p, &p, sp->p_bounds, d_N1 ,d_N2, d_dx, d_dy);
            pr.PrintDToPlotFile(pfile, diag_p,sp->x_lb, sp->y_lb, d_dx, d_dy, sp->p_bounds);
        }
    }
    else {
        // concentration
        if (log_c)
            pr.PrintDToPlotFile(cfile, c,sp-> x_lb, sp->y_lb, sp->dx, sp->dy, sp->c_bounds);

        // pressure
        if (log_p) {
            pr.PrintDToPlotFile(pfile, p,sp->x_lb, sp->y_lb, sp->dx, sp->dy, sp->p_bounds);
        }

        if (log_v)
            pr.PrintVelToPlotFile(ufile, vfile, v, sp->x_lb, sp->y_lb,sp->dx, sp->dy, sp->v_bounds);
    }

    // update gnuplot window
    if (sp->gnuplot) {
        gp->ContourPlot(cfile);
    }
}



void Impec2DSolver::Project(Array2d & arrd, Array2d const * const arrp, int bounds, int & d_N1 , int &d_N2, double &d_dx, double &d_dy) {

    // calculate the resolution of the correspondig diagonal grid
    d_N1 = round(sp->N1/sqrt(2));
    d_N2 = round(sp->N2/sqrt(2));


    // pj will cary out the projection/interpolation operation
    ParallelProjector pj(sp->x_lb,sp->y_lb, sp->x_ub, sp->y_ub, d_N1, d_N2,  arrp, bounds);


    // allocate temporary array
    arrd.resize(d_N1 + 2 * bounds,d_N2 + 2 *bounds);
    Idx2d idx;
    idx[0] = -bounds;
    idx[1] = -bounds;
    arrd.reindexSelf(idx);
    arrd = 0;


    // fill all the elements of diag_c with the corresponding interpolated values from the parallel grid
    pj.FillDiagonal(arrd, bounds);

    // write results to file
    d_dx = (sp->x_ub - sp->x_lb)/(sqrt(2.) * d_N1);
    d_dy = (sp->y_ub - sp->y_lb)/(sqrt(2.) *d_N2);
}


bool Impec2DSolver::TestComponent() {

    if (myid == 0)
        cout << "Executing test no. " << sp->testno << " for the ";

    if (sp->test_component == ELLIPTIC) {
        if (myid == 0)
            cout << "ELLIPTIC Solver " << endl;
        return es->TestNeumann(gp);
    }
    else if (sp->test_component == HYPERBOLIC) {
        if (myid == 0)
            cout << "Hyperbolic Solver " << endl;
        switch (sp->testno) {
        case 1:
            return TestAdvection();
            break;
        default:
            cerr << "Hyperbolic TEST:: Unknown test number: " << sp->testno << endl;
            return false;
        }
    }
    else {
        cout << " HYPERBOLIC Solver tests are not yet implemented!" << endl;
        return true;
    }
}

bool Impec2DSolver::Run() {
    bool stat;
    if (sp->is_test)
        stat =  this->TestComponent();
    else
        stat =  this->Solve();
    if (sp->gnuplot && myid == 0) {
        cout << endl << "Press ENTER to continue..." << endl;
        std::cin.clear();
        std::cin.ignore(std::cin.rdbuf()->in_avail());
        std::cin.get();
    }

    return stat;
}



/*
 *  Purpose:
*
*	This method calls "Advection()" several times for different resolutions according to the configuration file.
*	Advection returns the L2-norm which is then written to a file.
 *
 */

bool Impec2DSolver::TestAdvection() {

    bool stat = true;

    // start with the resolution as specified in the test-section in the configuration file
    sp->N1 = sp->init_x_cells;
    sp->N2 = sp->init_y_cells;


    ofstream norm_stream(sp->norm_file.c_str());

    // delete q because it was previously allocated in the constructor of Impec2DSolver
    delete []q;

    norm_stream.precision(30);
    if (norm_stream.is_open()) {
        for (int k = 0; k < sp->levels; k++) {

            // calculate the new resolution
            sp->dx = (sp->x_ub - sp->x_lb)/sp->N1;
            sp->dy = (sp->y_ub - sp->y_lb)/sp->N2;
            sp->M = sp->N1*sp->N2;

            // resize the arrays accordingly
            Idx2d idx;
            idx[0] = -sp->p_bounds;
            idx[1] = -sp->p_bounds;

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

            v = 0;
            c = 0;

            // allocate a new source term and initialize it to zero
            q = (double*) calloc(sp->M, sizeof(double));

            cout << "->Resolution:\t" << sp->N1 << "  x  " << sp->N2 << flush;
            double lpnorm = Advection();
            cout << "\t\tL" << sp->hyp_norm +1 << ":\t" << lpnorm << endl;

            // write norms to file
            norm_stream << log(sp->dx) << "\t" << log(lpnorm) << endl;
            free(q);

            // double the resolution
            sp->N1 *= 2;
            sp->N2 *= 2;
        }
    }
    else {
        cerr << "Cannot write norm-file! Convergence test aborted..." << endl;
        stat = false;
    }
    norm_stream.close();

    if (sp->levels > 1) {
        if (sp->gnuplot) {
            gp->PlotConvergence( sp->norm_file);
        }
    }
    q = NULL;
    return stat;
}

double Impec2DSolver::Advection() {

    Printer pr;
    stringstream ss;
    double (*fptr)(SimParams const*const, double,double) = NULL;
    double (*normptr)(SimParams const * const sp,const Array2d& c, double t, double (*fptr)(SimParams const*const,double, double)) = NULL;

    // i) choose profile according to the configuration file
    switch (sp->hyp_test) {
    case SINE:
        fptr=Sinewave;
        break;
    case TOPHAT:
        fptr=Tophat;
        break;
    case GAUSS:
        fptr=Gaussbell;
        break;
    case HSINE:
        fptr=VerticalSine;
        break;
    }
    if (!fptr) {
        cerr << "Advection:: Uninitialized function pointer!" << endl;
        exit(1);
    }

    // ii) choose profile according to the configuration file
    switch (sp->hyp_norm) {
    case L1:
        normptr=CalcAdvectionL1Norm;
        break;
    case L2:
        normptr=CalcAdvectionL2Norm;
        break;
    }
    if (!normptr) {
        cerr << "Advection:: Uninitialized norm pointer!" << endl;
        exit(1);
    }


    // iii) initialize concentration (without filling the boundaries)
    InitAdvection(fptr);
    ss <<  sp->num_result_file<<  "init" << sp->N1;
    string s = ss.str();
    pr.PrintDToPlotFile(s, c, sp->x_lb, sp->y_lb, sp->dx, sp->dy   , sp->p_bounds );


    // iv) initialize Velocity (const) including boundaries
    CalcAdvectionVelocities(v,c,0.);
    ss.str("");
    ss <<  sp->ex_u_file << sp->N1;
    string ufile = ss.str();
    ss.str("");
    ss <<  sp->ex_v_file << sp->N1;
    string vfile = ss.str();
    pr.PrintVelToPlotFile(ufile, vfile, v, sp->x_lb, sp->y_lb,sp->dx, sp->dy /*, sp->v_bounds*/);


    // v) get maximum stable timestep
    double dt_max =  hs->MaxStableAdvectionTimestep(c,v);

    // vi) time stepping
    bool state = true;;
    double t =sp->tstart;
    unsigned int tidx;
    int counter = 0;

    if (  find (sp->times.begin(), sp->times.end(), sp->tfinal) == sp->times.end())
        sp->times.push_back(sp->tfinal);
    sort (sp->times.begin(), sp->times.begin()+sp->times.size());
    tidx = 0;


    cout.setf(ios_base::fixed, ios_base::floatfield);
    cout.precision(9);

    // get start time
    struct timeval start, end;
    long mtime, seconds, useconds;
    gettimeofday(&start, NULL);
    double dt=dt_max;

    DiagState ds;
    while (state && t < sp->tfinal) {
        counter++;
        // match the speicified times
        if (tidx < sp->times.size() &&  t + dt >= sp->times[tidx]) {
            dt = sp->times[tidx] - t;
            tidx++;
        }
        else
            dt = dt_max;

        // solve the hyperbolic part
        state = state && hs->Solve(c,v,q,dt);
        if (!state) {
            cerr << "ERROR::HyperbolicSolver: [" << counter << "]\t" << "t:\t"<< "\tdt:\t" << dt << endl;
            break;
        }

        t += dt;

        // output the details to the console
        ds.ctr 	= counter;
        ds.dt 	= dt;
        ds.t	= t;
        ds.its  = 0;

        if (counter % 100)
            OutputState(ds);

        // log results to file
        if (sp->log_to_file) {
            if (IsLoggingTime(t))
                LogResultsMatch( t, true);

            // log_skip_steps = 0 encodes "log only steps given in the list"
            if ( IsLoggingStep(counter) || (sp->log_skip_steps != 0 &&   (counter == 1 || counter % sp->log_skip_steps == 0)))  {
                LogResultsMatch( counter, false);
            }

        }

    }
    cout.unsetf(ios_base::fixed);

    // timing
    gettimeofday(&end, NULL);
    seconds = end.tv_sec - start.tv_sec;
    useconds = end.tv_usec - start.tv_usec;
    mtime = ((seconds) * 1000 + useconds / 1000.0) + 0.5;
    cout << "\ttime [ms]: " << mtime << flush;

    ss.str("");
    ss << sp->num_result_file;
    ss << sp->N1 ;
    s = ss.str();
    pr.PrintDToPlotFile(s, c, sp->x_lb, sp->y_lb, sp->dx, sp->dy, sp->p_bounds);
    return normptr(sp,c,t,fptr);

}


void Impec2DSolver::InitAdvection(double (*fptr)(SimParams const*const,double, double)) {
    for (int i=0; i < sp->N1 ; i++) {
        for (int j=0; j < sp->N2 ; j++) {
            double x = sp->x_lb + (i + 0.5) * sp->dx;
            double y = sp->y_lb + (j + 0.5) * sp->dy;
            c(i, j) = fptr(sp,x,y);
        }
    }
}

void Impec2DSolver::CalcAdvectionVelocities(Array3d &v, const Array2d& c, double t) {
    double diag_vec_x = -(sp->x_ub-sp->x_lb);
    double diag_vec_y = (sp->y_ub-sp->y_lb);
    for (int i=-sp->v_bounds; i<=sp->N1 + sp->v_bounds; i++) {
        for (int j=-sp->v_bounds; j<=sp->N2 + sp->v_bounds; j++) {

            // constant velocity at the moment; velocity such that the domain is crossed
            // diagonally once in 1 Time Unit
            v(i, j,0) = diag_vec_x;
            v(i, j, 1) = diag_vec_y;
        }
    }
}
