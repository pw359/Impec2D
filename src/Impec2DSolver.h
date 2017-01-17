#ifndef IMPEC2DSOLVER_H_
#define IMPEC2DSOLVER_H_

#include <libconfig.h++>
#include <string>
#include "EllipticSolver.h"
#include "HyperbolicSolver.h"
#include "GenSimDefs.h"
#include "SimParams.h"
#include "GnuPlotter.h"

using namespace std;
using namespace libconfig;

// diagnostic state (for output and history)
struct DiagState {
    int ctr;
    int its;
    double t;
    double dt;
    double mass_err;
    double diff_err;
    double vel_err;
    double fillperdt;
};

class Impec2DSolver {

    // multidimensional arrays to store pressure (p), concentration (c) and velocity (v)
    // the pressure and the concentration are stored in the cell-centre whereas the velocity is defined on the interface
    Array2d p;
    Array2d c;
    Array3d v;

    // solver for the Elliptic Part (Abstract Class)
    EllipticSolver * es;

    // solver for the Hyperbolic Part (Abstract Class)
    HyperbolicSolver * hs;

    // discrete form of pressure equation:   A * x = -q + b;
    // in the paper of Bell and Shubin!
    double * q;

    // vertex based version of q
    double * q_vertex;

    GnuPlotter * gp;

    int myid;
    int num_procs;

public:
    SimParams * sp;
    Impec2DSolver(const string& configfile, int myid, int num_procs);
    virtual ~Impec2DSolver();
    bool Run();
private:
    double CalculateMassError(double t);
    double CalculateVelocityBoundaryError();
    void OutputState(const DiagState & ds);
    void WriteDiagnostics(const vector<DiagState> & diagstates);
    void InitAdvection(double (*fptr)(SimParams const*const,double, double));
    double CalculateDivergenceError();
    void CalcAdvectionVelocities(Array3d &v, const Array2d& c, double t);
    void InitializeConcentration();
    void Initialize();
    void InitializeSource();
    void Project(Array2d & arrd, Array2d const * const arrp, int bounds, int & d_N1 , int &d_N2, double &d_dx, double &d_dy);
    void ResizeArrays();
    void LogResultsMatch(double postfix, bool logTime);
    bool  IsLoggingStep(int counter);
    bool IsLoggingTime(double t);
    double GetTotalMass();
    bool Solve();
    bool TestComponent();
    bool TestAdvection();
    double Advection();

};

#endif /* IMPEC2DSOLVER_H_ */
