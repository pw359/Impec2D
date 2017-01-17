#ifndef HYPERBOLICSOLVER_H_
#define HYPERBOLICSOLVER_H_
#include <libconfig.h++>
#include "MGSolver.h"

using namespace libconfig;
#include "GenSimDefs.h"
#include "SimParams.h"
#include "TimeIntegrator.h"


enum SWEEP {X_SWEEP, Y_SWEEP};

class HyperbolicSolver {

private:
    TimeIntegrator * ti;
    MGSolver * mg;
protected:
    SimParams * sp;
    double total_mass_outflux;
    SWEEP prev_first_sweep;
    bool SourceSplit(Array2d & c,  double const * const q, const double dt);
    double Flux (double vel, double conc);

public:
    virtual bool Solve(Array2d & c, const Array3d& v, double const * const q, const double dt) = 0;
    HyperbolicSolver(SimParams *sp);
    virtual ~HyperbolicSolver();
    double GetTotalMassOutflux() {
        return total_mass_outflux;
    }
    void FillBoundaries(Array2d &c);
    double MaxStableTimestepCombined(const Array2d & c, const Array3d& v);
    double MaxStableAdvectionTimestep(const Array2d & c, const Array3d& v);
    double MaxStableTimestep(const Array2d & c,const  Array3d& v);
    bool DiffusionSplit(Array2d &c, const Array3d &v, const double dt);
    bool DiffusionSplitCN(Array2d &c, const Array3d &v, const double dt);
    void AssignMGSolver(MGSolver *mg);
private:
    void OutputState(int, double, double);
    void Initialize();
};

#endif /* HYPERBOLICSOLVER_H_ */
