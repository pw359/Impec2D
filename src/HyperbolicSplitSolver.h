#ifndef HYPERBOLICSPLITSOLVER_H_
#define HYPERBOLICSPLITSOLVER_H_


#include "GenSimDefs.h"
#include "SimParams.h"
#include "HyperbolicSolver.h"


class HyperbolicSplitSolver : public HyperbolicSolver {

protected:
    bool AdvectSplit(Array2d & c, const Array3d& v, double const * const q, const double dt);
    virtual bool Solve(Array2d & c, const Array3d& v, double const * const q, const double dt);
    bool BoundaryFluxCheckF(const Array2d & F);
    bool BoundaryFluxCheckG(const Array2d & G);

public:
    virtual void  CalculateIntercellFlux_F(Array2d & F,const Array2d & c, const Array3d& v,	const double dt) = 0;
    virtual void  CalculateIntercellFlux_G(Array2d & G,const Array2d & c, const Array3d& v,	const double dt) = 0;
    bool Advect_X(Array2d & c, const Array3d& v, const double dt);
    bool Advect_Y(Array2d & c, const Array3d& v, const double dt);
    HyperbolicSplitSolver(SimParams * sp);
    virtual ~HyperbolicSplitSolver();
};

#endif /* HYPERBOLICSPLITSOLVER_H_ */
