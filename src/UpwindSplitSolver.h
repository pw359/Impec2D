#ifndef UPWINDSPLITSOLVER_H_
#define UPWINDSPLITSOLVER_H_

#include "HyperbolicSplitSolver.h"
#include "GenSimDefs.h"
#include "SimParams.h"

class UpwindSplitSolver : public HyperbolicSplitSolver {

public:
    UpwindSplitSolver(SimParams *sp);
    virtual ~UpwindSplitSolver();
    virtual void  CalculateIntercellFlux_F( Array2d & F,const Array2d & c,	 const Array3d& v,	const double dt);
    virtual void  CalculateIntercellFlux_G( Array2d & G,const Array2d & c,	 const Array3d& v,	const double dt);
    double  CalculateIntercellFlux(double SPEED,double UL, double UR);

};

#endif /* UPWINDSPLITSOLVER_H_ */
