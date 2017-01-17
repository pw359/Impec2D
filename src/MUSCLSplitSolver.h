#ifndef MUSCLSPLITSOLVER_H_
#define MUSCLSPLITSOLVER_H_

#include "HyperbolicSplitSolver.h"
#include "GenSimDefs.h"
#include "SimParams.h"

enum SLOPE_LIMITER {SUPERBEE,MINBEE,LEERBEE};

class MUSCLSplitSolver : public HyperbolicSplitSolver {
public:
    MUSCLSplitSolver(SimParams *sp);
    virtual ~MUSCLSplitSolver();
    virtual void  CalculateIntercellFlux_F( Array2d & F,const Array2d & c,	 const Array3d& v,	const double dt);
    virtual void  CalculateIntercellFlux_G( Array2d & G,const Array2d & c,	 const Array3d& v,	const double dt);
private:
    double SlopeLimiter(SLOPE_LIMITER sl, double du_m, double du_p, double w);
    void  Reconstruct(double vL,double vR,double U_i, double U_ip1,double U_im1, double & UL_i, double & UR_i);
    double  GetSlope(double dU_prev, double dU_next);
    double ZetaL(double r, double w) ;
    double ZetaR(double r, double w);

};

#endif /* MUSCLSPLITSOLVER_H_ */



