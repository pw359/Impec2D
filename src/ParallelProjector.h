#ifndef PARALLELPROJECTOR_H_
#define PARALLELPROJECTOR_H_
#include "GenSimDefs.h"
#include <cmath>

class ParallelProjector {
	double angle;
	const Array2d  * p;

	// resolution of the diagonal grid
	double d_N1;
	double d_N2;

	// resolution of the parallel grid
	double p_N1;
	double p_N2;

	// domain bounds: [x_lb, x_ub] X [y_lb, y_ub]
	double x_lb;
	double y_lb;
	double x_ub;
	double y_ub;

	// grid spacing for diagonal grid
	double d_dx;
	double d_dy;

	// grid spacing for parallel grid
	double p_dx;
	double p_dy;
	int bounds;

public:
	ParallelProjector(	double x_lb,double y_lb,	double x_ub,double y_ub, int  d_N1, int d_N2, Array2d const * const parallel, int skipbounds =0);
	virtual ~ParallelProjector();
	double GetFromParallel(double x, double y );
	double GetFromDiagonal(double x, double y );
	void FillDiagonal(Array2d & d, int skip_bounds = 0);
private:
	void DiagonalToParallel(double d_x, double d_y, double & p_x, double & p_y) ;
	void ParallelToDiagonal(double& d_x, double& d_y, double  p_x, double p_y) ;
	double centred(double lb, double dx, int idx);
};



#endif /* PARALLELPROJECTOR_H_ */
