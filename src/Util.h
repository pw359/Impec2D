#ifndef UTIL_H_
#define UTIL_H_

#include "SimParams.h"
#include "GenSimDefs.h"


double Tophat(SimParams const * const sp, double x, double y);
double Gaussbell(SimParams const * const sp,double x, double y);
double DiagonalSinewave(SimParams const * const sp,double x, double y);
double HorizontalSine(SimParams const * const sp,double x, double y);
double VerticalSine(SimParams const * const sp,double x, double y);
double Sinewave(SimParams const * const sp,double x, double y);
double CalcAdvectionL2Norm(SimParams const * const sp,const Array2d& c, double t, double (*fptr)(SimParams const*const,double, double));
double CalcAdvectionL1Norm(SimParams const * const sp,const Array2d& c, double t, double (*fptr)(SimParams const*const,double, double));
void SetupDirectories(SimParams const * const sp);
void FillGhostCellsShiftedNeumannScalar(SimParams const * const sp,Array2d &arr);
void FillGhostCellsNeumannScalar(Array2d &arr, int bounds);
void FillGhostCellsPeriodicScalar(Array2d &arr, int n);
void FillGhostCellsNeumannVelocity(Array3d &v, int n);
double ExactVolFracCircle(double x1, double y1, double r, double dx, double dy);
void InitCircleExact(SimParams *sp, Array2d & arr, double r, double val) ;
void InitQuarterCircleDiagonal(SimParams *sp, double * s, double r, double val) ;
void InitQuarterCircleParallel(SimParams *sp, double * s, double r, double val) ;
void InitOneCellParallel(SimParams *sp, double * q,  double val);
void InitOneCellDiagonal	(SimParams *sp, double * q,  double val);
double CalculateAverageVelocity(const Array3d & v, SimParams*sp);
void InitHorizontalFront(SimParams *sp, double * s,  double t);

#endif /* UTIL_H_ */
