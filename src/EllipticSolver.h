#ifndef ELLIPTICSOLVER_H_
#define ELLIPTICSOLVER_H_

#include <vector>
#include "GenSimDefs.h"
#include "SimParams.h"
#include "GnuPlotter.h"
using namespace std;
using namespace libconfig;

enum Stencil {
	DIAG, PARA
};

class EllipticSolver {
	double mob1over4;
	bool firstTry;
	int total_its;
	int calls;
	double nu;


protected:
	int last_its;
  int last_parabolic_its;

  double secamt;
  double secamt2;
  double cosa2;
  double sina2;
  double costsina;
  double cosasint;
  double cost2;
  double sint2;
 
  double cosa;
  double sint;
  double sina;
  double cost;

  double alpha;
  double theta;
  double dxi;
  double deta;

  double x_lb;
  double x_ub;
  double y_lb;
  double y_ub;
  double dx;
  double dy;
  double dA;

	SimParams * sp;

public:

	EllipticSolver(SimParams * sp);
	virtual ~EllipticSolver();

	virtual bool SolveLinearMultigrid(Array2d & c_c, int c_bounds, Array2d &p_c, int p_bounds, double const * const q_v) = 0;

	bool Solve(Array2d & p, Array2d &c, double * q);
	void CalcVelocities(Array2d & p, Array3d &v, Array2d& c);
	void SetNu(double nu);
	int GetLastIts();
  int GetLastParabolicIts();
	double GetAvgIts();
	bool TestNeumann(GnuPlotter *gp);

protected:
	void Initialize();


	double K_11(double x, double y);
	double K_22(double x, double y);
	double KoverMu_11_im12j_harmonic(double c_im1j, double c_ij, int i, int j);
	double KoverMu_22_ijm12_harmonic(double c_ijm1, double c_ij, int i, int j);
	double KoverMu_11_im12j_mean_conc(double c_im1j, double c_ij, int i, int j);
	double KoverMu_22_ijm12_mean_conc(double c_ijm1, double c_ij, int i, int j);
	double KoverMu_11_im12j(double c_im1j, double c_ij, int i, int j);
	double KoverMu_22_ijm12(double c_ijm1, double c_ij, int i, int j);

	double mu(double c);

  double KoverMu_11_ip12jp12(int i, int j, double c_ij, double c_ip1j,	double c_ijp1, double c_ipjp1);
	double KoverMu_22_ip12jp12(int i, int j, double c_ij, double c_ip1j,	double c_ijp1, double c_ip1jp1);

	double GetParallelU1_im12j(const Array2d &c, const Array2d &p, int i,		int j);
	double GetParallelU2_ijm12(const Array2d &c, const Array2d &p, int i,		int j);
	double GetDiagonalU1_im12j(const Array2d &c, const Array2d &p, int i,	int j);
	double GetDiagonalU2_ijm12(const Array2d &c, const Array2d &p, int i,	int j);
	void   GetDiagonalW_im12j(double &w1_im12j, double &w2_im12j,	const Array2d &c, const Array2d &p, const int i, const int j);

	void M_im12j(double & M11_im12j, double & M12_im12j, double &M21_im12j,	double &M22_im12j, const Array2d&c, int i, int j);
	void XYToXiEta(const double v1, const double v2, double &w1, double &w2);
	void XYFromXiEta(double &v1, double& v2, const double w1, const double w2);

  double qTest1(double x, double y);
  double qTest2(double x, double y);
  double qTest3(double x, double y);
  double qTest4(double x, double y);
  double pTest (double x, double y);
private:

	void InitializeEllipticNeumannTest(Array2d &c, Array2d &p, Array2d &p_ex,
			double *q, Array3d & v_ex);
	void GetParameters();

};

#endif /* ELLIPTICSOLVER_H_ */
