#ifndef MGSOLVER_H_
#define MGSOLVER_H_
#include "_hypre_utilities.h"
#include "HYPRE_struct_ls.h"
#include "SimParams.h"
#include "mpi.h"
#include "GenSimDefs.h"
#include "EllipticSolver.h"


class MGSolver : public EllipticSolver {
    int myid, num_procs;		// process id and number of processes
    int nx,ny, Np, pi, pj;
    int ilower[2], iupper[2];
    int solver_id;
    int n_pre, n_post;			// number of pre-smoothing and post-smoothing iterations


    int **ell_dir;              // contains (i,j) indices of the cells for which Dirichlet Boundary conditions will be employed. (EllipticSolver)
    int **para_dir;             // same for the ParabolicSolver
    int N_ell_dir;              // number of cells with Dirichlet boundaries (EllipticSolver)
    int N_para_dir;             // ParabolicSolver

    int num_iterations;			// store number of iterations
    int max_iterations;			// maximum number of iterations
    double final_res_norm;


    HYPRE_StructGrid     grid;
    HYPRE_StructStencil  stencil;
    HYPRE_StructMatrix   A;
    HYPRE_StructVector   b;
    HYPRE_StructVector   x;
    HYPRE_StructSolver   solver;
    HYPRE_StructSolver   precond;

    int * stencil_indices;
    int nentries;

    double mu_res, mob1over4;
    double nu;                   // Weighting factor for the two stencils

public:
    MGSolver(SimParams *sp);
    virtual ~MGSolver();

    void ActivateProcessesForParabolic(const Array2d&c_c, const Array3d &v, double dt);
    void SetupGrid();
    void SetupStencil();
    void SetMatrixCoefficients_Elliptic(const Array2d & c);
    void SetMatrixCoefficients_Parabolic(const Array3d & v, Array2d const * const ls, const double dt);
    void SetRHS_Elliptic(double const * const q);
    void SetRHS_Parabolic(const Array2d &c, const Array3d &v, const double dt, Array2d const * const ls);
    void AssignRHS_Parabolic(double const * const q);
    void SetupDiffusionTensorF(double dx, double dy, const Array3d &v, int i, int j,	double &d11_im12j, double &d12_im12j);
    void SetupDiffusionTensorG(double dx, double dy, const Array3d &v, int i, int j,	double &d22_ijm12, double &d21_ijm12);
    void InitializeComponents();
    void DestroyComponents();
    void SolveSystem(double eps);

    void MapHYPREVectorToVertex(Array2d &p, int p_bounds);
    bool Solve_Elliptic_cc_to_c(const Array2d & c_c, int c_bounds, Array2d &p_c, int p_bounds, double const * const q_c);
    bool Solve_Parabolic_cc_to_c(Array2d & c_c, int c_bounds, const Array3d &v, int v_bounds,
                                 Array2d const*const ls, int ls_bounds, const double dt);

    virtual bool SolveLinearMultigrid(Array2d & c_c, int c_bounds, Array2d &p_c, int p_bounds, double const * const q_v);
    void GetCoefficients_Parabolic(int i, int j, double&W, double &NW, double & N, double & NE, double & E,double & SE, double &S, double &SW, double &C,const Array3d &v,
                                   const double dt,  Array2d const * const ls);
    void GetWeightedMatrixCoefficientsForCell(double nu, int i, int j, double&W, double &NW, double & N, double & NE, double & E,
            double & SE, double &S, double &SW, double &C,		const Array2d &c);
    void GetParallelCoefficientsForCell(int i, int j, double&W, double& S, double&C, double& N, double& E,	const Array2d &c);
    void GetDiagonalCoefficientsForCell(int i, int j, double&W, double &NW, double & N, double & NE, double & E,
                                        double & SE, double &S, double &SW, double &C,		const Array2d &c);

    void GetTransformedMobilityTensor(int i,int j, double& M11_ip12j,double&	M11_im12j,double& M22_ijp12,double& M22_ijm12,
                                      double& M12_ip12j,double& M21_ijp12,double& M12_im12j,double& M21_ijm12,const Array2d &c);


    void PrintHYPREVectorToFile(string filename, HYPRE_StructVector & v);
};

#endif /* MGSOLVER_H_ */
