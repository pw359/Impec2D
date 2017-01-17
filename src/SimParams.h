#ifndef SIMPARAMS_H_
#define SIMPARAMS_H_

#include <string>
#include <libconfig.h++>
#include <iostream>
#include <stdlib.h>
#include <vector>

using namespace std;
using namespace libconfig;


enum M_Type 			{ DIRECT, CG, CG_CUDA, MG};
enum HM_Type 			{ UPWINDSPLIT, WAFSPLIT, WAF, CTU, MUSCL};
enum S_Type				{ PARALLEL, DIAGONAL, HORIZONTAL};
enum Solver_Test		{ ELLIPTIC, HYPERBOLIC};
enum Hyp_Test  			{ SINE, TOPHAT, GAUSS, HSINE};
enum Hyp_Norm			{ L1, L2};
enum B_Type 			{ DIRICHLET, SHIFTEDNEUMANN,NEUMANN };
enum DIM_SPLITTING		{ ALTERNATING, X_FIRST, Y_FIRST, STRANG, SYMMETRIC_PARALLEL};
enum AVG_PROC		    { MEANCONC, HARMONIC};
enum TIME_INTEGRATOR 	{ TR,HEUN, FE, HEUNFE};
enum WELL				{ CIRCULAR, GAUSSIAN, ONECELL};
enum P_Type 		    { EXPLICIT, CRANKNICOLSON};
enum SYNCMODE 			{ SYNC, ASYNC};
enum INJ_MODE		    { PULSING, LINEAR, CONSTANT};

class SimParams {
public:

    // Keeps all simulation parameters specified in the input file
    Config cfg;

    // Simulation end time
    double tfinal;

    // Simulation start time
    double tstart;

    int c_bounds;

    // thickness of the boundary layer of the pressure array
    int p_bounds;

    // thickness of the boundary layer of the velocity array; it is always p_bounds-1
    int v_bounds;

    // number of cells in x-direction
    int N1;

    // number of cells in y-direction
    int N2;

    // total number of cells = N1 * N2
    int M;

    // domain bounds: [x_lb, x_ub] X [y_lb, y_ub]
    double x_lb;
    double y_lb;
    double x_ub;
    double y_ub;

    // spatial resolution
    double dx;
    double dy;

    // Method to solve the elliptic equation
    M_Type mtype;

    // Method to solve the elliptic equation
    P_Type ptype;

    // Method to solve the hyperbolic equation
    HM_Type hmtype;

    // source type: parallel or diagonal with respect to the grid;
    S_Type stype;


    // number of grid cells filled in one time step dt
    double I;
    double porosity;

    // true, if logging to file is activated -> creates plots every logskip steps starting with 0
    bool    log_to_file;
    bool    time_matching;
    int     log_skip_steps;
    string  log_dir;

    bool    log_delete_content_first;
    string  results_dir;

    vector<int>     counters;
    vector<double>  times;

    bool    log_mass_conservation;
    bool    check_mass_conservation;
    bool    check_noflow_boundary;
    double  max_mass_frac_error;

    double  inj_fluid_conc;
    double  cfl;
    double  b_val;

    B_Type  btype;
    double  eps;
    double  para_eps;
    double  mu_res;
    double  mu_inj;

    Solver_Test test_component;
    int     testno;
    string  num_result_file;
    string  ex_result_file ;
    string  norm_file ;

    string  ex_u_file;
    string  num_u_file;
    string  ex_v_file;
    string  num_v_file;
    int	    levels, init_x_cells, init_y_cells;
    bool    is_test;

    Hyp_Test    hyp_test;
    Hyp_Norm    hyp_norm;
    bool        gnuplot;
    DIM_SPLITTING dim_splitting;
    bool        check_divergence;
    AVG_PROC    avgproc;
    string      init_file;

    TIME_INTEGRATOR ti;
    double ti_eps;
    double nu;

    bool    diffusion;
    double  alpha_l;
    double  alpha_t;
    bool    levelset;
    int     ls_bounds;
    int     ls_flag_neigh;
    bool    ls_flagged_diff;
    int     ls_iter;

    double  well_rc;
    WELL    well_profile;
    double  q_freq;
    double  q_rand_max;
    int     front_const_width;
    int     front_rand_width;
    double  q_high_multiple;
    INJ_MODE im;
    double well_width;

    double      dt_diff_limit;
    bool        ls_mod_diff;
    double      sigma1, sigma2, sigma0;
    SYNCMODE    sm;

public:
    SimParams(const string & configfile);
    virtual ~SimParams();
};


#define sFORCE_SCHEME 			"force"
#define sFORCE_SPLIT_SCHEME		"force_split"
#define sFDUPWIND_SPLIT_SCHEME 	"fd_upwind_split"
#define sUPWIND_SPLIT_SCHEME 	"upwind_split"
#define sWAF_SPLIT_SCHEME 		"waf_split"
#define sWAF_SCHEME 			"waf"
#define sCTU_SCHEME 			"ctu"
#define sMUSCL_SCHEME		    "muscl"

#define sELL_DIRECT 			"direct"
#define sELL_CG					"cg"
#define sELL_MG					"mg"
#define sELL_CG_CUDA			"cg_cuda"
#define sELL_PARALLEL 		    "parallel"
#define sELL_DIAGONAL 			"diagonal"
#define sELL_HORIZONTAL 		"horizontal"
#define sN1						"x_cells"
#define sN2						"y_cells"
#define sNUM_GHOSTCELLS			"num_ghostcells"
#define sX_LB					"x_lb"
#define sY_LB				    "y_lb"
#define sX_UB					"x_ub"
#define sY_UB					"y_ub"
#define sT_START				"starttime"
#define sT_END					"endtime"
#define sNUM_INJ_CELLS			"I"
#define sLOG_FILE				"log_to_file"
#define sDIR_RESULTS			"results_dir"
#define sCHECK_MASS				"check_mass_conservation"
#define sCHECK_NOFLOW			"check_noflow_boundary"
#define sSOURCE_TYPE			"source_type"

#define sGNUPLOT_ACTIVE			"gnuplot"

#define sLOG_COUNTER_LIST 		"log_counter_list"
#define sLOG_TIME_LIST 			"log_time_list"
#define sHYPER_METHOD			"hyper_method"
#define sELL_METHOD 			"ell_method"

#define sLOG_SKIP_STEPS			"log_skip_steps"
#define sTIME_MATCHING	 	    "time_matching"
#define sLOG_MASS_CONS			"log_mass_conservation"
#define sDIR_LOGGING			"log_dir"
#define sDIM_SPLITTING			"dimensional_splitting"
#define sALTERNATING			"alternating_sweeps"
#define sXFIRST					"x_sweep_first"
#define sYFIRST				    "y_sweep_first"
#define sSTRANGSPLITTING		"strang"
#define sPARALLELSPLITTING		"parallel"

#define sPOROSITY					"porosity"
#define sCFL						"cfl"
#define sINJ_FLUID_CONCENTRATION	"inj_fluid_conc"
#define sELL_BVAL					"const_bound_val"
#define sIS_TEST					"is_test"
#define sHYP_TEST_SINE	    		"sine"
#define sHYP_TEST_HSINE				"horizontal_sine"
#define sHYP_TEST_TOPHAT			"tophat"
#define sHYP_TEST_GAUSS				"gaussbell"
#define sCHECK_DIV					"check_divergence"
#define sAVG_PROC			    	"average_procedure"
#define sHARMONIC					"harmonic"
#define sMEANCONC					"mean_concentration"
#define sNU							"nu"
#define sDIFF					    "diffusion"
#define sDIFF_ALPHA_L				"alpha_l"
#define sDIFF_ALPHA_T				"alpha_t"

#define sLEVELSET					"levelset"
#define sLEVELSET_BC				"ls_num_ghostcells"
#define sLEVELSET_SSIGMA0			"sigma0"
#define sLEVELSET_SSIGMA1			"sigma1"
#define sLEVELSET_SSIGMA2			"sigma2"
#define sLEVELSET_DIFF				"ls_flagged_diff"
#define sLEVELSET_NEIGH				"ls_flag_neigh"
#define sLEVELSET_SDITER			"ls_iter"

#define sWELL_PROFILE				"well_profile"
#define sWELL_RC					"well_rc"

#define sLEVELSET_MODDIFF			"ls_mod_diff"
#endif /* SIMPARAMS_H_ */

#define sPARABOLIC_METHOD           "para_method"
#define sPARABOLIC_EPS              "eps"
#define sSYNCMODE					"syncmode"
#define sASYNC						"async"
#define sSYNC						"sync"

