//=============================================================================
//  Impec2D: Implicit Pressure Explicit Concentration in two spatial dimensions
//  
//  This software solves the incompressible equations for miscible displacement 
//  for two symmetric configurations using a finite-volume solver. The under-
//  lying equations exhibit a physical instability. Small numerical errors,
//  e.g. arising from the particular choice of grid, can therefore lead to 
//  completely different solutions. Please have a look at my dissertation 
//  (goe.pdf) for more information.
//
//  Copyright:  
//
//  The software Impec2D was developed as part of my MPhil course on Scientific 
//  Computing in the Cavendish Laboratory, University of Cambridge. 
//  Copyright 2012, 2016 Peter Wirnsberger (peter.wirnsberger@gmail.com).
//
//  This file is part of Impec2D.
//
//  Impec2D is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  Impec2D is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with Impec2D.  If not, see <http://www.gnu.org/licenses/>.
//============================================================================


//============================================================================
/* SimParams.cpp
 *
 * This class keeps all the input parameters.
 */
//============================================================================



#include "SimParams.h"

SimParams::SimParams(const string & configfile) {

    // read the file
    try {
        cfg.readFile(configfile.c_str());

        // grid and boundary information
        const Setting& root = cfg.getRoot();
        const Setting &genset = root["application"]["general"];
        const Setting &ellipticset = root["application"]["elliptic"];
        const Setting &hyperbolicset = root["application"]["hyperbolic"];
        string sstype, ssync;
        if (!(genset.lookupValue(sN1, N1) && genset.lookupValue(sN2, N2)
                && genset.lookupValue(sNUM_GHOSTCELLS, p_bounds)
                && genset.lookupValue(sNUM_GHOSTCELLS, c_bounds)
                && genset.lookupValue(sX_LB, x_lb)
                && genset.lookupValue(sX_UB, x_ub)
                && genset.lookupValue(sY_UB, y_ub)
                && genset.lookupValue(sY_LB, y_lb)
                && genset.lookupValue(sT_START, tstart)
                && genset.lookupValue(sT_END, tfinal)
                && genset.lookupValue(sNUM_INJ_CELLS, I)
                && genset.lookupValue(sLOG_FILE, log_to_file)
                && genset.lookupValue(sDIR_RESULTS, results_dir)
                && genset.lookupValue(sCHECK_MASS, check_mass_conservation)
                && genset.lookupValue(sCHECK_NOFLOW, check_noflow_boundary)
                && genset.lookupValue(sIS_TEST, is_test)
                && genset.lookupValue(sSOURCE_TYPE, sstype)
                && genset.lookupValue(sGNUPLOT_ACTIVE, gnuplot)
                && genset.lookupValue(sCHECK_DIV, check_divergence)
                && genset.lookupValue(sLEVELSET, levelset)
                && genset.lookupValue(sSYNCMODE, ssync)

             )
           ) {

            std::cerr
                    << "Could not find one or more items in section <<general>>!"
                    << std::endl;
            exit(1);
        }

        if (ssync == sSYNC)
            sm = SYNC;
        else if (ssync == sASYNC)
            sm = ASYNC;
        else	{
            std::cerr << "Config::Synchronization: protocoll is not implemented!" << std::endl;
            exit(1);
        }

        if (check_mass_conservation) {
            if (!(genset.lookupValue("max_mass_frac_error", max_mass_frac_error))) {
                std::cerr
                        << "Could not find one or more options in section <<general>>! You need to specify <<max_mass_frac_error>>"
                        << std::endl;
                exit(1);
            }
        }

        // if <<log_file>> is true, then results will be written to << log_dir >>/<< quantity >>_counter.dat every log_skip_steps
        if (log_to_file) {
            if (!(genset.lookupValue(sLOG_SKIP_STEPS, log_skip_steps)
                    && genset.lookupValue(sTIME_MATCHING, time_matching)
                    && genset.lookupValue(sLOG_MASS_CONS, log_mass_conservation)
                    && genset.lookupValue(sDIR_LOGGING, log_dir))) {

                std::cerr
                        << "Could not find one or more logging-options in section <<general>>!"
                        << std::endl;
                exit(1);
            }
        }

        // source type
        if (sstype == sELL_PARALLEL)
            stype = PARALLEL;
        else if (sstype == sELL_DIAGONAL)
            stype = DIAGONAL;
        else if (sstype == sELL_HORIZONTAL)
            stype = HORIZONTAL;
        else {
            std::cerr << "Config::General: source type is not valid!"
                      << std::endl;
            exit(1);
        }

        // elliptic method
        string smethod, sbtype, savgproc;
        if (!(ellipticset.lookupValue(sELL_METHOD, smethod)
                && ellipticset.lookupValue("bound_type", sbtype)
                && ellipticset.lookupValue("eps", eps)
                && ellipticset.lookupValue("mu_res", mu_res)
                && ellipticset.lookupValue("mu_inj", mu_inj)
                && ellipticset.lookupValue(sELL_BVAL, b_val)
                && ellipticset.lookupValue(sAVG_PROC, savgproc)
             )) {
            std::cerr
                    << "Could not find one or more items in <<elliptic>> section!"
                    << std::endl;
            exit(1);
        }

        // boundary type
        if (sbtype == "neumann")
            btype = NEUMANN;
        else {
            std::cerr << "Config::Elliptic: Boundary type is not supported at the moment!" << std::endl;
            exit(1);
        }

        if (smethod == sELL_DIRECT)
            mtype = DIRECT;
        else if (smethod == sELL_CG)
            mtype = CG;
        else if (smethod == sELL_MG)
            mtype = MG;
        else if (smethod == sELL_CG_CUDA)
            mtype = CG_CUDA;
        else {
            std::cerr << "Config::Elliptic: Method is unknown!" << std::endl;
            exit(1);
        }

        if (savgproc == sHARMONIC)
            avgproc = HARMONIC;
        else if (savgproc == sMEANCONC)
            avgproc = MEANCONC;
        else {
            std::cerr << "Config::Elliptic: Averaging procedure is unknown!" << std::endl;
            exit(1);
        }

        // hyperbolic method
        string shmethod, sdimsplitting;
        if (!(hyperbolicset.lookupValue(sHYPER_METHOD, shmethod)
                && hyperbolicset.lookupValue(sPOROSITY, porosity)
                && hyperbolicset.lookupValue(sCFL, cfl)
                && hyperbolicset.lookupValue(sINJ_FLUID_CONCENTRATION,
                                             inj_fluid_conc)
                && hyperbolicset.lookupValue(sINJ_FLUID_CONCENTRATION,
                                             inj_fluid_conc)
                && hyperbolicset.lookupValue(sNU,	nu)
                && hyperbolicset.lookupValue(sDIFF,	diffusion)
                && hyperbolicset.lookupValue(sDIFF_ALPHA_L,	alpha_l)
                && hyperbolicset.lookupValue(sDIFF_ALPHA_T, alpha_t)
                && hyperbolicset.lookupValue("dt_diff_limiter", dt_diff_limit)
                && hyperbolicset.lookupValue(sDIM_SPLITTING, sdimsplitting)
             )) {
            std::cerr
                    << "Could not find one or more items in <<hyperbolic>> section!"
                    << std::endl;
            exit(1);
        }

        if (shmethod == sWAF_SPLIT_SCHEME) {
            hmtype = WAFSPLIT;
        }
        else if (shmethod == sWAF_SCHEME) {
            hmtype = WAF;
        }
        else if (shmethod == sCTU_SCHEME) {
            hmtype = CTU;
        }
        else if (shmethod == sMUSCL_SCHEME) {
            hmtype = MUSCL;
        }
        else if (shmethod == sUPWIND_SPLIT_SCHEME) {
            hmtype = UPWINDSPLIT;
        }
        else {
            std::cerr << "Config::Hyperbolic: Method is unknown!" << std::endl;
            exit(1);
        }

        // options for dimensional splitting
        if (sdimsplitting == sALTERNATING) {
            dim_splitting = ALTERNATING;
        }
        else if (sdimsplitting == sXFIRST) {
            dim_splitting = X_FIRST;
        }
        else if (sdimsplitting == sYFIRST) {
            dim_splitting = Y_FIRST;
        }
        else if (sdimsplitting == sSTRANGSPLITTING) {
            dim_splitting = STRANG;
        }
        else if (sdimsplitting == sPARALLELSPLITTING) {
            dim_splitting = SYMMETRIC_PARALLEL;
        }

        else {
            std::cerr << "Config::Hyperbolic: Dimenstional Splitting specification is unknown!" << std::endl;
            exit(1);
        }

        // counter and timestep matching
        const Setting &counters_list = genset[sLOG_COUNTER_LIST];
        for (int i = 0; i < counters_list.getLength(); ++i)
            counters.push_back(counters_list[i]);

        const Setting &times_list = genset[sLOG_TIME_LIST];
        for (int i = 0; i < times_list.getLength(); ++i)
            times.push_back(times_list[i]);


        // test
        if (is_test) {
            const Setting& root = cfg.getRoot();
            const Setting &testset = root["application"]["test"];

            string comp_to_test, sprofile;
            int ihypnorm;
            if (!(testset.lookupValue("component", comp_to_test)
                    && testset.lookupValue("testno", testno)
                    && testset.lookupValue("num_result_file", num_result_file)
                    && testset.lookupValue("ex_result_file", ex_result_file)
                    && testset.lookupValue("norm_file", norm_file)
                    && testset.lookupValue("initial_x_cells", init_x_cells)
                    && testset.lookupValue("initial_y_cells", init_y_cells)
                    && testset.lookupValue("num_u_file", num_u_file)
                    && testset.lookupValue("ex_u_file", ex_u_file)
                    && testset.lookupValue("num_v_file", num_v_file)
                    && testset.lookupValue("ex_v_file", ex_v_file)
                    && testset.lookupValue("levels", levels)
                    && testset.lookupValue("profile", sprofile)
                    && testset.lookupValue("hyp_norm", ihypnorm))) {
                std::cerr
                        << "Could not find one or more items in section <<general>>!"
                        << std::endl;
            }

            if (comp_to_test == "elliptic")
                test_component = ELLIPTIC;
            else if (comp_to_test == "hyperbolic")
                test_component = HYPERBOLIC;
            else {
                std::cerr
                        << "TestComponent> Config::General: source type is not valid!"
                        << std::endl;
                exit(1);
            }


            if (sprofile == sHYP_TEST_SINE)
                hyp_test = SINE;
            else if (sprofile ==sHYP_TEST_TOPHAT)
                hyp_test = TOPHAT;
            else if (sprofile ==sHYP_TEST_GAUSS)
                hyp_test = GAUSS;
            else if (sprofile ==sHYP_TEST_HSINE)
                hyp_test = HSINE;


            else {
                std::cerr
                        << "TestComponent> Config:: Test profile is not valid!"
                        << std::endl;
                exit(1);
            }

            if (ihypnorm == 1)
                hyp_norm = L1;
            else if (ihypnorm ==2)
                hyp_norm = L2;
            else {
                std::cerr
                        << "TestComponent> Config:: Norm is not supported!"
                        << std::endl;
                exit(1);
            }
        }

        if (!(genset.lookupValue("init_file", init_file)	)) {
            init_file = "";
        }
        else {
            cout << "Reading initial data from " << init_file << endl;
        }

        string sti;
        if (!(hyperbolicset.lookupValue("time_integration", sti)	&&
                hyperbolicset.lookupValue("ti_eps", ti_eps)	)) {
            std::cerr
                    << "Could not find property <<time_integration>> or <<ti_eps>>!"
                    << std::endl;
            exit(1);
        }

        ti = TR;

        // levelet
        /*
        const Setting &levelset = root["application"]["levelset"];
        if (	!(levelset.lookupValue(sLEVELSET_BC, ls_bounds) &&
                  levelset.lookupValue(sLEVELSET_DIFF, ls_flagged_diff) &&
                  levelset.lookupValue(sLEVELSET_MODDIFF, ls_mod_diff) &&
                  levelset.lookupValue(sLEVELSET_NEIGH, ls_flag_neigh) &&
                  levelset.lookupValue(sLEVELSET_SDITER, ls_iter) &&
                  levelset.lookupValue(sLEVELSET_SSIGMA1, sigma1) &&
                  levelset.lookupValue(sLEVELSET_SSIGMA2, sigma2) &&
                  levelset.lookupValue(sLEVELSET_SSIGMA0, sigma0)
              )) {
            std::cerr
                    << "Could not find one or more items in section <<levelset>>!"
                    << std::endl;
            exit(1);
        }
        */

        // wells
        string swp, smode;
        const Setting &wellset = root["application"]["wells"];
        if (!(wellset.lookupValue(sWELL_PROFILE, swp)&&
                wellset.lookupValue(sWELL_RC, well_rc) &&
                wellset.lookupValue("q_freq",q_freq) &&
                wellset.lookupValue("q_rand_max", q_rand_max) &&
                wellset.lookupValue("front_const_width",front_const_width) &&
                wellset.lookupValue("well_width",well_width) &&
                wellset.lookupValue("front_rand_width",front_rand_width) &&
                wellset.lookupValue("q_high_multiple",q_high_multiple) &&
                wellset.lookupValue("mode", smode)

             )) {
            std::cerr
                    << "Could not find one or more items in section <<wells>>!"
                    << std::endl;
            exit(1);
        }


        if (smode == "pulsing")
            im = PULSING;
        else if (smode == "linear")
            im = LINEAR;
        else if (smode == "constant")
            im = CONSTANT;
        else {
            std::cerr
                    << "Unknown injection mode: " << smode
                    << std::endl;
            exit(1);
        }


        if (swp == "circle")
            well_profile = CIRCULAR;
        else if (swp == "gauss")
            well_profile = GAUSSIAN;
        else if (swp == "one_cell")
            well_profile = ONECELL;
        else {
            std::cerr
                    << "Unknown profile <<wells>>!"
                    << std::endl;
            exit(1);
        }

        // parabolic
        string sparamethod;
        const Setting &parabolicset = root["application"]["parabolic"];
        if (!(  parabolicset.lookupValue(sPARABOLIC_METHOD, sparamethod) &&
                parabolicset.lookupValue(sPARABOLIC_EPS, para_eps)
             )) {
            std::cerr
                    << "Could not find one or more items in section <<parabolic>>!"
                    << std::endl;
            exit(1);
        }


        if (sparamethod == "crank-nicolson")
            ptype = CRANKNICOLSON;
        else if (sparamethod == "explicit")
            ptype = EXPLICIT;
        else {
            std::cerr
                    << "Unknown parabolic solver <<" << sparamethod <<">>!"
                    << std::endl;
            exit(1);
        }

        // Crank-Nicolson requires a multigrid-solver
        if (ptype == CRANKNICOLSON && mtype != MG) {
            std::cerr << "SimParams::InvalidOptions> The Crank-Nicolson scheme requires a Multigrid-Solver. Change ell_method to MG!" << endl;
            exit(1);
        }

    } catch (const FileIOException &fioex) {
        std::cerr << "I/O error while reading file." << std::endl;
        exit(1);
    }
    catch (const ParseException &pex) {
        std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
                  << " - " << pex.getError() << std::endl;
        exit(1);
    }

    // velocity boundaries
    v_bounds = p_bounds;


    // calculate resolution
    dx = (x_ub - x_lb) / N1;
    dy = (y_ub - y_lb) / N2;

    // total number of cells
    M = N1 * N2;

    // sigmas are in units of dx
    sigma0 *= dx;
    sigma1 *= dx;
    sigma2 *= dx;

    if (stype == PARALLEL && (dx != dy  || N1 != N2)) {
        cout << "ERROR: This resolution is not supported for the parallel configuration. The module ParallelProjector assumes square grid blocks and equal resolution in both dimensions." << endl;
        exit(1);
    }
    if ((ls_mod_diff && ! levelset) || (ls_flagged_diff && ! levelset) ) {
        cout<< "ERROR: Cannot have flagged or modified diffusion without the levelset option!" << endl;
        exit(1);
    }
}

SimParams::~SimParams() {}

