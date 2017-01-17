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
/*
 * GnuPlotter
 *
 * This class allows for output visualisation with gnuplot by using the C++ 
 * version of the interface originally developed by N. Devillard. Please 
 * have a look at the file gnuplot_i.hpp and the link below to find out
 * more about all author contributions and detailed information on how 
 * the communication works (last accessed on 17/1/2017).
 *
 * http://ndevilla.free.fr/gnuplot/  
 *
 */
//============================================================================


#include "GnuPlotter.h"

//Gnuplot class handles POSIX-Pipe-communication with Gnuplot
#include "gnuplot_i.hpp"
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
using namespace std;

GnuPlotter::GnuPlotter(const string temp_dir) {
    this->temp_dir = temp_dir;
    p = new Gnuplot("");
}

GnuPlotter::GnuPlotter(bool levelset,const string temp_dir) {
    this->temp_dir = temp_dir;
    p = new Gnuplot("");
    p->cmd("set term x11 size 700,700");
    cmds = cmds + "\n" + "reset";
    cmds = cmds + "\n" + "set xtics 0,1";
    cmds = cmds + "\n" + "set ytics 0,1";
    cmds = cmds + "\n" + "set pm3d map";
    cmds = cmds + "\n" + "set title 'c(x,y) Contour Plot'";
    cmds = cmds + "\n" + "splot FILE";
}

GnuPlotter::~GnuPlotter() {
    delete p;
}


void GnuPlotter::ContourPlot(const string & data) {
    string plotcmd = cmds;
    string token = "FILE";
    string datafile = "\"" + data + "\"";
    int found = plotcmd.find(token);
    while (found != string::npos) {
        plotcmd.replace(found,token.length(), datafile);
        found = plotcmd.find(token);
    }
    p->cmd(plotcmd.c_str());
}

void GnuPlotter::ReplaceTokenByFile(string & s, const string & token, const string & file) {
    string datafile = "\"" + file + "\"";
    int found = s.find(token);
    while (found != string::npos) {
        s.replace(found,token.length(), datafile);
        found = s.find(token);
    }
}

void GnuPlotter::PlotConvergence(const string& data) {
    string ccmds;
    ccmds = ccmds + "\n" + "reset";
    ccmds = ccmds + "\n" + "set xlabel 'log(dx)'";
    ccmds = ccmds + "\n" + "set ylabel 'log(err)'";
    ccmds = ccmds + "\n" + "set size 0.7,0.7";
    ccmds = ccmds + "\n" + "f(x)=a*x+b";
    ccmds = ccmds + "\n" + "fit f(x) FILE via a,b";
    ccmds = ccmds + "\n" + "ctitle=sprintf(\"fit, f(x)=%f*x+%f\",a,b)";
    ccmds = ccmds + "\n" + "plot FILE, f(x) t ctitle w lines";

    string token = "FILE";
    string datafile = "\"" + data + "\"";

    int found = ccmds.find(token);
    while (found != string::npos) {
       ccmds.replace(found,token.length(), datafile);
       found = ccmds.find(token);
    }
    p->cmd(ccmds.c_str());



//    ifstream myfile ("conv.plt");
/*    if (myfile.is_open())   {
        while ( myfile.good() )   {
            getline (myfile,line);
            ccmds = ccmds + "\n" + line;
        }
        myfile.close();
        string token = "FILE";
        string datafile = "\"" + data + "\"";

        int found = ccmds.find(token);
        while (found != string::npos) {
            ccmds.replace(found,token.length(), datafile);
            found = ccmds.find(token);
        }
        p->cmd(ccmds.c_str());
    }
    else
        cout << "Unable to open file";
*/
}

