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
/* Printer
 *
 * This class provides the interface for file I/O in various formats.
 */
//============================================================================


#include "Printer.h"
#include <iostream>
#include <string>
using namespace std;
#include <fstream>


Printer::Printer() {}
Printer::~Printer() {}

/*
 * Output double array to console.
 */
void Printer::PrintDToConsole(const Array2d & a, int skip_bounds) {
    const int startx = a.base()[0]+skip_bounds;
    const int starty = a.base()[1]+skip_bounds;
    const int endx = startx +a.extent()[0] - 2*skip_bounds;
    const int endy = starty + a.extent()[1] - 2*skip_bounds;
    for (int i = startx; i < endx; i++) {
        for (int j = starty; j < endy; j++)
            cout << "\t" << a(i, j);
        cout << endl;
    }
}

/*
 * Write velocity array to file.
 */

void Printer::PrintVelToConsole( Array3d & a, int skip_bounds ) {

#if USEBOOST
    const int startx = a.index_bases(,0)+skip_bounds;
    const int starty = a.index_bases()[1]+skip_bounds;
    const int endx = startx + a.shape(,0) - 2*skip_bounds;
    const int endy = starty + a.shape()[1] - 2*skip_bounds;
#else
    const int startx = a.base()[0]+skip_bounds;
    const int starty = a.base()[1]+skip_bounds;
    const int endx = startx +a.extent()[0] - 2*skip_bounds;
    const int endy = starty + a.extent()[1] - 2*skip_bounds;
#endif
    for (int l = 0; l < a.shape()[2]; l++) {
        cout << "Component: " << l << endl;
        for (int i = startx; i < endx; i++) {
            for (int j = starty; j < endy; j++)
                cout << "\t" << a(i, j, l);
            cout << endl;
        }
    }
}

/*
 * Write the two velocity components to individual files.
 */

void Printer::PrintVelToPlotFile( std::string & ufile, std::string & vfile,  Array3d & a,
                                  double x_lb, double y_lb, double dx, double dy,  int skip_bounds) {

    ofstream u_file(ufile.c_str());
    ofstream v_file(vfile.c_str());

    if (u_file.is_open() && v_file.is_open()) {
        string SEP = " ";

        const int startx = a.base()[0]+skip_bounds;
        const int starty = a.base()[1]+skip_bounds;
        const int endx = startx +a.extent()[0] - 2*skip_bounds;
        const int endy = starty + a.extent()[1] - 2*skip_bounds;

        for (int i = startx; i < endx; i++) {
            for (int j = starty; j < endy-1; j++) {
                double x_u = x_lb + i* dx;
                double y_u = y_lb + (j+0.5) * dy;
                u_file << x_u << SEP << y_u << SEP << a(i, j,0);
                u_file << endl;
            }
            u_file << endl;
        }
        for (int i = startx; i < endx-1; i++) {
            for (int j = starty; j < endy; j++) {
                double x_v = x_lb + (i+0.5) * dx;
                double y_v = y_lb + j * dy;
                v_file << x_v << SEP << y_v << SEP << a(i, j, 1);
                v_file << endl;
            }
            v_file << endl;
        }
    } else
        cerr << "PRINTER-ERROR: Unable to write to file: " << ufile<< " or " << v_file<< endl;
    u_file.close();
    v_file.close();

}


/*
 * Write double array to a file.
 */

void Printer::PrintDToFile(string & filename, Array2d & a, int skip_bounds) {
    ofstream file(filename.c_str());
    if (file.is_open()) {

        const int startx = a.base()[0]+skip_bounds;
        const int starty = a.base()[1]+skip_bounds;
        const int endx = startx +a.extent()[0] - 2*skip_bounds;
        const int endy = starty + a.extent()[1] - 2*skip_bounds;

        for (int i = startx; i < endx; i++) {
            for (int j = starty; j < endy; j++)
                file << "\t" << a(i, j);
            file << endl;
        }
    } else
        cerr << "PRINTER-ERROR: Unable to write to file: " << filename << endl;
    file.close();
}


/*
 * Write double array to a file.
 */

void Printer::PrintDToPlotFile(string & filename, Array2d & a, double x_lb, double y_lb,
                               double dx, double dy, double skip_bounds) {
    ofstream file(filename.c_str());
    if (file.is_open()) {
        string SEP = " ";

        const int startx = a.base()[0]+skip_bounds;
        const int starty = a.base()[1]+skip_bounds;
        const int endx = startx +a.extent()[0] - 2*skip_bounds;
        const int endy = starty + a.extent()[1] - 2*skip_bounds;
        for (int i = startx; i < endx; i++) {
            for (int j = starty; j < endy; j++) {
                double x = x_lb + (i + 0.5) * dx;
                double y = y_lb + (j+ 0.5) * dy;
                file << x << SEP << y << SEP << a(i, j);
                file << endl;
            }
            file << endl;
        }
    } else
        cerr << "PRINTER-ERROR: Unable to write to file: " << filename << endl;
    file.close();
}
