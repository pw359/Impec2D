#ifndef PRINTER_H_
#define PRINTER_H_
#include "GenSimDefs.h"
#include <string>

class Printer {
public:
    Printer();
    virtual ~Printer();
    void PrintDToConsole( const Array2d & a, int skip_bounds = 0);
    void PrintDToFile(std::string & filename, Array2d & a, int skip_bounds = 0);
    void PrintDToPlotFile(std::string & filename, Array2d & a,
                          double x_lb, double y_lb, double dx, double dy, double skip_bounds=0) ;
    void PrintVelToConsole( Array3d & a, int skip_bounds = 0);
    void PrintVelToPlotFile( std::string & ufile, std::string & vfile,  Array3d & a,
                             double x_lb, double y_lb, double dx, double dy,  int skip_bounds= 0);

};

#endif /* PRINT_H_ */
