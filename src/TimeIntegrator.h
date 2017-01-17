#ifndef TIMEINTEGRATOR_H_
#define TIMEINTEGRATOR_H_
#include "SimParams.h"
#include "GenSimDefs.h"

class TimeIntegrator {
protected:
    SimParams * sp;
    double total_mass_outflux;
public:
    TimeIntegrator(SimParams * sp);
    virtual ~TimeIntegrator();
    double GetTotalMassOutflux();
    virtual int Integrate(Array2d & c, double const * const q, const double dt) = 0;
    int IntegrateSubcycle(Array2d & c, double const * const q, const double dt);
};

#endif /* TIMEINTEGRATOR_H_ */
