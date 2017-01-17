#ifndef TRINTEGRATOR_H_
#define TRINTEGRATOR_H_

#include "TimeIntegrator.h"
#include "SimParams.h"
#include "GenSimDefs.h"

class TRIntegrator : public TimeIntegrator {
public:
    TRIntegrator(SimParams *sp);
    virtual ~TRIntegrator();
    virtual int Integrate(Array2d & c, double const * const q, const double dt);
};

#endif /* FEINTEGRATOR_H_ */
