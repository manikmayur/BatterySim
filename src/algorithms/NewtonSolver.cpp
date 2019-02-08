#include "algorithms/algorithms.h"
#include <math.h>


double NewtonSolver(double (* f) (const double), double xstart, double C=0.0)
{
    //Solve f(x) = C by Newton iteration.
    // xstart    starting point for Newton iteration
    // C         constant
    double ff, dfdx, f0 = f(xstart) - C;
    double emax, x0 = xstart;
    double step, dx = 1.0e-6;
    unsigned int n = 0;
    while (n < 200)
    {
        ff = f(x0 + dx) - C;
        dfdx = (ff - f0)/dx;
        step = - f0/dfdx;

        // avoid taking steps too large
        if (std::abs(step) > 0.1)
            step = 0.1*step/std::abs(step);

        x0 += step;
        emax = 0.00001;  // 0.01 mV tolerance
        if (std::abs(f0) < emax && n > 8)
            return x0;
        f0 = f(x0) - C;
        n += 1;
    }
    return -1.0;
}
