/*
 * algorithms.h
 *
 *  Created on: 31.01.2019
 *      Author: Manik
 */

#ifndef ALGORITHMS_H_
#define ALGORITHMS_H_

#include <vector>
#include <math.h>
#include <boost/numeric/odeint.hpp>


double zeroin(				// An estimate to the root
	const double ax,		// Specify the interval the root
	const double bx,		// to be sought in
	double (* f) (const double),	// Function under investigation
	const double tol);		// Acceptable tolerance

double NewtonSolver(double (* f) (const double), double xstart, double C);
double fzerotx(double (*f) (const double), double ax, double bx);
double interpolate(double x, std::vector<double> xList, std::vector<double> fList);
double integrate(double (*f) (const double, const double), double x0, double t0, double t1, size_t N=100);
double integrate(double (*f) (const double), double x0, double t0, double t1, size_t N=100);

/*inline double integral(double (* f) (const double), double t0, double t1)
{
	auto my_system = [](const double &x,
						const double &dxdt,
						const double t) {dxdt[0] = f(t);};
	std::vector<double> x0(1); // Initial condition, vector of 1 element (scalar problem)
    x0[0] = 0;
    // Observer, prints time and state when called (during integration)
    auto my_observer = [](const double &x, const double t )
    		{std::cout<< t << " " << x[0] << "\n";};
    // Integration parameters
    double dt = (t1-t0)/100;

    // Run integrator
    odeint::integrate(my_system, x0, t0, t1, dt, my_observer);
	return 0;
}*/

#endif /* ALGORITHMS_H_ */
