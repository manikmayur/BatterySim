/*
 * integration.cpp
 *
 *  Created on: 09.02.2019
 *      Author: Manik
 */
#include "algorithms/algorithms.h"
namespace odeint = boost::numeric::odeint;

 // Defining a shorthand for the type of the mathematical state
 typedef std::vector< double > state_type;

 // This is a ode system functor
 // System to be solved: dx/dt = f(x,t)
 struct systemInt
 {
	 systemInt(double (*f)(const double)) : f(f) {}
	 void operator()(const state_type &x , state_type &dxdt , const double t) const
	 {dxdt[0] = f(t); }

 private:
	 double (*f)(const double);
 };

 struct systemInt2
 {
	 systemInt2(double (*f)(const double, const double)) : f(f) {}
	 void operator()(const state_type &x , state_type &dxdt , const double t) const
	 {dxdt[0] = f(t,x[0]); }

 private:
	 double (*f)(const double, const double);
 };

 // This is a observer functor for output processing
 struct observer
 {
	 observer(double *r) : result(r) {}
	 void operator()(const state_type &x, const double t)
	 {
		 *result = x[0];
	 }
 private:
	 double *result;
 };

double integrate(double (*f) (const double, const double), double x0, double t0, double t1, size_t N)
{
	double result;
	double dt = (t1-t0)/N;
	systemInt2 odeSystem(f);
	observer odeObserver(&result);
	state_type x(1);
	x[0] = x0;
	// Run integrator
	odeint::integrate(odeSystem, x, t0, t1, dt, odeObserver);
	return result;
}

double integrate(double (*f) (const double), double x0, double t0, double t1, size_t N)
{
	double result;
	double dt = (t1-t0)/N;
	systemInt odeSystem(f);
	observer odeObserver(&result);
	state_type x(1);
	x[0] = x0;
	// Run integrator
	odeint::integrate(odeSystem, x, t0, t1, dt, odeObserver);
	return result;
}
