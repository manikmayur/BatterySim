/*
 * algorithms.h
 *
 *  Created on: 31.01.2019
 *      Author: Manik
 */

#ifndef ALGORITHMS_H_
#define ALGORITHMS_H_


double zeroin(				// An estimate to the root
	const double ax,		// Specify the interval the root
	const double bx,		// to be sought in
	double (* f) (const double),	// Function under investigation
	const double tol);		// Acceptable tolerance

double NewtonSolver(double (* f) (const double), double xstart, double C);
double fzerotx(double (*f) (const double), double ax, double bx);


#endif /* ALGORITHMS_H_ */
