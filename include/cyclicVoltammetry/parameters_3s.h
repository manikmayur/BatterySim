/*
 * parameters.h
 *
 *  Created on: 7 Mar 2018
 *      Author: Manik
 */

#ifndef PARAMETERS_3S_H_
#define PARAMETERS_3S_H_

#include <string>
#include <cmath>


double const F = 96485; // [C/mol] "Faradays constant"
double const R = 8.314; // [J/(mol.K)] "Universal gas constant"
double const v = 0.1; // [V/s] "Voltammetric scan rate"
double const c_bulk = 1e-3; // [mol/L] "Reactant bulk concentration"
double const c_bulk_p1 = 0; // [mol/L] "Product bulk concentration"
double const c_bulk_p2 = 0; // [mol/L] "Product bulk concentration"
double const c_ref = 1e-3; // [mol/L] "Reactant bulk concentration"
double const DA = 1.0e-9; // [m^2/s] "Reactant diffusion coefficient"
double const DB = 1.0e-9; // [m^2/s] "Product diffusion coefficient"
double const DC = 1.0e-9; // [m^2/s] "Product diffusion coefficient"
double const Cdl = 0.2*0; // [F/m^2] "Double layer capacity"
double const T = 298.15; // [K] Temperature
double const E1 = -0.5; // [V] "Start potential"
double const E2 = 0.5; // [V] "Switching potential"
double const L = 6*sqrt(DA*2*std::abs(E1-E2)/v); // [m] "Outer bound on diffusion layer"
unsigned int const n_scp = 3; // "Number of scans before measurement"
unsigned int const n_sc = 1; // "Number of scans, measurement"
double const i0 = 10.0; // [A/m^2] "Exchange current density"
double const k0 = i0/F; // "Reaction rate"
double const alpha_a = 0.5; // [1] "Anodic transfer coefficient"
double const alpha_c = 0.5; // [1] "Cathodic transfer coefficient"
int const vA = 1; // [1] "Valence of the species A"
int const vB = -1; // [1] "Valence of the species B"
double const tp = (E2-E1)/v ; // [s] "Peak time"
unsigned int const n = 1; // [1] "Number of electrons transferred"
double const Eeq = 0.2; // [V] "Equilibrium potential"



#endif /* PARAMETERS_3S_H_ */
