/*
 * parameters.h
 *
 *  Created on: 7 Mar 2018
 *      Author: Manik
 */

#ifndef PARAMETERS_2S_H_
#define PARAMETERS_2S_H_

#include <string>
#include <math.h>

int const F = 96485; // [C/mol] "Faraday constant"
double const R = 8.314; // [J/(mol.K)] "Universal gas constant"
double const T = 298.15; // [K] Temperature
double const P = 101325.0; // [Pa] Pressure
int const nSpecies = 2; // [C/mol] "Number of species"
// Cantera parameters
std::string const inputFile = "cantera/Ferrocene_CV.xml"; // XML file name
std::string const reactionSurfName = "WE_surface";
std::string const electrolytePhaseName = "electrolyte";
std::string const electrodePhaseName = "conductor";
std::string const cA_name = "Ferrocene[elyte]";
std::string const cB_name = "Ferrocene+[elyte]";
// Transport paramters
double const cA_bulk = 1e-3; // [mol/L] "Reactant bulk concentration"
double const cB_bulk = 1e-6; // [mol/L] "Product bulk concentration"
double const c_ref = 1e-3; // [mol/L] "Reactant bulk concentration"
double const DA = 4e-10; // [m^2/s] "Reactant diffusion coefficient"
double const DB = 4e-10; // [m^2/s] "Product diffusion coefficient"
double const Cdl = 0.2*0; // [F/m^2] "Double layer capacity"
// Cyclic voltammetry parameters
double const v = 0.01; // [V/s] "Voltammetric scan rate"
double const E1 = 0; // [V] "Start potential"
double const E2 = 1; // [V] "Switching potential"
double const tp = (E2-E1)/v ; // [s] "Peak time"
double const L = 6*sqrt(DA*2*std::abs(E1-E2)/v); // [m] "Outer bound on diffusion layer"
unsigned int const n_scp = 3; // "Number of scans before measurement"
unsigned int const n_sc = 1; // "Number of scans, measurement"
// Butler-Volmer parameters
unsigned int const n = 1; // [1] "Number of electrons transferred"
double const Eeq = 0.2; // [V] "Equilibrium potential"
double const i0 = 10.0; // [A/m^2] "Exchange current density"
double const k0 = i0/F; // "Reaction rate"
double const alpha_a = 0.5; // [1] "Anodic transfer coefficient"
double const alpha_c = 0.5; // [1] "Cathodic transfer coefficient"
int const vA = 1; // [1] "Valence of the species A"
int const vB = -1; // [1] "Valence of the species B"

#define Vcell(t) ((t<=tp) ? v*t + E1 : -v*(t-tp) + E2)

#endif /* PARAMETERS_2S_H_ */
