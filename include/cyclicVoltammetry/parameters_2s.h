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
#include <iostream>
#include "yaml-cpp/yaml.h"

int const Faraday = 96485; // [C/mol] "Faraday constant"
double const gasConstant = 8.314; // [J/(mol.K)] "Universal gas constant"
double const T = 298.15; // [K] Temperature
double const P = 101325.0; // [Pa] Pressure
int const nSpecies = 2; // [C/mol] "Number of species"
std::string const paramFile = "data/parameters_2s.yaml"; // XML file name
YAML::Node const static params = YAML::LoadFile(paramFile);

// Load cantera parameters
std::string const p_inputFile = params["inputFile"].as<std::string>(); // XML file name
std::string const p_nameReactionSurfName = params["reactionSurfName"].as<std::string>();
std::string const p_nameElectrolytePhase = params["electrolytePhaseName"].as<std::string>();
std::string const p_nameElectrodePhase = params["electrodePhaseName"].as<std::string>();
std::string const p_nameSpeciesA = params["cA_name"].as<std::string>();
std::string const p_nameSpeciesB = params["cB_name"].as<std::string>();

// Transport paramters
double const p_c0A = params["cA_bulk"].as<double>(); // [mol/L] "Reactant bulk concentration"
double const p_c0B = params["cB_bulk"].as<double>(); // [mol/L] "Product bulk concentration"
double const p_c0 = params["c_ref"].as<double>(); // [mol/L] "Reactant bulk concentration"
double const p_DA = params["DA"].as<double>(); // [m^2/s] "Reactant diffusion coefficient"
double const p_DB = params["DB"].as<double>(); // [m^2/s] "Product diffusion coefficient"
double const p_cDL = params["Cdl"].as<double>(); // [F/m^2] "Double layer capacity"
// Cyclic voltammetry parameters
double const p_scanRate = params["v"].as<double>(); // [V/s] "Voltammetric scan rate"
double const p_startVoltage = params["E1"].as<double>(); // [V] "Start potential"
double const p_switchVoltage = params["E2"].as<double>(); // [V] "Switching potential"
double const p_tP = std::abs((p_switchVoltage-p_startVoltage)/p_scanRate) ; // [s] "Peak time"
double const p_L = 6*sqrt(2*p_DA*p_tP); // [m] "Outer bound on diffusion layer"
unsigned int const p_nScp = params["n_scp"].as<unsigned int>(); // "Number of scans before measurement"
unsigned int const p_nSc = params["n_sc"].as<unsigned int>(); // "Number of scans, measurement"
// Butler-Volmer parameters
unsigned int const p_nElectrons = params["n"].as<unsigned int>(); // [1] "Number of electrons transferred"
double const p_Eeq = params["Eeq"].as<double>(); // [V] "Equilibrium potential"
double const p_i0 = params["i0"].as<double>(); // [A/m^2] "Exchange current density"
double const p_k0 = p_i0/Faraday; // "Reaction rate"
double const p_alphaA = params["alpha_a"].as<double>(); // [1] "Anodic transfer coefficient"
double const p_alphaC = params["alpha_c"].as<double>(); // [1] "Cathodic transfer coefficient"
int const p_vA = params["vA"].as<int>(); // [1] "Valence of the species A"
int const p_vB = params["vB"].as<int>(); // [1] "Valence of the species B"

#define Vcell(t) ((t<=p_tP) ? p_scanRate*t + p_startVoltage : -p_scanRate*(t-p_tP) + p_switchVoltage)

#endif /* PARAMETERS_2S_H_ */
