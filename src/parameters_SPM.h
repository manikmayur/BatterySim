/*
 * parameters_SPM.h
 *
 *  Created on: Jan 27, 2019
 *      Author: Manik
 */

#ifndef PARAMETERS_SPM_H_
#define PARAMETERS_SPM_H_

#include <string>
#include <math.h>
#include <iostream>
#include "yaml-cpp/yaml.h"

int const Faraday = 96485; // [C/mol] "Faraday constant"
double const gasConstant = 8.314; // [J/(mol.K)] "Universal gas constant"
double const T = 298.15; // [K] Temperature
double const Tref = 298.15; // [K] Reference temperature
double const P = 101325.0; // [Pa] Pressure
int const nSpecies = 2; // [C/mol] "Number of species"

std::string const paramFile = "src/parameters_SPM.yaml"; // XML file name
YAML::Node const static params = YAML::LoadFile(paramFile);

// Load cantera parameters
std::string const p_inputFile = params["inputFile"].as<std::string>(); // XML file name
std::string const p_nameCathodePhase = params["cathodePhaseName"].as<std::string>();
std::string const p_nameCathodeIntSpecies = params["cathodeIntSpeciesName"].as<std::string>();
std::string const p_nameCathodeVacSpecies = params["cathodeVacSpeciesName"].as<std::string>();
std::string const p_nameAnodePhase = params["anodePhaseName"].as<std::string>();
std::string const p_nameAnodeIntSpecies = params["anodeIntSpeciesName"].as<std::string>();
std::string const p_nameAnodeVacSpecies = params["anodeVacSpeciesName"].as<std::string>();
std::string const p_nameElectrolytePhase = params["electrolytePhaseName"].as<std::string>();
std::string const p_nameConductorPhase = params["conductorPhaseName"].as<std::string>();
std::string const p_nameCathodeSurf = params["cathodeSurfaceName"].as<std::string>();
std::string const p_nameAnodeSurf = params["anodeSurfaceName"].as<std::string>();

// Load parameters
double const p_I1C = params["I1C"].as<double>(); // [A] "1C discharge current"
double const p_cR = params["cR"].as<double>(); // [1] "C-rate"
double const p_Iapp = p_cR*p_I1C; // [A] "Applied current"

// Cell geometry parameters
double const p_L_sep = params["L_sep"].as<double>(); // [m] "Separator thickness"
double const p_L_ca = params["L_ca"].as<double>();                      // [m] "Cathode thickness"
double const p_L_an = params["L_an"].as<double>();                    // [m] "Anode thickness"
double const p_L_cell = p_L_ca + p_L_sep + p_L_an; // [m] "Cell thickness"

// Cathode material parameters
double const p_rP_ca = params["rP_ca"].as<double>(); // [m] "Particle radius cathode"
double const p_xLimax_ca = params["xLimax_ca"].as<double>();  // [1] "Maximum cathode stoichiometry"
double const p_xLimin_ca = params["xLimin_ca"].as<double>(); // [1] "Minimum cathode stoichiometry"
double const p_DLiref_ca = params["DLiref_ca"].as<double>();  // [m^2/s] "Solid phase Li-diffusivity LCO"
double const p_Ediff_ca = params["Ediff_ca"].as<double>(); // [kJ/mol] "Cathode diffusion activation energy"
double const p_DLi_ca = p_DLiref_ca*std::exp(p_Ediff_ca/gasConstant*(1/T-1/Tref)); // [m^2/s] "Solid phase Li-diffusivity LCO"
double const p_csMax_ca = params["csMax_ca"].as<double>(); // [mol/m^3] "Max solid phase concentration cathode"
double const p_S_ca = params["S_ca"].as<double>();

// Anode material parameters
double const p_rP_an = params["rP_an"].as<double>();    // [m] "Particle radius anode"
double const p_xLimax_an = params["xLimax_an"].as<double>();    // [1] "Maximum anode stoichiometry"
double const p_xLimin_an = params["xLimin_an"].as<double>();    // [1] "Minimum anode stoichiometry"
double const p_DLiref_an = params["DLiref_an"].as<double>(); // [m^2/s] "Solid phase Li-diffusivity LMO"
double const p_Ediff_an = params["Ediff_an"].as<double>(); // [kJ/mol] "Anode diffusion activation energy"
double const p_DLi_an = p_DLiref_an*std::exp(p_Ediff_an/gasConstant*(1/T-1/Tref)); // [m^2/s] "Solid phase Li-diffusivity LMO"
double const p_csMax_an = params["csMax_an"].as<double>(); //[mol/m^3] "Max solid phase concentration anode"
double const p_S_an = params["S_an"].as<double>();

// Electrolyte material parameters
double const p_cE = params["cE"].as<double>(); // [mol/m^3] "Electrolyte concentration"
double const p_cE_ref = params["cE_ref"].as<double>(); // [mol/m^3] "Electrolyte reference concentration"
/*p_theta1 = (-5.636e-7*p_Iapp - 7.283e-6)*(p_Tref - 273.15)^3 +...
    (5.676e-5*p_Iapp + 6.453e-4)*(p_Tref - 273.15)^2 + ...
    (-2.221e-3*p_Iapp - 1.635e-2)*(p_Tref - 273.15) + ...
    (2.437e-2*p_Iapp + 1.428e-1);
p_theta2 = (-6.824e-6*p_Iapp + 1.372e-5)*(p_Tref - 273.15)^3 +
    (6.054e-4*p_Iapp - 1.216e-3)*(p_Tref - 273.15)^2 +
    (-1.497e-2*p_Iapp + 3.025e-2)*(p_Tref - 273.15) +
    (7.179e-2*p_Iapp - 1.456e-1);
p_Rel = p_theta1 + p_theta2*(p_T - p_Tref); */
#endif /* PARAMETERS_SPM_H_ */
