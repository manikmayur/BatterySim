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

enum modelType {SPM,SPMT,SPMe,SPMeT};

int const Faraday = 96485; // [C/mol] "Faraday constant"
double const R = 8.314; // [J/(mol.K)] "Universal gas constant"
double const Tref = 298.15; // [K] Reference temperature
double const P = 101325.0; // [Pa] Pressure

std::string const paramFile = "data/parameters_SPM.yaml"; // XML file name
YAML::Node const static params = YAML::LoadFile(paramFile);

// Load model type
modelType const p_model = static_cast<modelType>(params["model"].as<int>()); // Model type

size_t const p_NT = params["tSteps"].as<double>();
size_t const p_NR = params["rSteps"].as<double>();
double const RTOL = params["aTol"].as<double>();    /* scalar relative tolerance */
double const ATOL = params["rTol"].as<double>();      /* scalar absolute tolerance */

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

// Load operating parameters
double const p_I1C = params["I1C"].as<double>(); // [A] "1C discharge current"
double const p_cR = params["cR"].as<double>(); // [1] "C-rate"
double const p_Iapp = p_cR*p_I1C; // [A] "Applied current"
double const p_tInit = params["tStart"].as<double>();
double const p_tTotal = 3600/abs(p_cR); // [s] "Total runtime"
double const p_Tamb = params["T_amb"].as<double>(); // [K] "Ambient temperature"

// Cell parameters
double const p_L_sep = params["L_sep"].as<double>(); // [m] "Separator thickness"
double const p_L_ca = params["L_ca"].as<double>();   // [m] "Cathode thickness"
double const p_L_an = params["L_an"].as<double>();   // [m] "Anode thickness"
double const p_L_cell = p_L_ca + p_L_sep + p_L_an;   // [m] "Cell thickness"
double const p_rhoCell = params["rhoCell"].as<double>(); // [kg/m3] "Cell density"
double const p_volCell = params["volCell"].as<double>(); // [m3] "Cell volume"
double const p_cpCell = params["cpCell"].as<double>(); // [J/(kg.K)] "Cell heat capacity"
double const p_hA = params["hA"].as<double>(); // [J/(s.K)] "Cell transfer coefficient"

// Cathode material parameters
double const p_rP_ca = params["rP_ca"].as<double>(); // [m] "Particle radius cathode"
double const p_xLimax_ca = params["xLimax_ca"].as<double>();  // [1] "Maximum cathode stoichiometry"
double const p_xLimin_ca = params["xLimin_ca"].as<double>(); // [1] "Minimum cathode stoichiometry"
double const p_DLiref_ca = params["DLiref_ca"].as<double>();  // [m^2/s] "Solid phase Li-diffusivity LCO"
double const p_Ediff_ca = params["Ediff_ca"].as<double>(); // [kJ/mol] "Cathode diffusion activation energy"
inline double const p_DLi_ca(double T) {return p_DLiref_ca*std::exp(p_Ediff_ca/R*(1/T-1/Tref));}; // [m^2/s] "Solid phase Li-diffusivity LCO"
double const p_csMax_ca = params["csMax_ca"].as<double>(); // [mol/m^3] "Max solid phase concentration cathode"
double const p_S_ca = params["S_ca"].as<double>();

// Anode material parameters
double const p_rP_an = params["rP_an"].as<double>();    // [m] "Particle radius anode"
double const p_xLimax_an = params["xLimax_an"].as<double>();    // [1] "Maximum anode stoichiometry"
double const p_xLimin_an = params["xLimin_an"].as<double>();    // [1] "Minimum anode stoichiometry"
double const p_DLiref_an = params["DLiref_an"].as<double>(); // [m^2/s] "Solid phase Li-diffusivity LMO"
double const p_Ediff_an = params["Ediff_an"].as<double>(); // [kJ/mol] "Anode diffusion activation energy"
inline double p_DLi_an(double T) {return p_DLiref_an*std::exp(p_Ediff_an/R*(1/T-1/Tref));}; // [m^2/s] "Solid phase Li-diffusivity LMO"
double const p_csMax_an = params["csMax_an"].as<double>(); //[mol/m^3] "Max solid phase concentration anode"
double const p_S_an = params["S_an"].as<double>();

// Electrolyte material parameters
double const p_cE = params["cE"].as<double>(); // [mol/m^3] "Electrolyte concentration"
double const p_cE_ref = params["cE_ref"].as<double>(); // [mol/m^3] "Electrolyte reference concentration"

double const p_theta1 = (-5.636e-7*abs(p_Iapp) - 7.283e-6)*pow((p_Tamb - 273.15),3)
   + (5.676e-5*abs(p_Iapp) + 6.453e-4)*pow((p_Tamb - 273.15),2)
   + (-2.221e-3*abs(p_Iapp) - 1.635e-2)*(p_Tamb - 273.15)
   + (2.437e-2*abs(p_Iapp) + 1.428e-1);

double const p_theta2 = (-6.824e-6*abs(p_Iapp) + 1.372e-5)*pow((p_Tamb - 273.15),3)
   + (6.054e-4*abs(p_Iapp) - 1.216e-3)*pow((p_Tamb - 273.15),2)
   + (-1.497e-2*abs(p_Iapp) + 3.025e-2)*(p_Tamb - 273.15)
   + (7.179e-2*abs(p_Iapp) - 1.456e-1);

inline double p_Rel(double T) {return p_theta1 + p_theta2*(T - p_Tamb);}

// Active material properties
double const p_DLi = params["DLi"].as<double>(); // 6e-14 #[m^2/s] Diffusion coefficient of Li in active material (Laresgoiti 2015 JPS)
double const p_cPmax = params["cPmax"].as<double>();// %[mol/m^3] Maximum concentration of Li in active material (Laresgoiti 2015 JPS)
double const p_cPinit = 0.5*p_cPmax;
double const p_alphaP = params["alphaP"].as<double>();
double const p_modEP = params["modEP"].as<double>();// [Pa] Youngs modulus of active material (Laresgoiti 2015 JPS)
double const p_nuP = params["nuP"].as<double>(); // % [1] Poisson's ratio of active material (Laresgoiti 2015 JPS)
double const p_omegaP = params["omegaP"].as<double>(); // [cm^3/mol] Partial molar volume of Li in active material (Laresgoiti 2015 JPS)
double const p_RP = params["RP"].as<double>(); // [m] Radius of the active material (Laresgoiti 2015 JPS)

// SEI properties
double const p_alphaS = params["alphaS"].as<double>();
double const p_modES = params["modES"].as<double>(); // [Pa] Youngs modulus of SEI (Laresgoiti 2015 JPS)
double const p_nuS = params["nuS"].as<double>(); // [1] Poisson's ratio of SEI (Laresgoiti 2015 JPS)
double const p_omegaS = params["omegaS"].as<double>(); // [cm^3/mol] Partial molar volume of Li in SEI (Laresgoiti 2015 JPS)
double const p_RS = p_RP + 0.2e-6; // [m] Radius of the SEI (Laresgoiti 2015 JPS)

#endif /* PARAMETERS_SPM_H_ */
