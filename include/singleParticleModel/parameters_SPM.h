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
enum domainType {AL,CA,EL,AN,CU};

int const F = 96485; // [C/mol] "Faraday constant"
double const R = 8.314; // [J/(mol.K)] "Universal gas constant"
double const Tref = 298.15; // [K] Reference temperature
double const Pref = 101325.0; // [Pa] Pressure

std::string const paramFile = "data/parameters_SPM.yaml"; // XML file name
YAML::Node const static params = YAML::LoadFile(paramFile);

// Load model parameters
modelType const p_model = static_cast<modelType>(params["model"].as<int>()); // Model type

size_t const p_NT = params["tSteps"].as<double>();
size_t const p_NR = params["rSteps"].as<double>();
double const RTOL = params["aTol"].as<double>();    /* scalar relative tolerance */
double const ATOL = params["rTol"].as<double>();      /* scalar absolute tolerance */
YAML::Node const grids = params["grid"];

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
double const p_tTotal = params["tTotal"].as<double>(); // [s] "Total runtime"
double const p_Tamb = params["T_amb"].as<double>(); // [K] "Ambient temperature"
double const p_Vcut = params["V_cutoff"].as<double>(); // [V] "Cutoff voltage"

double const p_h = 1; //[W/m2K]
double const p_brugg = 0.4;

// Cell parameters
double const p_Lsep = params["L_sep"].as<double>(); // [m] "Separator thickness"
double const p_Lca = params["L_ca"].as<double>();   // [m] "Cathode thickness"
double const p_Lan = params["L_an"].as<double>();   // [m] "Anode thickness"
double const p_Lcell = p_Lca + p_Lsep + p_Lan;   // [m] "Cell thickness"
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
double const p_eps_ca = 0.3;

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
double const p_tp = params["tP"].as<double>(); // [1] Transference number

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

typedef struct
{
	domainType domType;
	size_t NX,idx0,idxL;
	double L;
	double dx;
	double rho;
	double por;
	double epsF;
	double cP;
	double kappaL;
	double sigmaL;
	double diffL;
	double kappaS;
	double sigmaS;
	double diffS;
	double k;
	double actEK;
	double actED;
	double rP;
	double cLiMax;
	double cLiInit;
	double aLi;
	double Cdl;
} domain;

domain const al =
{
	.domType = domainType::AL,
	.NX = grids[0].as<size_t>(),
	.idx0 = 0,
	.idxL = al.NX-1,
	.L = params["Aluminium"]["thickness"].as<double>(),
	.dx = al.L/al.NX,
	.rho = params["Aluminium"]["density"].as<double>(),
	.por = params["Aluminium"]["porosity"].as<double>(),
	.epsF = 0.0,
	.cP = params["Aluminium"]["cP"].as<double>(),
	.kappaL = 0.0,
	.sigmaL = 0.0,
	.diffL= 0.0,
	.kappaS = params["Aluminium"]["kappaS"].as<double>(),
	.sigmaS = params["Aluminium"]["sigmaS"].as<double>(),
	.diffS= 0.0,
	.k = 0.0,
	.actEK = 0.0,
	.actED = 0.0,
	.rP= 0.0,
	.cLiMax = 0.0,
	.cLiInit = 0.0,
	.aLi = 0.0,
	.Cdl = 0.0
};

domain const ca =
{
	.domType = domainType::CA,
	.NX = grids[1].as<size_t>(),
	.idx0 = al.idxL+1,
	.idxL = al.idxL+ca.NX,
	.L = params["Cathode"]["thickness"].as<double>(),
	.dx = ca.L/ca.NX,
	.rho = params["Cathode"]["density"].as<double>(),
	.por = params["Cathode"]["porosity"].as<double>(),
	.epsF = params["Cathode"]["fillerFrac"].as<double>(),
	.cP = params["Cathode"]["cP"].as<double>(),
	.kappaL = 0.0,
	.sigmaL = 0.0,
	.diffL= p_DLi,
	.kappaS = params["Cathode"]["kappaS"].as<double>(),
	.sigmaS = params["Cathode"]["sigmaS"].as<double>(),
	.diffS = params["Cathode"]["diffS"].as<double>(),
	.k = params["Cathode"]["rateConst"].as<double>(),
	.actEK = params["Cathode"]["actEK"].as<double>(),
	.actED = params["Cathode"]["actED"].as<double>(),
	.rP = params["Cathode"]["rP"].as<double>(),
	.cLiMax = params["Cathode"]["cLiMax"].as<double>(),
	.cLiInit = ((1-params["InitSOC"].as<double>()/100)*
			(params["Cathode"]["xLiMax"].as<double>()-params["Cathode"]["xLiMin"].as<double>())
			+params["Cathode"]["xLiMin"].as<double>())*ca.cLiMax,
	.aLi = params["Cathode"]["aLi"].as<double>(),
	.Cdl = params["Cathode"]["Cdl"].as<double>()
};

domain const el =
{
	.domType = domainType::EL,
	.NX = grids[2].as<size_t>(),
	.idx0 = ca.idxL+1,
	.idxL = ca.idxL+el.NX,
	.L = params["Separator"]["thickness"].as<double>(),
	.dx = el.L/el.NX,
	.rho = params["Separator"]["density"].as<double>(),
	.por = params["Separator"]["porosity"].as<double>(),
	.epsF = 0.0,
	.cP = params["Separator"]["cP"].as<double>(),
	.kappaL = 0.0,
	.sigmaL = 0.0,
	.diffL = p_DLi,
	.kappaS = params["Separator"]["kappaS"].as<double>(),
	.sigmaS = params["Separator"]["sigmaS"].as<double>(),
	.diffS = 0.0,
	.k = 0.0,
	.actEK = 0.0,
	.actED = 0.0,
	.rP = 0.0,
	.cLiMax = 0.0,
	.cLiInit = 0.0,
	.aLi = 0.0,
	.Cdl = 0.0
};

domain const an =
{
	.domType = domainType::AN,
	.NX = grids[3].as<size_t>(),
	.idx0 = el.idxL+1,
	.idxL = el.idxL+an.NX,
	.L = params["Anode"]["thickness"].as<double>(),
	.dx = an.L/an.NX,
	.rho = params["Anode"]["density"].as<double>(),
	.por = params["Anode"]["porosity"].as<double>(),
	.epsF = params["Anode"]["fillerFrac"].as<double>(),
	.cP = params["Anode"]["cP"].as<double>(),
	.kappaL = 0.0,
	.sigmaL = 0.0,
	.diffL = p_DLi,
	.kappaS = params["Anode"]["kappaS"].as<double>(),
	.sigmaS = params["Anode"]["sigmaS"].as<double>(),
	.diffS = params["Anode"]["diffS"].as<double>(),
	.k = params["Anode"]["rateConst"].as<double>(),
	.actEK = params["Anode"]["actEK"].as<double>(),
	.actED = params["Anode"]["actED"].as<double>(),
	.rP = params["Anode"]["rP"].as<double>(),
	.cLiMax = params["Anode"]["cLiMax"].as<double>(),
	.cLiInit = (params["InitSOC"].as<double>()/100*
			(params["Anode"]["xLiMax"].as<double>()-params["Anode"]["xLiMin"].as<double>())
			+params["Anode"]["xLiMin"].as<double>())*an.cLiMax,
	.aLi = params["Anode"]["aLi"].as<double>(),
	.Cdl = params["Anode"]["Cdl"].as<double>()
};

domain const cu =
{
	.domType = domainType::CU,
	.NX = grids[4].as<size_t>(),
	.idx0 = an.idxL+1,
	.idxL = an.idxL+cu.NX,
	.L = params["Copper"]["thickness"].as<double>(),
	.dx = an.L/an.NX,
	.rho = params["Copper"]["density"].as<double>(),
	.por = params["Copper"]["porosity"].as<double>(),
	.epsF = 0.0,
	.cP = params["Copper"]["cP"].as<double>(),
	.kappaL = 0.0,
	.sigmaL = 0.0,
	.diffL = 0.0,
	.kappaS = params["Copper"]["kappaS"].as<double>(),
	.sigmaS = params["Copper"]["sigmaS"].as<double>(),
	.diffS = 0.0,
	.k = 0.0,
	.actEK = 0.0,
	.actED = 0.0,
	.rP = 0.0,
	.cLiMax = 0.0,
	.cLiInit = 0.0,
	.aLi = 0.0,
	.Cdl = 0.0
};

size_t const NXTOTAL = al.NX+ca.NX+el.NX+an.NX+cu.NX;

inline domain getDomain(size_t jx)
{
	if (jx >= al.idx0 && jx <= al.idxL) return al;
	if (jx >= ca.idx0 && jx <= ca.idxL) return ca;
	if (jx >= el.idx0 && jx <= el.idxL) return el;
	if (jx >= an.idx0 && jx <= an.idxL) return an;
	if (jx >= cu.idx0 && jx <= cu.idxL) return cu;
	else std::cout<<"getDomain: Unknown grid index "<<jx<<std::endl;
	return al;
}
inline double dx(size_t ix)
{
	return getDomain(ix).dx;
}
// Thermal conductivity solid
inline double kappaS(size_t ix)
{
	return getDomain(ix).kappaS;
}
// Electric conductivity solid
inline double sigmaS(size_t ix)
{
	return getDomain(ix).sigmaS*(1-getDomain(ix).por-getDomain(ix).epsF);
}
// Electric conductivity liquid
inline double sigmaL(double por, double ce, double T)
{
	return std::pow(por,p_brugg)*1e-4*ce*std::pow(-10.5+0.668*1e-3*ce+0.494*1e-6*std::pow(ce,2)
			+(0.074-1.78*1e-5*ce-8.86*1e-10*std::pow(ce,2))*T+(-6.96*1e-5+2.8*1e-8*ce)*std::pow(T,2),2);
}
// Diffusivity Li in solid
inline double diffS(size_t ix, double T)
{
	return getDomain(ix).diffS*std::exp(-getDomain(ix).actED/R*(1/T-1/Tref));
}
// Diffusivity Li in liquid
inline double diffL(double por,double Ce, double T)
{
	return std::pow(por,p_brugg)*1e-4*std::pow(10,((-4.43-54/(T-229-Ce*5e-3)-Ce*0.22e-3)));
}
inline double rateConst(size_t ix, double T)
{
	return getDomain(ix).k*std::exp(-getDomain(ix).actEK/R*(1/T-1/Tref));
}

#endif /* PARAMETERS_SPM_H_ */
