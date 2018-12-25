/*
 * calc_itotCantera.cpp
 *
 *  Created on: 15 May 2018
 *      Author: Manik
 */

#include <iostream>
#include <string>
#include "calc_itotCantera.h"

using namespace Cantera;

void calc_itotCantera2s(doublereal cA, doublereal cB, doublereal t, doublereal &itot, doublereal phis) {
	std::string s = "";
	double T=298.15, P=101325.0;
	std::string inputFile = "/Users/Manik/CloudStation/Development/workspace/TestSunDIALS/Work_LiO2_organic_LiBaLu_CV_TDPA.xml"; // XML file name

	try {
	// Electrode phases

	// import the reference bulk phases
	ThermoPhase* cond = (newPhase(inputFile,"conductor"));
	ThermoPhase* elyte = (newPhase(inputFile,"elyte"));

	std::vector<ThermoPhase*> phase_ca;
	phase_ca.push_back(elyte);
	phase_ca.push_back(cond);

	cond->setState_TP(T,P);
	elyte->setState_TP(T,P);
	//printf(elyte->report().c_str());

	// Import the surface
	Interface surf(inputFile, "C_surface", phase_ca);

	surf.setState_TP(T,P);
	vector_fp xT(elyte->nSpecies());
	std::vector<std::string> nS(elyte->nSpecies());
	elyte->getMoleFractions(xT.data());
	nS = elyte->speciesNames();
	//Tab = table(nS, xT)

	// Reference phase
	ThermoPhase* Li_metal = (newPhase(inputFile,"lithium"));
	std::vector<ThermoPhase*> phase_ref;
	phase_ref.push_back(Li_metal);
	phase_ref.push_back(elyte);
	phase_ref.push_back(cond);

	//Edge ref_surf(inputFile, "Li_surface", phase_ref);
	//thermoInitCathodeROP_demo(xmlFile);

	// Initialization
	doublereal cS, cSC, cSA, rho_elyte, mmw_elyte;

	rho_elyte = elyte->density();
	mmw_elyte = elyte->meanMolecularWeight();

	cS = xT[1]*rho_elyte/mmw_elyte;
	cSC = xT[2]*rho_elyte/mmw_elyte;
	cSA = xT[3]*rho_elyte/mmw_elyte;

	xT[1] = cS/(cSC+cSA+cA+cB); // Solvent
	xT[2] = cSC/(cSC+cSA+cA+cB); // Salt cation
	xT[3] = cSA/(cSC+cSA+cA+cB); // Salt anion
	xT[4] = cA/(cSC+cSA+cA+cB); // RM oxidised
	xT[5] = cB/(cSC+cSA+cA+cB); // RM reduced

	elyte->setMoleFractions(xT.data());
	cond->setElectricPotential(phis);

	vector_fp wdot(cond->nSpecies() + elyte->nSpecies() + surf.nSpecies());
	surf.getNetProductionRates(wdot.data());
/*
	for (size_t kk = 0; kk < surf.nTotalSpecies(); kk++)
		std::cout << surf.kineticsSpeciesName(kk)<<" wdot= "<<wdot[kk]<<std::endl;*/
	}
	catch (CanteraError& err){
		printf(err.what());
		Cantera::appdelete();
	}
}

void calc_itotCantera3s(doublereal cA, doublereal cB, doublereal cC, doublereal t, doublereal &itot, doublereal phis) {
	std::string s = "";
	double T=298.15, P=101325.0;
	std::string inputFile = "/Users/Manik/CloudStation/Development/workspace/TestSunDIALS/Work_LiO2_organic_LiBaLu_CV_TDPA.xml"; // XML file name

	try {
	// Electrode phases

	// import the reference bulk phases
	ThermoPhase* cond = (newPhase(inputFile,"conductor"));
	ThermoPhase* elyte = (newPhase(inputFile,"elyte"));

	std::vector<ThermoPhase*> phase_ca;
	phase_ca.push_back(elyte);
	phase_ca.push_back(cond);

	cond->setState_TP(T,P);
	elyte->setState_TP(T,P);
	printf(elyte->report().c_str());

	// Import the surface
	Interface surf(inputFile, "C_surface", phase_ca);

	surf.setState_TP(T,P);
	vector_fp xT(elyte->nSpecies());
	std::vector<std::string> nS(elyte->nSpecies());
	elyte->getMoleFractions(xT.data());
	nS = elyte->speciesNames();

	// Reference phase
	ThermoPhase* Li_metal = (newPhase(inputFile,"lithium"));
	std::vector<ThermoPhase*> phase_ref;
	phase_ref.push_back(Li_metal);
	phase_ref.push_back(elyte);
	phase_ref.push_back(cond);

	//Edge ref_surf(inputFile, "Li_surface", phase_ref);
	//thermoInitCathodeROP_demo(xmlFile);

	// Initialization
	doublereal cS, cSC, cSA, rho_elyte, mmw_elyte;

	rho_elyte = elyte->density();
	mmw_elyte = elyte->meanMolecularWeight();

	cS = xT[1]*rho_elyte/mmw_elyte;
	cSC = xT[2]*rho_elyte/mmw_elyte;
	cSA = xT[3]*rho_elyte/mmw_elyte;

	xT[1] = cS/(cSC+cSA+cA+cB+cC); // Solvent
	xT[2] = cSC/(cSC+cSA+cA+cB+cC); // Salt cation
	xT[3] = cSA/(cSC+cSA+cA+cB+cC); // Salt anion
	xT[4] = cA/(cSC+cSA+cA+cB+cC); // RM oxidised
	xT[5] = cB/(cSC+cSA+cA+cB+cC); // RM reduced
	xT[6] = cC/(cSC+cSA+cA+cB+cC); // RM reduced

	elyte->setMoleFractions(xT.data());
	cond->setElectricPotential(phis);

	vector_fp wdot(cond->nSpecies() + elyte->nSpecies() + surf.nSpecies());
	surf.getNetProductionRates(wdot.data());

	for (size_t kk = 0; kk < surf.nTotalSpecies(); kk++)
		std::cout << surf.kineticsSpeciesName(kk)<<" wdot= "<<wdot[kk]<<std::endl;
	}
	catch (CanteraError& err){
		printf(err.what());
		Cantera::appdelete();
	}
}
