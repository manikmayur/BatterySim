/*
 * calc_itotCantera.cpp
 *
 *  Created on: 15 May 2018
 *      Author: Manik
 */

#include <iostream>
#include "calc_itotCantera.h"

void calc_itotCantera2s(doublereal cA, doublereal cB, doublereal t, doublereal rop[nSpecies+1], doublereal phis) {
	std::string s = "";

	try
	{
		// Electrode phases

		// import the reference bulk phases
		Cantera::ThermoPhase* elyte = (Cantera::newPhase(inputFile,"electrolyte"));
		Cantera::ThermoPhase* cond = (Cantera::newPhase(inputFile,"conductor"));

		std::vector<Cantera::ThermoPhase*> phase_ca;
		phase_ca.push_back(elyte);
		phase_ca.push_back(cond);

		// Import the surface
		Cantera::Interface surf(inputFile, "WE_surface", phase_ca);

		// Set the state
		cond->setState_TP(T,P);
		elyte->setState_TP(T,P);
		surf.setState_TP(T,P);

		// Electrolyte phase data
		doublereal cTotal = 0.0, rho_elyte, mmw_elyte;

		size_t nS = elyte->nSpecies();
		Cantera::vector_fp xT(nS);
		elyte->getMoleFractions(xT.data());

		rho_elyte = elyte->density();
		mmw_elyte = elyte->meanMolecularWeight();

		// Update mole fractions
		for (size_t kk = 0; kk < nS; kk++)
		{
			if (elyte->speciesName(kk).compare(cA_name)==0)
				cTotal += cA;
			else if (elyte->speciesName(kk).compare(cB_name)==0)
				cTotal += cB;
			else
				cTotal += xT[kk]*rho_elyte/mmw_elyte;
		}
		for (size_t kk = 0; kk < nS; kk++)
		{
			if (elyte->speciesName(kk).compare(cA_name)==0)
				xT[kk] = cA/cTotal;
			else if (elyte->speciesName(kk).compare(cB_name)==0)
				xT[kk] = cB/cTotal;
			else
				xT[kk] = (xT[kk]*rho_elyte/mmw_elyte)/cTotal; // Solvent
		}
		elyte->setMoleFractions(xT.data());

		// Update electrode potential
		cond->setElectricPotential(phis);

		// Calculate reaction rates
		Cantera::vector_fp wdot(cond->nSpecies() + elyte->nSpecies() + surf.nSpecies());
		surf.getNetProductionRates(wdot.data());
		rop[0] = wdot[surf.kineticsSpeciesIndex(cA_name)]*1e3;
		rop[1] = wdot[surf.kineticsSpeciesIndex(cB_name)]*1e3;
		rop[2] = wdot[surf.kineticsSpeciesIndex("electron")]*1e3;
		std::cout<<surf.kineticsSpeciesIndex(cA_name)<<" "<<rop[0]<<" "<<surf.kineticsSpeciesIndex(cB_name)<<" "<<rop[1]<<std::endl;
	}
	catch (Cantera::CanteraError& err){
		printf(err.what());
		Cantera::appdelete();
	}
}
