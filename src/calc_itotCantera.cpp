/*
 * calc_itotCantera.cpp
 *
 *  Created on: 15 May 2018
 *      Author: Manik
 */

#include <iostream>
#include "calc_itotCantera.h"

Cantera::Interface *surf;

void initCantera()
{
	// import the reference bulk phases
	Cantera::ThermoPhase* elyte = (Cantera::newPhase(p_inputFile, p_nameElectrolytePhase));
	Cantera::ThermoPhase* cond = (Cantera::newPhase(p_inputFile, p_nameElectrodePhase));

	std::vector<Cantera::ThermoPhase*> phaseListCA;
	phaseListCA.push_back(elyte);
	phaseListCA.push_back(cond);

	// Import the surface
	surf = Cantera::importInterface(p_inputFile, p_nameReactionSurfName, phaseListCA);
}

void calc_ropCantera2S(doublereal cA, doublereal cB, doublereal t, doublereal rop[nSpecies+1], doublereal phis)
{
	std::string s = "";
	size_t idxElyte = 0, idxElode = 0;

	try
	{
		// Initialize phases and phase index
		for (size_t k = 0; k < surf->nPhases(); k++)
		{
			surf->thermo(k).setState_TP(T,P);
			if (surf->thermo(k).name().compare(p_nameElectrolytePhase)==0)
			{
				idxElyte = k;
			}
			if (surf->thermo(k).name().compare(p_nameElectrodePhase)==0)
			{
				idxElode = k;
			}
		}

		// Electrolyte phase data
		doublereal cTotal = 0.0, rhoElyte, mmwElyte;

		size_t nS = surf->thermo(idxElyte).nSpecies();
		Cantera::vector_fp xT(nS);
		surf->thermo(idxElyte).getMoleFractions(xT.data());

		rhoElyte = surf->thermo(idxElyte).density();
		mmwElyte = surf->thermo(idxElyte).meanMolecularWeight();

		// Update mole fractions
		for (size_t kk = 0; kk < nS; kk++)
		{
			if (surf->thermo(idxElyte).speciesName(kk).compare(p_nameSpeciesA)==0)
				cTotal += cA;
			else if (surf->thermo(idxElyte).speciesName(kk).compare(p_nameSpeciesB)==0)
				cTotal += cB;
			else
				cTotal += xT[kk]*rhoElyte/mmwElyte;
		}
		for (size_t kk = 0; kk < nS; kk++)
		{
			if (surf->thermo(idxElyte).speciesName(kk).compare(p_nameSpeciesA)==0)
				xT[kk] = cA/cTotal;
			else if (surf->thermo(idxElyte).speciesName(kk).compare(p_nameSpeciesB)==0)
				xT[kk] = cB/cTotal;
			else
				xT[kk] = (xT[kk]*rhoElyte/mmwElyte)/cTotal; // Solvent
		}
		surf->thermo(idxElyte).setMoleFractions(xT.data());

		// Update electrode potential
		surf->thermo(idxElode).setElectricPotential(phis);

		// Calculate reaction rates
		Cantera::vector_fp wdot(surf->thermo(idxElode).nSpecies() + surf->thermo(idxElyte).nSpecies() + surf->nSpecies());
		surf->getNetProductionRates(wdot.data());
		rop[0] = wdot[surf->kineticsSpeciesIndex(p_nameSpeciesA)]*1e3;
		rop[1] = wdot[surf->kineticsSpeciesIndex(p_nameSpeciesB)]*1e3;
		rop[2] = wdot[surf->kineticsSpeciesIndex("electron")]*1e3;
	}
	catch (Cantera::CanteraError& err)
	{
		printf(err.what());
		Cantera::appdelete();
	}
}
