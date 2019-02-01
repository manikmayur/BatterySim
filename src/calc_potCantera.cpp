/*
 * calc_itotCantera.cpp
 *
 *  Created on: 15 May 2018
 *      Author: Manik
 */

#include <src/calc_potCantera.h>
#include <iostream>
#include "parameters_SPM.h"
typedef struct
{
	Cantera::Interface * surface;
	phaseType pT;
} surfPhase;

static surfPhase newSurf;
static Cantera::Interface *surfaceCA;
static Cantera::Interface *surfaceAN;

void initCanteraSPM()
{
	// import the reference bulk phases
	Cantera::ThermoPhase* anode = Cantera::newPhase(p_inputFile, p_nameAnodePhase);
	Cantera::ThermoPhase* cathode = Cantera::newPhase(p_inputFile, p_nameCathodePhase);
	Cantera::ThermoPhase* cond = Cantera::newPhase(p_inputFile, p_nameConductorPhase);
	Cantera::ThermoPhase* elyte = Cantera::newPhase(p_inputFile, p_nameElectrolytePhase);

	std::vector<Cantera::ThermoPhase*> phaseListCathode;
	std::vector<Cantera::ThermoPhase*> phaseListAnode;
	phaseListCathode.push_back(cathode);
	phaseListCathode.push_back(cond);
	phaseListCathode.push_back(elyte);

	phaseListAnode.push_back(anode);
	phaseListAnode.push_back(cond);
	phaseListAnode.push_back(elyte);

	// Import the surfaces
	surfaceCA = Cantera::importInterface(p_inputFile, p_nameCathodeSurf, phaseListCathode);
	surfaceAN = Cantera::importInterface(p_inputFile, p_nameAnodeSurf, phaseListAnode);
}

double calc_ilocCantera(double phiS)
{
	double rop = 0.0;
	size_t idx = 0;
	double elodeArea = (newSurf.pT == phaseType::CA) ? p_S_ca:p_S_an;

	idx = getPhaseIdxbyName(p_nameConductorPhase, newSurf.surface);

	// Set the electrode and electrolyte potential
	newSurf.surface->thermo(idx).setElectricPotential(phiS);

	// Calculate reaction rates
	Cantera::vector_fp wdot(newSurf.surface->nTotalSpecies());
	newSurf.surface->getNetProductionRates(wdot.data());

	// Get the net reaction rate at the cathode-side interface
	rop = wdot[newSurf.surface->kineticsSpeciesIndex("electron")]*1e3; // [mol/m2/s]

	return rop*Faraday*elodeArea; // [A]
}

size_t getPhaseIdxbyName(std::string phaseName, Cantera::Interface *surface)
{
	size_t idx;
	for (size_t k = 0; k < surface->nPhases(); k++)
	{
		if (surface->thermo(k).name().compare(phaseName)==0)
		{
			idx = k;
			return idx;
		}
	}
	return -1;
}

double calc_potCantera(double phiL, std::string surfName, Cantera::compositionMap speciesMoleFrac)
{
	std::string s = "";
	size_t idxElyte, idxElode;
	double phiS = 1.0;

	try
	{
		if (surfName.compare(p_nameCathodeSurf)==0)
		{
			newSurf.surface = surfaceCA;
			newSurf.pT = phaseType::CA;
		}
		else if (surfName.compare(p_nameAnodeSurf)==0)
		{
			newSurf.surface = surfaceAN;
			newSurf.pT = phaseType::AN;
		}
		for (size_t k = 0; k < newSurf.surface->nPhases(); k++)
		{
			newSurf.surface->thermo(k).setState_TP(T,P);
			if (newSurf.surface->thermo(k).name().compare(p_nameElectrolytePhase)==0)
			{
				idxElyte = k;
			}
			if (newSurf.surface->thermo(k).name().compare(p_nameCathodePhase)==0 || newSurf.surface->thermo(k).name().compare(p_nameAnodePhase)==0)
			{
				idxElode = k;
			}
		}

		newSurf.surface->thermo(idxElode).setMoleFractionsByName(speciesMoleFrac);

		// Set the electrode and electrolyte potential
		newSurf.surface->thermo(idxElyte).setElectricPotential(phiL);
		phiS = NewtonSolver(&calc_ilocCantera, 0.1, p_Iapp);

		return phiS;
	}
	catch (Cantera::CanteraError& err)
	{
		printf(err.what());
		Cantera::appdelete();
	}
	return 0.0;
}
