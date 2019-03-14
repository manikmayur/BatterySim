/*
 * calc_itotCantera.cpp
 *
 *  Created on: 15 May 2018
 *      Author: Manik
 */

#include "cantera/canteraFunctions.h"
#include "singleParticleModel/parameters_SPM.h"
#include <iostream>
typedef struct
{
	Cantera::Interface *surface;
	phaseType pT;
} surfPhase;

static surfPhase newSurf;
static Cantera::Interface *surfaceCA;
static Cantera::Interface *surfaceAN;

void initCanteraSPM()
{
	// Import the reference bulk phases
    try
    {
    Cantera::ThermoPhase* cathode = Cantera::newPhase(p_inputFile, p_nameCathodePhase);
	Cantera::ThermoPhase* anode = Cantera::newPhase(p_inputFile, p_nameAnodePhase);
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
    catch (Cantera::CanteraError& err)
    {
            // handle exceptions thrown by Cantera
            std::cout << err.what() << std::endl;
            std::cout << " terminating... " << std::endl;
            Cantera::appdelete();
    }
}

double calc_ilocCantera(double phiS)
{
	double rop = 0.0;
	size_t idx = 0;
	double elodeArea = (newSurf.pT == phaseType::phCA) ? p_S_ca:p_S_an;

	idx = getPhaseIdxbyName(p_nameConductorPhase, newSurf.surface);

	// Set the electrode and electrolyte potential
	newSurf.surface->thermo(idx).setElectricPotential(phiS);

	// Calculate reaction rates
	Cantera::vector_fp wdot(newSurf.surface->nTotalSpecies());
	newSurf.surface->getNetProductionRates(wdot.data());

	// Get the net reaction rate at the cathode-side interface
	rop = wdot[newSurf.surface->kineticsSpeciesIndex("electron")]*1e3; // [mol/m2/s]

	return rop*F*elodeArea; // [A]
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

double calc_potCantera(std::string surfName, double xLi, double phiL, double I, double T)
{
	size_t idxElyte, idxElode;
	double phiS, guess = 1.0;
	Cantera::compositionMap speciesMoleFrac;
	try
	{
		if (surfName.compare(p_nameCathodeSurf)==0)
		{
			newSurf.surface = surfaceCA;
			newSurf.pT = phaseType::phCA;
			speciesMoleFrac[p_nameCathodeIntSpecies] = xLi;
			speciesMoleFrac[p_nameCathodeVacSpecies] = 1-xLi;
		}
		else if (surfName.compare(p_nameAnodeSurf)==0)
		{
			newSurf.surface = surfaceAN;
			newSurf.pT = phaseType::phAN;
			speciesMoleFrac[p_nameAnodeIntSpecies] = xLi;
			speciesMoleFrac[p_nameAnodeVacSpecies] = 1-xLi;
		}
		for (size_t k = 0; k < newSurf.surface->nPhases(); k++)
		{
			newSurf.surface->thermo(k).setState_TP(T,Pref);
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
		phiS = NewtonSolver(&calc_ilocCantera, guess, -I);

		return phiS;
	}
	catch (Cantera::CanteraError& err)
	{
		printf(err.what());
		Cantera::appdelete();
	}
	return 0.0;
}

double calc_entropyCantera(std::string surfName, double xLi, double T)
{
	size_t idxElode;
	Cantera::compositionMap speciesMoleFrac;

	try
	{
		if (surfName.compare(p_nameCathodeSurf)==0)
		{
			newSurf.surface = surfaceCA;
			newSurf.pT = phaseType::phCA;
			speciesMoleFrac[p_nameCathodeIntSpecies] = xLi;
			speciesMoleFrac[p_nameCathodeVacSpecies] = 1-xLi;
		}
		else if (surfName.compare(p_nameAnodeSurf)==0)
		{
			newSurf.surface = surfaceAN;
			newSurf.pT = phaseType::phAN;
			speciesMoleFrac[p_nameAnodeIntSpecies] = xLi;
			speciesMoleFrac[p_nameAnodeVacSpecies] = 1-xLi;
		}
		for (size_t k = 0; k < newSurf.surface->nPhases(); k++)
		{
			newSurf.surface->thermo(k).setState_TP(T,Pref);
			if (newSurf.surface->thermo(k).name().compare(p_nameCathodePhase)==0 || newSurf.surface->thermo(k).name().compare(p_nameAnodePhase)==0)
			{
				idxElode = k;
			}
		}
		Cantera::vector_fp entropy(newSurf.surface->thermo(idxElode).nSpecies());
		newSurf.surface->thermo(idxElode).setMoleFractionsByName(speciesMoleFrac);
		newSurf.surface->thermo(idxElode).getEntropy_R(entropy.data());
		entropy[0] = entropy[0]*Cantera::GasConstant
				- Cantera::GasConstant*std::log(xLi/(1-xLi));
		return entropy[0]*1e-3 ; // Convert from J/kmol/K to J/mol/K
	}
	catch (Cantera::CanteraError& err)
	{
		printf(err.what());
		Cantera::appdelete();
	}
	return 0.0;
}
