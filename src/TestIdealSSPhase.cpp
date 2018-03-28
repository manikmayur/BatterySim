/*
 * TestIdealSSPhase.cpp
 *
 *  Created on: 23.02.2018
 *      Author: mmayur
 */

#include "cantera/thermo/ConstDensityTabulatedThermo.h"
#include "cantera/thermo/IdealSolidSolnPhaseTabulatedThermo.h"
#include <iostream>
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/Interface.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/kinetics/importKinetics.h"
#include "cantera/Edge.h"
#include <string>

using namespace Cantera;

void thermoIDSSElectrolyte_demo(std::string xmlFile);
void thermoIDSSCathodeROP_demo(std::string xmlFile);
void thermoOxygenSolubility(std::string xmlFile);

int main() {
	std::string xmlFile;
	std::cout<<"Testing Cantera"<<std::endl;
	xmlFile = "Work_LiO2_organic_LiBaLu_November_1M_parallel_reactions.xml";
	try {
		//thermoIDSSElectrolyte_demo(xmlFile);
		//thermoIDSSCathodeROP_demo(xmlFile);
		thermoOxygenSolubility(xmlFile);
	}
	catch (CanteraError& err){
		std::cout<<err.what()<< std::endl;
	}
	return 0;
}

void thermoIDSSElectrolyte_demo(std::string xmlFile) {
	std::string s = "";
	double T=298.15;
	double P=101325;
	double rho_mix[4] = {1100.0, 1121.0, 1163.0, 1217.0};
	double conc_LiTFSI[4] = {0.0, 0.1, 0.5, 1.0};
	double conc_O2[4] = {0.00033, 0.00023, 0.000209, 0.000252};
	double sum_conc = 0.0;
	try {
		ThermoPhase* tp1 = (newPhase(xmlFile,"elyte"));
		std::vector<ThermoPhase*> phaseList;
		phaseList.push_back(tp1);
		tp1->setState_TP(T,P);
		vector_fp conc(tp1->nSpecies(),0.0);
		vector_fp xx(tp1->nSpecies(),0.0);
		vector_fp mw(tp1->nSpecies(),0.0);
		tp1->getMolecularWeights(mw.data());;
		//
		for (unsigned int j = 0; j<4; j++) {
			sum_conc = 0.0;
			std::cout<<"\nCase: "<<j<<", conc LiTFSI:"<<conc_LiTFSI[j]<<std::endl;
			conc[1]=conc_LiTFSI[j];
			conc[2]=conc_LiTFSI[j];
			conc[3]=conc_O2[j];
			conc[0]=rho_mix[j]/mw[0]-conc[1]*mw[1]/mw[0]-conc[2]*mw[2]/mw[0]-conc[3]*mw[3]/mw[0];
			for (unsigned int jj = 0; jj<4; jj++)
				sum_conc += conc[jj];
			xx[0]=conc[0]/sum_conc;
			xx[1]=conc[1]/sum_conc;
			xx[2]=conc[2]/sum_conc;
			xx[3]=conc[3]/sum_conc;
			tp1->setMoleFractions(xx.data());
			for (size_t k = 0; k < tp1->nSpecies(); k++) {
				std::cout << tp1->speciesName(k)<< " conc = " << conc[k]<<", x = "<<xx[k]<< std::endl;
			}
			std::cout<<tp1->report();
			tp1->getConcentrations(conc.data());
			for (size_t k = 0; k < tp1->nSpecies(); k++) {
				std::cout << tp1->speciesName(k)<< " conc = " << conc[k]<< "MW = "<<mw[k]<< std::endl;
			}
		}

	}catch (CanteraError& err){
		std::cout<<err.what()<< std::endl;
	}
}

void thermoIDSSCathodeROP_demo(std::string xmlFile) {
	std::string s = "";
	double ns[2];
	double T=298.15;
	double P=101325;
	double x, dG, dH, dS, dG0;
	int totSpecies=0;
	std::vector<ThermoPhase*> phaseList;
	try {
		ThermoPhase* tp = (newPhase(xmlFile,"elyte"));
		//ThermoPhase* tp1 = (newPhase(xmlFile,"LiO2"));
		//ThermoPhase* tp2 = (newPhase(xmlFile,"Li2O2"));
		ThermoPhase* tp3 = (newPhase(xmlFile,"conductor"));

		phaseList.push_back(tp);
		//phaseList.push_back(tp1);
		//phaseList.push_back(tp2);
		phaseList.push_back(tp3);
		Interface surf(xmlFile, "C_surface", phaseList);
		for (size_t k = 0; k < phaseList.size(); k++) {
			phaseList[k]->setState_TP(T,P);
			totSpecies+=phaseList[k]->nSpecies();
		}
		surf.setState_TP(T,P);
		totSpecies+=surf.nSpecies();
		//std::cout<<tp->report();
		//std::cout<<surf.report();
		surf.getDeltaGibbs(&dG);
		surf.getDeltaEnthalpy(&dH);
		surf.getDeltaEntropy(&dS);
		std::cout<<"No: "<<surf.nReactions()<<" "<<surf.reactionString(0)<<std::endl;
		//surf.setMultiplier(0,f);
		vector_fp wdot(surf.nTotalSpecies());
		//surf.getNetProductionRates(DATA_PTR(wdot)); //Cantera 2.2
		surf.getNetProductionRates(wdot.data());
		surf.getDeltaSSGibbs(&dG0);
		surf.getDeltaEnthalpy(&dH);
		surf.getDeltaEntropy(&dS);
		std::cout<<"dG= "<<dG0/1000<<" dH= "<<dH/1000<<" dS= "<<dS/1000<<std::endl;
		for (size_t kk = 0; kk < surf.nTotalSpecies(); kk++)
			std::cout << surf.kineticsSpeciesName(kk)<< " wdot = " << wdot[kk] << std::endl;
		//for (size_t k = 0; k < phaseList.size(); k++) {
		//std::cout << "Name: "<<phaseList[k]->name()<<" density: "<<phaseList[k]->density()<< std::endl;


		//std::cout << phaseList[k]->speciesName(kk)<< " wdot = " << wdot[k+kk] << std::endl;
		//}

	} catch (CanteraError& err){
		std::cout<<err.what()<< std::endl;
	}
}

void thermoOxygenSolubility(std::string xmlFile) {
	std::string s = "";
	double T=298.15;
	double P=101325;
	double dG, dH, dS, dG0, mmu0;
	int totSpecies=0;
	std::vector<ThermoPhase*> phaseList;
	try {
		ThermoPhase* tp = (newPhase(xmlFile,"gas_cathode"));
		//ThermoPhase* tp1 = (newPhase(xmlFile,"LiO2"));
		//ThermoPhase* tp2 = (newPhase(xmlFile,"Li2O2"));
		ThermoPhase* tp1 = (newPhase(xmlFile,"elyte"));

		phaseList.push_back(tp);
		//phaseList.push_back(tp1);
		//phaseList.push_back(tp2);
		phaseList.push_back(tp1);
		Interface surf(xmlFile, "O_surface", phaseList);
		for (size_t k = 0; k < phaseList.size(); k++) {
			phaseList[k]->setState_TP(T,P);
			totSpecies+=phaseList[k]->nSpecies();
		}
		surf.setState_TP(T,P);
		totSpecies+=surf.nSpecies();
		surf.getDeltaGibbs(&dG);
		surf.getDeltaEnthalpy(&dH);
		surf.getDeltaEntropy(&dS);
		std::cout<<"No of reactions: "<<surf.nReactions()<<" #Species: "<<surf.nTotalSpecies()<<" #Phases: "<<surf.nPhases()<<std::endl;
		std::cout<<surf.reactionString(0)<<std::endl;
		//surf.setMultiplier(0,f);
		vector_fp wdot(surf.nTotalSpecies());
		vector_fp Kc(surf.nReactions());
		surf.getNetProductionRates(wdot.data());
		surf.getDeltaSSGibbs(&dG0);
		surf.getDeltaEnthalpy(&dH);
		surf.getDeltaEntropy(&dS);
		//! Equilibrium constant for all reactions including the voltage term
		/*!
		 *   Kc = exp(-(mu0P-mu0R)/RT)
		 *
		 *   where deltaG is the electrochemical potential difference between
		 *   products minus reactants.
		 */
		surf.getEquilibriumConstants(Kc.data());
		std::cout<<"Kc= "<<Kc[0]<<" dG= "<<dG0/1000<<" dH= "<<dH/1000<<" dS= "<<dS/1000<<std::endl;
		for (size_t k = 0; k < surf.nPhases(); k++) {
			std::cout << "Name: "<<surf.thermo(k).name()<<" #Species: "<< surf.thermo(k).nSpecies()<<" density: "<<surf.thermo(k).density()<<std::endl;
			vector_fp mu0(surf.thermo(k).nSpecies());
			vector_fp mu(surf.thermo(k).nSpecies());
			surf.thermo(k).getStandardChemPotentials(mu0.data());
			surf.thermo(k).getChemPotentials(mu.data());
			for (size_t kk = 0; kk < surf.thermo(k).nSpecies(); kk++) {
				mmu0 = mu0[kk] + Faraday * surf.thermo(k).electricPotential()*surf.thermo(k).charge(kk);
				mmu0 -= surf.thermo(0).RT() * surf.thermo(k).logStandardConc(k);
				std::cout <<"Species: "<< surf.thermo(k).speciesName(kk)<<" mmu0: "<<mmu0<<" mu0: "<<mu0[kk]<<" mu: "<<mu[kk]<<std::endl;
			}
		}
		for (size_t kk = 0; kk < surf.nTotalSpecies(); kk++){
			std::cout << surf.kineticsSpeciesName(kk)<< " wdot = " << wdot[kk]<< std::endl;
		}

	} catch (CanteraError& err){
		std::cout<<err.what()<< std::endl;
	}
}
