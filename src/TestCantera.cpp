//============================================================================
// Name        : TestCantera.cpp
// Author      : Manik Mayur
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

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

void thermoInitAnode_demo(std::string xmlFile);
void thermoInitCathodeROP_demo(std::string xmlFile);
void thermoInitAnodeROP_demo(std::string xmlFile);
void thermoInitRefdemo(std::string xmlFile);

int main() {
	std::string xmlFile;
	std::cout<<"Testing Cantera"<<std::endl;
	//xmlFile = "Work_LiO2_organic_LiBaLu_November_1M_parallel_reactions.xml";
	xmlFile = "Intercalation_LCO_Graphite.xml";
	//xmlFile = "Final_Kupper_2016_JElectrochemSoc_LFP_C6_revised.xml";
	try {
		//thermoInitAnode_demo(xmlFile);
		//thermoInitAnodeROP_demo(xmlFile);
		//thermoInitRefdemo(xmlFile);
		thermoInitCathodeROP_demo(xmlFile);
		//thermoInitElectrolyte_demo(xmlFile);
		//simple_demo2();
	}
	catch (CanteraError& err){
		std::cout<<err.what()<< std::endl;
	}
	return 0;
}

void thermoInitAnode_demo(std::string xmlFile) {
	std::string s = "";
	double x=0.0;
	double c[4];
	double ns[2];
	double T=298.15;
	double P=101325;
	double minTemp, maxTemp, refPressure;
	double mu = 0.0;
	int type = 0;
	try {
		ThermoPhase* tp = (newPhase(xmlFile,"anode"));
		IdealSolidSolnPhaseTabulatedThermo* Li_ion = dynamic_cast<IdealSolidSolnPhaseTabulatedThermo*> (tp);
		//ConstDensityTabulatedThermo* Li_ion = dynamic_cast<ConstDensityTabulatedThermo*> (tp);
		//std::cout<<Li_ion->report();
		//SpeciesThermo& sp = Li_ion->speciesThermo(); //Cantera 2.2
		Li_ion->setState_TP(T,P);

		MultiSpeciesThermo& sp = Li_ion->speciesThermo();
		size_t index = Li_ion->m_kk_mod;
		type = sp.reportType(index);

		//std::cout<<Li_ion->report();
		//Li_ion->getChemPotentials(&mu);
		//printf("\nCantera: c0=%f, mu=%f\n",Li_ion->standardConcentration(Li_ion->m_kk_int), mu);
		//printf("\nInterpolation: At x=0.5, s=%f, sc=%f, xlx=%f\n",GasConstant*Li_ion->mean_X(Li_ion->entropy_R()), Li_ion->interp_s(x)*0.5e3-GasConstant*log(0.5), GasConstant*Li_ion->sum_xlogx());
		// Ask for mole fraction
		printf("Enter a mole fraction:\n");
		std::cin>>x;
		ns[0]=x;
		ns[1]=1-x;
		Li_ion->setMoleFractions(ns);
		std::cout<<Li_ion->report();

		sp.reportParams(index, type, c, minTemp, maxTemp, refPressure);
		printf("\nSpecies Thermo: At x=%f, h=%f, s=%f, E0=%f\n",x, c[1],c[2], 1e-3/96485*(-c[1]+298.15*c[2]));

		printf("\nValues from data: At x=%f, h=%f, s=%f\n",x, Li_ion->interp_h(x)*1e3,Li_ion->interp_s(x)*1e3);
		//printf("\nMethods of Thermoclass: At x=%f, s=%f, scalc=%f, xlx=%f\n",x, GasConstant*Li_ion->mean_X(Li_ion->getEntropy_R()), Li_ion->interp_s(x)*1e3+GasConstant*log(x/(1-x)), GasConstant*Li_ion->sum_xlogx());
		delete tp;
	} catch (const std::exception& ex) {
		std::cout<<ex.what()<<std::endl;
	}
}

void printData(ThermoPhase* tp) {
	double c[4];
	double minTemp, maxTemp, refPressure;
	double T=298.15;
	vector_fp Xo(tp->nSpecies());
	vector_fp G(tp->nSpecies());
	vector_fp H(tp->nSpecies());
	vector_fp S(tp->nSpecies());
	/* Cantera 2.2
	tp->getMoleFractions(DATA_PTR(Xo));
	tp->getGibbs_RT(DATA_PTR(G));
	tp->getEnthalpy_RT(DATA_PTR(H));
	tp->getEntropy_R(DATA_PTR(S));
	 */
	tp->getMoleFractions(Xo.data());
	tp->getGibbs_RT(G.data());
	tp->getEnthalpy_RT(H.data());
	tp->getEntropy_R(S.data());
	for (size_t k = 0; k < tp->nSpecies(); k++) {
		int type = tp->speciesThermo(k).reportType(k);
		tp->speciesThermo(k).reportParams(k,type,c,minTemp, maxTemp, refPressure);
		std::cout << tp->speciesName(k) <<" X= "<< Xo[k] <<" DH= " << c[1]<<" DS= " << c[2]<< std::endl;
		std::cout << tp->speciesName(k) <<" DH= "<< H[k]*GasConstant*T <<" DS= " << S[k]*GasConstant<<" DG= " << G[k]*GasConstant*T<< std::endl;
	}
}
void thermoInitCathodeROP_demo(std::string xmlFile) {
	std::string s = "";
	double ns[2];
	double T=298.15;
	double P=101325;
	double x, dG, dH, dS, dG0, f;
	try {
		ThermoPhase* tp = (newPhase(xmlFile,"electrolyte"));
		ThermoPhase* tp1 = (newPhase(xmlFile,"cathode"));
		ThermoPhase* tp2 = (newPhase(xmlFile,"conductor"));
		std::vector<ThermoPhase*> phaseList;
		phaseList.push_back(tp);
		phaseList.push_back(tp1);
		phaseList.push_back(tp2);
		Edge surf(xmlFile, "interface_cathode", phaseList);
		//Interface surf(xmlFile, "interface_cathode", phaseList);
		std::cout<<surf.phaseIndex("cathode")<<" "<<surf.phaseIndex("electrolyte")<<" "<<surf.phaseIndex("electron_cathode");
		tp->setState_TP(T,P);
		tp1->setState_TP(T,P);
		tp2->setState_TP(T,P);
		surf.setState_TP(T,P);
		x=0.1;
		ns[0]=x;
		ns[1]=1-x;
		tp->setMoleFractions(ns);
		for (size_t k = 0; k < phaseList.size(); k++) {
			phaseList[k]->setState_TP(T,P);
		}
		surf.setState_TP(T,P);
		surf.getDeltaGibbs(&dG);
		surf.getDeltaSSGibbs(&dG0);
		surf.getDeltaEnthalpy(&dH);
		surf.getDeltaEntropy(&dS);
		std::cout<<"No: "<<surf.nReactions()<<" "<<surf.reactionString(0)<<std::endl;
		vector_fp wdot(surf.nTotalSpecies());
		surf.getNetProductionRates(wdot.data());
		std::cout<<"dG= "<<dG0/1000<<" dH= "<<dH/1000<<" dS= "<<dS/1000<<std::endl;
		f = std::exp(0.5/(GasConstant*T)*(-dG0));
		//surf.setMultiplier(0,f);
		//surf.getNetProductionRates(wdot.data());
		//for (size_t k = 0; k < phaseList.size(); k++) {
			//std::cout << "Name: "<<phaseList[k]->name()<<" density: "<<phaseList[k]->density()<< std::endl;
			for (size_t kk = 0; kk < surf.nTotalSpecies(); kk++)
			std::cout << surf.kineticsSpeciesName(kk)<< " wdot = " << wdot[kk] << std::endl;
		//}
	} catch (CanteraError& err){
		std::cout<<err.what()<< std::endl;
	}
}

void thermoInitAnodeROP_demo(std::string xmlFile) {
	std::string s = "";
	std::string phName, spName;
	double ns[2];
	double T=298.15;
	double P=101325;
	double x,dG[20],dH[20],dS[20],dG0[20],dH0[20],dS0[20];
	double minTemp, maxTemp, refPressure;
	double c[4];
	double Ec1, Ec2;
	int type;

	ThermoPhase* tp = (newPhase(xmlFile,"anode"));
	ThermoPhase* tp1 = (newPhase(xmlFile,"electrolyte"));
	ThermoPhase* tp2 = (newPhase(xmlFile,"electron"));

	std::vector<ThermoPhase*> phaseList;

	phaseList.push_back(tp);
	phaseList.push_back(tp1);
	phaseList.push_back(tp2);

	Edge surf(xmlFile, "edge_anode_electrolyte", phaseList);

	//Interface surf(xmlFile, "edge_anode_electrolyte", phaseList);
	//Interface refSurf(xmlFile, "edge_lithium_electrolyte", refPhaseList);
	std::cout<<surf.phaseIndex("anode")<<" "<<surf.phaseIndex("electrolyte")<<" "<<surf.phaseIndex("electron_anode");
	tp->setState_TP(T,P);
	tp1->setState_TP(T,P);
	tp2->setState_TP(T,P);
	surf.setState_TP(T,P);

	MultiSpeciesThermo& sp = tp->speciesThermo();
	type = sp.reportType(0);
	sp.reportParams(0, type, c, minTemp, maxTemp, refPressure);
	printf("\nSpecies Thermo: At x=%f, h=%f, s=%f, E0=%f\n",tp->moleFraction(0), c[1],c[2], 1e-3/96485*(-c[1]+298.15*c[2]));
	x=0.4;
	ns[0]=x;
	ns[1]=1-x;
	tp->setMoleFractions(ns);
	std::cout<<tp->report();

	//std::cout<<surf.report();
	for(size_t i = 0; i < surf.nReactions(); i++){
		std::cout<<surf.reaction(i)->id << std::endl;
	}
	for (size_t k=0; k < tp->nSpecies(); k++) {
		type = sp.reportType(k);
		sp.reportParams(k, type, c, minTemp, maxTemp, refPressure);
		printf("\nSpecies Thermo: At x=%f, h=%f, s=%f, E0=%f\n",x, c[1],c[2], 1e-3/96485*(-c[1]+298.15*c[2]));
	}

	std::cout<<"No: "<<surf.nReactions()<<" "<<surf.reactionString(0)<<std::endl;
	for (size_t k = 0; k < surf.nPhases(); k++) {
		phName = surf.thermo(k).id();
		surf.thermo(k).getEnthalpy_RT(dH0);
		surf.thermo(k).getEntropy_R(dS0);
		surf.thermo(k).getPartialMolarEnthalpies(dH);
		surf.thermo(k).getPartialMolarEntropies(dS);
		//dH = surf.thermo(k).enthalpy_mole();
		//dS = surf.thermo(k).entropy_mole();
		for (size_t n = 0; n < surf.thermo(k).nSpecies(); n++) {
			spName = surf.thermo(k).speciesName(n);
			std::cout<< "Phase: " <<  phName << " species: "<< spName << " dH0 = " << dH0[n]*GasConstant*T << " dS0 = " << dS0[n]*GasConstant <<std::endl;
			std::cout<< "Phase: " <<  phName << " species: "<< spName << " partial dH = " << dH[n] << " partial dS = " << dS[n] <<std::endl;
		}
	}
	vector_fp wdot(tp->nSpecies() + tp1->nSpecies() +tp2->nSpecies() +surf.nSpecies());
	surf.getNetProductionRates(wdot.data());
	surf.getDeltaGibbs(dG);
	surf.getDeltaEnthalpy(dH);
	surf.getDeltaEntropy(dS);
	surf.getDeltaSSGibbs(dG0);
	surf.getDeltaSSEnthalpy(dH0);
	surf.getDeltaSSEntropy(dS0);

	//f = std::exp(0.5/(GasConstant*T)*(-dG0[0]));
	//surf.setMultiplier(0,f);
	surf.getNetProductionRates(wdot.data());
	std::cout<<"dG0= "<< dG0[0] <<" dH0= "<< dH0[0] <<" dS0= "<< dS0[0] << std::endl;
	std::cout<<"dG= " << dG[0] << " dH= "<< dH[0] <<" dS= "<< dS[0] << std::endl;
	std::cout<<"dGc= "<< dH[0] - T*dS[0] << " Ec=" << -dG[0]/96485e3 << " E0=" << -dG0[0]/96485e3 << std::endl;
	Ec1 = -dG[0]/96485e3;
	/*
	for (size_t n = 0; n < surf.nPhases(); n++) {
		std::cout<< "Name: " << surf.thermo(n).id()<< " dS: " << surf.thermo(n).getPartialMolarEntropies() << std::endl;
	}
	for (size_t k = 0; k < tp1->nSpecies(); k++) {
		std::cout<< tp1->speciesName(k) <<" wdot = " << wdot[k] << std::endl;
	}
	for (size_t k = 0; k < tp->nSpecies(); k++) {
		std::cout << tp->speciesName(k)<< " wdot = " << wdot[k+tp1->nSpecies()] << std::endl;
	}*/
	for (size_t k = 0; k < tp2->nSpecies(); k++) {
		std::cout << tp2->speciesName(k)<< " wdot = " << wdot[k+tp1->nSpecies()+tp->nSpecies()] <<" speciesIdx = " <<k+tp1->nSpecies()+tp->nSpecies() << std::endl;
	}
}

void thermoInitRefdemo(std::string xmlFile) {

	std::string s = "";
	std::string phName, spName;
	double ns[2];
	double T=298.15;
	double P=101325;
	double x,dG[20],dH[20],dS[20],dG0[20],dH0[20],dS0[20];
	double minTemp, maxTemp, refPressure;
	double c[4];
	double Ec1, Ec2;
	int type;

	std::vector<ThermoPhase*> refPhaseList;

	ThermoPhase* tp1 = (newPhase(xmlFile,"electrolyte"));
	ThermoPhase* tp2 = (newPhase(xmlFile,"electron"));
	ThermoPhase* tp = (newPhase(xmlFile,"lithium_metal"));

	refPhaseList.push_back(tp);
	refPhaseList.push_back(tp1);
	refPhaseList.push_back(tp2);

	Edge refSurf(xmlFile, "edge_reference_electrode", refPhaseList);

	tp->setState_TP(T,P);
	tp1->setState_TP(T,P);
	tp2->setState_TP(T,P);
	refSurf.setState_TP(T,P);

	for(size_t i = 0; i < refSurf.nReactions(); i++){
		std::cout<<refSurf.reaction(i)->id << std::endl;
	}

	// Reference electrode reaction
	std::cout<< "Reference electrode reaction" << std::endl;
	std::cout<<"No: "<<refSurf.nReactions()<<" "<<refSurf.reactionString(0)<<std::endl;
	for (size_t k = 0; k < refSurf.nPhases(); k++) {
		phName = refSurf.thermo(k).id();
		refSurf.thermo(k).getEnthalpy_RT(dH0);
		refSurf.thermo(k).getEntropy_R(dS0);
		refSurf.thermo(k).getPartialMolarEnthalpies(dH);
		refSurf.thermo(k).getPartialMolarEntropies(dS);
		//dH = refSurf.thermo(k).enthalpy_mole();
		//dS = refSurf.thermo(k).entropy_mole();
		for (size_t n = 0; n < refSurf.thermo(k).nSpecies(); n++) {
			spName = refSurf.thermo(k).speciesName(n);
			std::cout<< "Phase: " <<  phName << " species: "<< spName << " dH0 = " << dH0[n]*GasConstant*T << " dS0 = " << dS0[n]*GasConstant <<std::endl;
			std::cout<< "Phase: " <<  phName << " species: "<< spName << " partial dH = " << dH[n] << " partial dS = " << dS[n] <<std::endl;
		}
	}
	vector_fp wdot(tp->nSpecies() + tp1->nSpecies() +tp2->nSpecies() +refSurf.nSpecies());
	refSurf.getNetProductionRates(wdot.data());
	refSurf.getDeltaGibbs(dG);
	refSurf.getDeltaEnthalpy(dH);
	refSurf.getDeltaEntropy(dS);
	refSurf.getDeltaSSGibbs(dG0);
	refSurf.getDeltaSSEnthalpy(dH0);
	refSurf.getDeltaSSEntropy(dS0);

	std::cout<<"dG0= "<< dG0[0] <<" dH0= "<< dH0[0] <<" dS0= "<< dS0[0] << std::endl;
	std::cout<<"dG= " << dG[0] << " dH= "<< dH[0] <<" dS= "<< dS[0] << std::endl;
	std::cout<<"dGc= "<< dH[0] - T*dS[0] << " Ec=" << -dG[0]/96485e3 << " E0=" << -dG0[0]/96485e3 << std::endl;

	for (size_t k = 0; k < tp2->nSpecies(); k++) {
		std::cout << tp2->speciesName(k)<< " wdot = " << wdot[k+tp1->nSpecies()+tp->nSpecies()] <<" speciesIdx = " <<k+tp1->nSpecies()+tp->nSpecies() << std::endl;
	}
}
