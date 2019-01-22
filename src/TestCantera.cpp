//============================================================================
// Name        : TestCantera.cpp
// Author      : Manik Mayur
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <cantera/thermo/BinarySolutionTabulatedThermo.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/Interface.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/kinetics/importKinetics.h"
#include "cantera/Edge.h"
#include "cantera/reactionpaths.h"
#include "calc_itotCantera.h"

void thermoInitAnode_demo(std::string inFile);
void thermoInitCathodePhases(std::string inFile);
void thermoInitAnodeROP_demo(std::string inFile);
void thermoInitRefdemo(std::string inFile);
void writeRxnPathDiagram(double time, Cantera::ReactionPathBuilder& b,
		Cantera::Kinetics& reaction, std::ostream& logfile, std::ostream& outfile);
void thermoTestSPM(std::string inFile);
void printThermoReactionO2(std::string inFile);
void printThermoReactionTDPA(std::string inFile);
void printThermoReactionLi2O2(std::string inFile);
void printThermoReactionLi(std::string inFile);
void printThermoReactionPEM(std::string inFile);
void printThermoReactionFerrocene(std::string inFile) ;
void printThermoKinetics(std::string inFile, std::string surfName, std::vector<Cantera::ThermoPhase*> phaseList, double T, double P);

int main() {
    std::string inFile;
    std::cout<<"Testing Cantera"<<std::endl;
    try {
        //thermoInitAnode_demo(inFile);
        //thermoInitAnodeROP_demo(inFile);
        //thermoInitRefdemo(inFile);
        //thermoInitCathodePhases(inFile);
        //thermoTestSPM(inFile);
        //thermoInitReactionROP_demo(inFile);
        //printThermoReactionO2("cantera\/Work_LiO2_organic_LiBaLu_CV_TDPA.cti");
        //printThermoReactionLi2O2("cantera\/Work_LiO2_organic_LiBaLu_CV_TDPA.cti");
        //printThermoReactionLi2O2("cantera\/Work_LiO2_organic_LiBaLu_November_2M_parallel_reactions.cti");
        //printThermoReactionTDPA("cantera\/Work_LiO2_organic_LiBaLu_CV_TDPA.cti");
        //printThermoReactionLi("cantera\/Work_LiO2_organic_LiBaLu_CV_TDPA.cti");
        //printThermoReactionPEM("cantera\/Final_Gruebl_2018_PEMFC_origin.cti");
        //thermoInitElectrolyte_demo(inFile);
        //simple_demo2();
    	printThermoReactionFerrocene("cantera\/Ferrocene_CV.xml");
    }
    catch (Cantera::CanteraError& err){
        std::cout<<err.what()<< std::endl;
    }
    return 0;
}

void thermoInitAnode_demo(std::string inFile) {
    std::string s = "";
    double x=0.0;
    double c[4];
    double ns[2];
    double T=298.15;
    double P=101325;
    double minTemp, maxTemp, refPressure;
    int type = 0;
    try {
        Cantera::ThermoPhase* tp = (Cantera::newPhase(inFile,"anode"));
        Cantera::BinarySolutionTabulatedThermo* Li_ion = dynamic_cast<Cantera::BinarySolutionTabulatedThermo*> (tp);
        Li_ion->setState_TP(T,P);

        Cantera::MultiSpeciesThermo& sp = Li_ion->speciesThermo();
        //size_t index = Li_ion->m_kk_mod;
        //type = sp.reportType(index);

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

        //sp.reportParams(index, type, c, minTemp, maxTemp, refPressure);
        printf("\nSpecies Thermo: At x=%f, h=%f, s=%f, E0=%f\n",x, c[1],c[2], 1e-3/96485*(-c[1]+298.15*c[2]));

        //printf("\nValues from data: At x=%f, h=%f, s=%f\n",x, Li_ion->interp_h(x)*1e3,Li_ion->interp_s(x)*1e3);
        //printf("\nMethods of Thermoclass: At x=%f, s=%f, scalc=%f, xlx=%f\n",x, GasConstant*Li_ion->mean_X(Li_ion->getEntropy_R()), Li_ion->interp_s(x)*1e3+GasConstant*log(x/(1-x)), GasConstant*Li_ion->sum_xlogx());
        delete tp;
    } catch (const std::exception& ex) {
        std::cout<<ex.what()<<std::endl;
    }
}

void printData(Cantera::ThermoPhase* tp) {
    double c[4];
    double minTemp, maxTemp, refPressure;
    double T=298.15;
    Cantera::vector_fp Xo(tp->nSpecies());
    Cantera::vector_fp G(tp->nSpecies());
    Cantera::vector_fp H(tp->nSpecies());
    Cantera::vector_fp S(tp->nSpecies());
    tp->getMoleFractions(Xo.data());
    tp->getGibbs_RT(G.data());
    tp->getEnthalpy_RT(H.data());
    tp->getEntropy_R(S.data());
    for (size_t k = 0; k < tp->nSpecies(); k++) {
        int type = tp->speciesThermo(k).reportType(k);
        tp->speciesThermo(k).reportParams(k,type,c,minTemp, maxTemp, refPressure);
        std::cout << tp->speciesName(k) <<" X= "<< Xo[k] <<" DH= " << c[1]<<" DS= " << c[2]<< std::endl;
        std::cout << tp->speciesName(k) <<" DH= "<< H[k]*Cantera::GasConstant*T <<" DS= " << S[k]*Cantera::GasConstant<<" DG= " << G[k]*Cantera::GasConstant*T<< std::endl;
    }
}

void thermoInitCathodePhases(std::string inFile) {
    std::string s = "";
    double T=298.15;
    double P=101325;
    double dG, dH, dS, dG0;
    try {
    	Cantera::ThermoPhase* tp = (Cantera::newPhase(inFile,"elyte"));
        std::cout<<tp->report()<<std::endl;
        /*ThermoPhase* tp1 = (newPhase(inFile,"cathode"));
        //ThermoPhase* tp2 = (newPhase(inFile,"conductor"));
        //std::vector<ThermoPhase*> phaseList;
        //phaseList.push_back(tp);
        //phaseList.push_back(tp1);
        phaseList.push_back(tp2);
        //Edge surf(inFile, "interface_cathode", phaseList);
        Interface surf(inFile, "C_surface", phaseList);

        //Interface surf(inFile, "interface_cathode", phaseList);
        //std::cout<<surf.phaseIndex("cathode")<<" "<<surf.phaseIndex("electrolyte")<<" "<<surf.phaseIndex("electron_cathode");
        tp->setState_TP(T,P);
        //tp1->setState_TP(T,P);
        tp2->setState_TP(T,P);
        surf.setState_TP(T,P);
        //tp->setMoleFractions(ns);
        for (size_t k = 0; k < phaseList.size(); k++) {
            phaseList[k]->setState_TP(T,P);
        }
        surf.setState_TP(T,P);
        surf.getDeltaGibbs(&dG);
        surf.getDeltaSSGibbs(&dG0);
        surf.getDeltaEnthalpy(&dH);
        surf.getDeltaEntropy(&dS);
        vector_fp kc(surf.nReactions());
        vector_fp fdot(surf.nReactions());
        vector_fp rdot(surf.nReactions());
        surf.getEquilibriumConstants(kc.data());
        surf.getFwdRateConstants(fdot.data());
        surf.getRevRateConstants(rdot.data());
        std::cout<<"No: "<<surf.nReactions()<<" "<<surf.reactionString(0)<<std::endl;
        std::cout<<"dG= "<<dG0/1000<<" dH= "<<dH/1000<<" dS= "<<dS/1000<<" Kc = " <<kc[0]<<std::endl;

        vector_fp wdot(surf.nTotalSpecies());
        surf.getNetProductionRates(wdot.data());

        //surf.setMultiplier(0,f);
        //for (size_t k = 0; k < phaseList.size(); k++) {
        //std::cout << "Name: "<<phaseList[k]->name()<<" density: "<<phaseList[k]->density()<< std::endl;
        for (size_t kk = 0; kk < surf.nTotalSpecies(); kk++)
            std::cout << surf.kineticsSpeciesName(kk)<< " wdot = " << wdot[kk]<< std::endl;

        // create a reaction path diagram builder
        ReactionPathBuilder b;
        std::ofstream rplog("rp1.log");   // log file
        std::ofstream rplot("rp1.dot");   // output file
        b.init(rplog, surf);         // initialize

        // main loop
        writeRxnPathDiagram(0, b, surf, rplog, rplot);
        for (int i = 1; i <= nsteps; i++) {
            tm = i*dt;
            sim.advance(tm);
            writeRxnPathDiagram(tm, b, gas, rplog, rplot);
        }

        //}*/
    } catch (Cantera::CanteraError& err){
        std::cout<<err.what()<< std::endl;
    }
}

void thermoInitAnodeROP_demo(std::string inFile) {
    std::string s = "";
    std::string phName, spName;
    double ns[2];
    double T=298.15;
    double P=101325;
    double x,dG[20],dH[20],dS[20],dG0[20],dH0[20],dS0[20];
    double minTemp, maxTemp, refPressure;
    double c[4];
    int type;

    Cantera::ThermoPhase* tp = (Cantera::newPhase(inFile,"anode"));
    Cantera::ThermoPhase* tp1 = (Cantera::newPhase(inFile,"electrolyte"));
    Cantera::ThermoPhase* tp2 = (Cantera::newPhase(inFile,"electron"));

    std::vector<Cantera::ThermoPhase*> phaseList;

    phaseList.push_back(tp);
    phaseList.push_back(tp1);
    phaseList.push_back(tp2);

    Cantera::Edge surf(inFile, "edge_anode_electrolyte", phaseList);

    //Interface surf(inFile, "edge_anode_electrolyte", phaseList);
    //Interface refSurf(inFile, "edge_lithium_electrolyte", refPhaseList);
    std::cout<<surf.phaseIndex("anode")<<" "<<surf.phaseIndex("electrolyte")<<" "<<surf.phaseIndex("electron_anode");
    tp->setState_TP(T,P);
    tp1->setState_TP(T,P);
    tp2->setState_TP(T,P);
    surf.setState_TP(T,P);

    Cantera::MultiSpeciesThermo& sp = tp->speciesThermo();
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
            std::cout<< "Phase: " <<  phName << " species: "<< spName << " dH0 = " << dH0[n]*Cantera::GasConstant*T << " dS0 = " << dS0[n]*Cantera::GasConstant <<std::endl;
            std::cout<< "Phase: " <<  phName << " species: "<< spName << " partial dH = " << dH[n] << " partial dS = " << dS[n] <<std::endl;
        }
    }
    Cantera::vector_fp wdot(tp->nSpecies() + tp1->nSpecies() +tp2->nSpecies() +surf.nSpecies());
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

    for (size_t k = 0; k < tp2->nSpecies(); k++) {
        std::cout << tp2->speciesName(k)<< " wdot = " << wdot[k+tp1->nSpecies()+tp->nSpecies()] <<" speciesIdx = " <<k+tp1->nSpecies()+tp->nSpecies() << std::endl;
    }
}

void thermoInitRefdemo(std::string inFile) {

    std::string s = "";
    std::string phName, spName;
    double T=298.15;
    double P=101325;
    double dG[20],dH[20],dS[20],dG0[20],dH0[20],dS0[20];

    std::vector<Cantera::ThermoPhase*> refPhaseList;

    Cantera::ThermoPhase* tp1 = (Cantera::newPhase(inFile,"electrolyte"));
    Cantera::ThermoPhase* tp2 = (Cantera::newPhase(inFile,"electron"));
    Cantera::ThermoPhase* tp = (Cantera::newPhase(inFile,"lithium_metal"));

    refPhaseList.push_back(tp);
    refPhaseList.push_back(tp1);
    refPhaseList.push_back(tp2);

    Cantera::Edge refSurf(inFile, "edge_reference_electrode", refPhaseList);

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
            std::cout<< "Phase: " <<  phName << " species: "<< spName << " dH0 = " << dH0[n]*Cantera::GasConstant*T << " dS0 = " << dS0[n]*Cantera::GasConstant <<std::endl;
            std::cout<< "Phase: " <<  phName << " species: "<< spName << " partial dH = " << dH[n] << " partial dS = " << dS[n] <<std::endl;
        }
    }
    Cantera::vector_fp wdot(tp->nSpecies() + tp1->nSpecies() +tp2->nSpecies() +refSurf.nSpecies());
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

void writeRxnPathDiagram(double time, Cantera::ReactionPathBuilder& b,
		Cantera::Kinetics& reaction, std::ostream& logfile, std::ostream& outfile)
{
    // create a new empty diagram
	Cantera::ReactionPathDiagram d;
    d.show_details = false; // show the details of which reactions contribute to the flux
    d.threshold = 0.001; // set the threshold for the minimum flux relative value
    d.bold_color = "orange"; // color for bold lines
    d.normal_color = "steelblue"; // color for normal-weight lines
    d.dashed_color = "gray"; // color for dashed lines
    d.dot_options = "center=1;size=\"6,9\";ratio=auto"; // options for the 'dot' program
    d.bold_min = 0.0; // minimum relative flux for bold lines
    d.dashed_max = 0.01; // maximum relative flux for dashed lines
    d.label_min = 0.01; // minimum relative flux for labels
    d.scale = -1; // autoscale
    d.flow_type = Cantera::OneWayFlow; //OneWayFlow; // set to either NetFlow or OneWayFlow
    d.arrow_width = -2.0; // arrow width. If < 0, then scale with flux value
    d.title = fmt::format("time = {} (s)", time); // title
    b.build(reaction, "E", logfile, d); // build the diagram following elemental nitrogen
    d.exportToDot(outfile); // write an input file for 'dot'
}

void thermoTestSPM(std::string inFile) {
    std::string s = "";
    double T=298.15;
    double P=101325;
    double dG, dH, dS, dG0;
    try {
    	Cantera::ThermoPhase* tp = (Cantera::newPhase(inFile,"anode"));
    	Cantera::ThermoPhase* tp1 = (Cantera::newPhase(inFile,"electron"));
    	Cantera::ThermoPhase* tp2 = (Cantera::newPhase(inFile,"electrolyte"));
    	Cantera::ThermoPhase* tp3 = (Cantera::newPhase(inFile,"cathode"));
        std::vector<Cantera::ThermoPhase*> phaseList;
        phaseList.push_back(tp);
        phaseList.push_back(tp1);
        phaseList.push_back(tp2);
        Cantera::Interface surf(inFile, "edge_anode_electrolyte", phaseList);

        //Interface surf(inFile, "interface_cathode", phaseList);
        //std::cout<<surf.phaseIndex("cathode")<<" "<<surf.phaseIndex("electrolyte")<<" "<<surf.phaseIndex("electron_cathode");
        tp3->setState_TP(T,P);
        surf.setState_TP(T,P);
        for (size_t k = 0; k < phaseList.size(); k++) {
            phaseList[k]->setState_TP(T,P);
        }
        std::cout<<tp->report()<<std::endl;
        std::cout<<tp1->report()<<std::endl;
        std::cout<<tp2->report()<<std::endl;
        std::cout<<tp3->report()<<std::endl;
        /*
        //std::cout<<tp2->molarDensity()<<" : "<<tp2->density()<<std::endl;
        for (size_t kk = 0; kk < tp->nSpecies(); kk++)
                    std::cout << tp->speciesName(kk)<<" : "<<tp->standardConcentration(kk)<<" : "<< tp->concentration(kk)<< std::endl;
        for (size_t kk = 0; kk < tp2->nSpecies(); kk++)
                            std::cout<< tp2->speciesName(kk)<<" : "<<tp2->standardConcentration(kk)<<" : "<<tp2->concentration(kk)<< std::endl;
        surf.setState_TP(T,P);
        surf.getDeltaGibbs(&dG);
        surf.getDeltaSSGibbs(&dG0);
        surf.getDeltaEnthalpy(&dH);
        surf.getDeltaEntropy(&dS);
        vector_fp kc(surf.nReactions());
        vector_fp fdot(surf.nReactions());
        vector_fp rdot(surf.nReactions());
        surf.getEquilibriumConstants(kc.data());
        surf.getFwdRateConstants(fdot.data());
        surf.getRevRateConstants(rdot.data());
        std::cout<<"No: "<<surf.nReactions()<<" "<<surf.reactionString(0)<<std::endl;
        std::cout<<"dG= "<<dG0/1000<<" dH= "<<dH/1000<<" dS= "<<dS/1000<<" Kc = " <<kc[0]<<std::endl;
        std::cout<<" rf= "<<fdot[0]<<" rr= "<<rdot[0]<<std::endl;
        surf.getFwdRatesOfProgress(fdot.data());
        surf.getRevRatesOfProgress(rdot.data());
        std::cout<<" ropf= "<<fdot[0]<<" ropr= "<<rdot[0]<<std::endl;

        vector_fp wdot(surf.nTotalSpecies());
        surf.getNetProductionRates(wdot.data());

        //for (size_t k = 0; k < phaseList.size(); k++) {
        //std::cout << "Name: "<<phaseList[k]->name()<<" density: "<<phaseList[k]->density()<< std::endl;
        for (size_t kk = 0; kk < surf.nTotalSpecies(); kk++)
            std::cout << surf.kineticsSpeciesName(kk)<< " wdot = " << wdot[kk]<< std::endl;
            */

    } catch (Cantera::CanteraError& err){
        std::cout<<err.what()<< std::endl;
    }
}

void printThermoReactionTDPA(std::string inFile) {
    std::string surfName;
    double T=298.15;
    double P=101325;
    try {
        /*surfName = "O_surface";
        ThermoPhase* tp2 = (newPhase(inFile,"gas_cathode"));*/
        surfName = "C_surface";
        Cantera::ThermoPhase* tp2 = (Cantera::newPhase(inFile,"conductor"));
        /*surfName = "Li2O2_precipitation_from_solution";
        ThermoPhase* tp2 = (newPhase(inFile,"Li2O2"));*/
        /*surfName = "LiO2_precipitation";
        ThermoPhase* tp2 = (newPhase(inFile,"conductor"));
        ThermoPhase* tp3 = (newPhase(inFile,"LiO2"));
        ThermoPhase* tp4 = (newPhase(inFile,"Li2O2"));*/
        Cantera::ThermoPhase* tp1 = (Cantera::newPhase(inFile,"elyte"));
        std::vector<Cantera::ThermoPhase*> phaseList;
        phaseList.push_back(tp1);
        phaseList.push_back(tp2);
        // Set phase standard state
        for (size_t k = 0; k < phaseList.size(); k++) {
            phaseList[k]->setState_TP(T,P);
        }
        printThermoKinetics(inFile, surfName, phaseList, T, P);
    } catch (Cantera::CanteraError& err){
        std::cout<<err.what()<< std::endl;
    }
}

void printThermoReactionFerrocene(std::string inFile) {
    std::string surfName;
    std::ifstream dataFile;
    std::string line;
    std::vector<double> data(5);

    double T=298.15;
    double P=101325;
    dataFile.open("data/outData.dat");
	FILE * pFile;
	pFile = fopen ("compData.dat","w");
    try
    {
        surfName = "WE_surface";
        Cantera::ThermoPhase* tp1 = (Cantera::newPhase(inFile,"electrolyte"));
        Cantera::ThermoPhase* tp2 = (Cantera::newPhase(inFile,"conductor"));
        std::vector<Cantera::ThermoPhase*> phaseList;
        phaseList.push_back(tp1);
        phaseList.push_back(tp2);
        // Set phase standard state
        for (size_t k = 0; k < phaseList.size(); k++)
        {
            phaseList[k]->setState_TP(T,P);
        }
        printThermoKinetics(inFile, surfName, phaseList, T, P);
        //
        while (std::getline(dataFile, line))
        {
        	double d0, d1, d2, d3, d4, rop[3];
        	std::sscanf(line.c_str(),"%lf\t%lf\t%lf\t%lf\t%lf", &d0, &d1, &d2, &d3, &d4);
        	calc_itotCantera2s(d1, d2, d0, rop, d3);
        	fprintf (pFile, "%.2e\t%12.3e\t%12.3e\t%12.3e\t%12.3e\n",d0, d1, d2, d3, rop[2]*96485);
        }
		/* Free memory */
		dataFile.close();
		fclose (pFile);
    }
    catch (Cantera::CanteraError& err)
    {
        std::cout<<err.what()<< std::endl;
    }
}

void printThermoReactionLi2O2(std::string inFile) {
    //std::string surfName[] = {"Li2O2_surface_RM", "Li2O2_precipitation_from_solution",
    		//"Li2O2_precipitation_surface", "LiO2_precipitation", "Li2O2_surface"};
    std::string surfName[] = {"C_surface", "O_surface", "Li2O2_surface"};
    double T=298.15;
    double P=101325;
	std::vector<Cantera::ThermoPhase*> phaseList;
    try {
    	for (size_t i=0; i<surfName->size(); i++)
    	{
    		if (surfName[i] == "O_surface")
    		{
    			std::cout<<"Printing "<<surfName[i]<<std::endl;
    			phaseList.clear();
    			Cantera::ThermoPhase* tp1 = (Cantera::newPhase(inFile,"gas_cathode"));
    			Cantera::ThermoPhase* tp2 = (Cantera::newPhase(inFile,"elyte"));
    			phaseList.push_back(tp1);
    			phaseList.push_back(tp2);
    		}
    		if (surfName[i] == "Li2O2_surface_RM")
    		{
    			std::cout<<"Printing "<<surfName[i]<<std::endl;
    			phaseList.clear();
    			Cantera::ThermoPhase* tp1 = (Cantera::newPhase(inFile,"elyte"));
    			Cantera::ThermoPhase* tp2 = (Cantera::newPhase(inFile,"Li2O2"));
    			Cantera::ThermoPhase* tp3 = (Cantera::newPhase(inFile,"conductor"));
    			phaseList.push_back(tp1);
    			phaseList.push_back(tp2);
    			phaseList.push_back(tp3);
    		}
    		if (surfName[i] == "Li2O2_precipitation_from_solution")
    		{
    			std::cout<<"Printing "<<surfName[i]<<std::endl;
    			phaseList.clear();
    			Cantera::ThermoPhase* tp1 = (Cantera::newPhase(inFile,"elyte"));
    			Cantera::ThermoPhase* tp2 = (Cantera::newPhase(inFile,"Li2O2"));
    			phaseList.push_back(tp1);
    			phaseList.push_back(tp2);
    		}
    		if (surfName[i] == "Li2O2_precipitation_surface")
    		{
    			std::cout<<"Printing "<<surfName[i]<<std::endl;
				phaseList.clear();
				Cantera::ThermoPhase* tp1 = (Cantera::newPhase(inFile,"elyte"));
				Cantera::ThermoPhase* tp2 = (Cantera::newPhase(inFile,"LiO2"));
				Cantera::ThermoPhase* tp3 = (Cantera::newPhase(inFile,"Li2O2"));
    			phaseList.push_back(tp1);
    			phaseList.push_back(tp2);
    			phaseList.push_back(tp3);
    		}
    		if (surfName[i] == "LiO2_precipitation")
    		{
    			std::cout<<"Printing "<<surfName[i]<<std::endl;
				phaseList.clear();
				Cantera::ThermoPhase* tp1 = (Cantera::newPhase(inFile,"elyte"));
				Cantera::ThermoPhase* tp2 = (Cantera::newPhase(inFile,"LiO2"));
				Cantera::ThermoPhase* tp3 = (Cantera::newPhase(inFile,"conductor"));
    			phaseList.push_back(tp1);
    			phaseList.push_back(tp2);
    			phaseList.push_back(tp3);
    		}
    		if (surfName[i] == "Li2O2_surface")
    		{
    			std::cout<<"Printing "<<surfName[i]<<std::endl;
				phaseList.clear();
				Cantera::ThermoPhase* tp1 = (Cantera::newPhase(inFile,"elyte"));
				Cantera::ThermoPhase* tp2 = (Cantera::newPhase(inFile,"Li2O2"));
    			Cantera::ThermoPhase* tp3 = (Cantera::newPhase(inFile,"conductor"));
    			phaseList.push_back(tp1);
    			phaseList.push_back(tp2);
    			phaseList.push_back(tp3);
    		}
    		else
    			continue;
    		// Set phase standard state
    		std::cout<<"Setting phase standard state... "<<std::endl;
    		for (size_t k = 0; k < phaseList.size(); k++) {
    			phaseList[k]->setState_TP(T,P);
    		}
    		printThermoKinetics(inFile, surfName[i], phaseList, T, P);
    	}

    } catch (Cantera::CanteraError& err){
        std::cout<<err.what()<< std::endl;
    }
}

void printThermoReactionO2(std::string inFile) {
    std::string surfName;
    double T=298.15;
    double P=101325;
    try {
        surfName = "O_surface";
        Cantera::ThermoPhase* tp2 = (Cantera::newPhase(inFile,"gas_cathode"));
        Cantera::ThermoPhase* tp1 = (Cantera::newPhase(inFile,"elyte"));
        std::vector<Cantera::ThermoPhase*> phaseList;
        phaseList.push_back(tp1);
        phaseList.push_back(tp2);
        Cantera::Interface surf(inFile, surfName, phaseList);
        // Set phase standard state
        for (size_t k = 0; k < phaseList.size(); k++) {
            phaseList[k]->setState_TP(T,P);
        }
        surf.setState_TP(T,P);
        printThermoKinetics(inFile, surfName, phaseList, T, P);
    } catch (Cantera::CanteraError& err){
        std::cout<<err.what()<< std::endl;
    }
}

void printThermoReactionLi(std::string inFile) {
    std::string surfName;
    double T=298.15;
    double P=101325;
    try {
        surfName = "Li_surface";
        Cantera::ThermoPhase* tp1 = (Cantera::newPhase(inFile,"lithium"));
        Cantera::ThermoPhase* tp2 = (Cantera::newPhase(inFile,"elyte"));
        Cantera::ThermoPhase* tp3 = (Cantera::newPhase(inFile,"conductor"));
        std::vector<Cantera::ThermoPhase*> phaseList;
        phaseList.push_back(tp1);
        phaseList.push_back(tp2);
        phaseList.push_back(tp3);
        Cantera::Interface surf(inFile, surfName, phaseList);
        // Set phase standard state
        for (size_t k = 0; k < phaseList.size(); k++) {
            phaseList[k]->setState_TP(T,P);
        }
        std::cout<<tp3->speciesIndex("electron")<<std::endl;
        printThermoKinetics(inFile, surfName, phaseList, T, P);
    } catch (Cantera::CanteraError& err){
        std::cout<<err.what()<< std::endl;
    }
}

void printThermoReactionPEM(std::string inFile) {
    std::string surfName;
    double T=343.15;
    double P=101325;
    try {
        surfName = "cathode_reaction_surface";
        Cantera::ThermoPhase* tp1 = (Cantera::newPhase(inFile,"gas_cathode"));
        Cantera::ThermoPhase* tp2 = (Cantera::newPhase(inFile,"Platinum"));
        Cantera::ThermoPhase* tp3 = (Cantera::newPhase(inFile,"Nafion"));
        std::vector<Cantera::ThermoPhase*> phaseList;
        phaseList.push_back(tp1);
        phaseList.push_back(tp2);
        phaseList.push_back(tp3);
        Cantera::Interface surf(inFile, surfName, phaseList);
        // Set phase standard state
        for (size_t k = 0; k < phaseList.size(); k++) {
            phaseList[k]->setState_TP(T,P);
        }
        printThermoKinetics(inFile, surfName, phaseList, T, P);

    } catch (Cantera::CanteraError& err){
        std::cout<<err.what()<< std::endl;
    }
}

void printThermoKinetics(std::string inFile, std::string surfName, std::vector<Cantera::ThermoPhase*> phaseList, double T, double P) {

	Cantera::Interface surf(inFile, surfName, phaseList);
    surf.setState_TP(T,P);
    // Phase thermodynamics
    std::cout<<"\nPrinting species thermodynamics..."<<std::endl;
    for (size_t k = 0; k < surf.nPhases(); k++) {
        std::cout << "Name: "<<surf.thermo(k).name()<<" #Species: "<< surf.thermo(k).nSpecies()
        		<<" density: "<<surf.thermo(k).density()
				<<" Epot: "<<surf.thermo(k).electricPotential() <<std::endl;
        Cantera::vector_fp mu0(surf.thermo(k).nSpecies());
        Cantera::vector_fp muE(surf.thermo(k).nSpecies());
        Cantera::vector_fp mu(surf.thermo(k).nSpecies());
        Cantera::vector_fp h0(surf.thermo(k).nSpecies());
        Cantera::vector_fp s0(surf.thermo(k).nSpecies());
        Cantera::vector_fp ac(surf.thermo(k).nSpecies());
        surf.thermo(k).getStandardChemPotentials(mu0.data());
        surf.thermo(k).getElectrochemPotentials(muE.data());
        surf.thermo(k).getChemPotentials(mu.data());
        surf.thermo(k).getEnthalpy_RT(h0.data());
        surf.thermo(k).getEntropy_R(s0.data());
        surf.thermo(k).getActivityConcentrations(ac.data());
        surf.thermo(1).setElectricPotential(0);
        double mmu0;
        for (size_t kk = 0; kk < surf.thermo(k).nSpecies(); kk++) {
            mmu0 = mu0[kk] + Cantera::Faraday * surf.thermo(k).electricPotential()*surf.thermo(k).charge(kk);
            mmu0 -= surf.thermo(0).RT() * surf.thermo(k).logStandardConc(kk);
            std::cout<<std::setprecision (15)<<"Species: "<< surf.thermo(k).speciesName(kk)
                     <<", c0 (kmol/m3) = "<<surf.thermo(k).standardConcentration(kk)<<", c (kmol/m3) = "<<surf.thermo(k).concentration(kk)
                     <<", ac (kmol/m3) = "<< ac[kk]<<", mu (J/mol) = "<<mmu0/1e3
                     <<", mu0 (J/mol) = "<<mu0[kk]/1e3<<", muE (J/mol) = "<<muE[kk]/1e3<<" h0 (J/mol): "<< h0[kk]*Cantera::GasConstant*T/1e3
                     <<" s0 (J/mol/K): "<<s0[kk]*Cantera::GasConstant/1e3<<std::endl;
        }
        std::cout<<std::endl;
    }

    // Reaction thermodynamics
    std::cout<<surf.report();
    std::cout<<"Printing reaction thermodynamics..."<<std::endl;
    Cantera::vector_fp dG0(surf.nReactions()), dG(surf.nReactions());
    Cantera::vector_fp dmu(surf.nReactions());
    /*vector_fp dH(surf.nReactions());
    vector_fp dS(surf.nReactions());*/
    Cantera::vector_fp kc(surf.nReactions()), fdot(surf.nReactions());
    Cantera::vector_fp rdot(surf.nReactions()), frop(surf.nReactions()), rrop(surf.nReactions());
    /*surf.getDeltaSSGibbs(dG0.data());
    surf.getDeltaGibbs(dG.data());
    surf.getDeltaElectrochemPotentials(dmu.data());
    surf.getDeltaEnthalpy(dH.data());
    surf.getDeltaEntropy(dS.data());*/
    surf.getEquilibriumConstants(kc.data());
    surf.getFwdRateConstants(fdot.data());
    surf.getRevRateConstants(rdot.data());
    surf.getFwdRatesOfProgress(frop.data());
    surf.getRevRatesOfProgress(rrop.data());
    for (size_t k = 0; k < surf.nReactions(); k++) {
        std::cout<<std::setprecision (15)<<"No: "<<k<<" "<<surf.reactionString(k)<<" Type:"<<surf.reactionType(k)<<"\n\n"
        /*<<"dG(J/mol)= "<<dG[k]/1e3<<" dG0(J/mol)= "<<dG0[k]/1e3<<" dmu(J/mol)= "<<dmu[k]/1e3<<"\n"
        <<"dH(J/mol)= "<<dH[k]/1e3<<" dS(J/(mol.K))= "<<dS[k]/1e3<<"\n"
        <<"G_cal(J/mol)=dH-TdS= " <<1e3*(dH[k]-T*dS[k])<<"\n"
        <<"Keq = " <<kc[k]<<" Keq_cal = exp(-dG0/RT) = "<<std::exp(-dG0[k]/(GasConstant*T))
        <<" Keq_cal = exp(-dG/RT) = "<<std::exp(-dG[k]/(GasConstant*T))<<"\n"
		<<" Keq_cal = exp(-dmu/RT) = "<<std::exp(-dmu[k]/(GasConstant*T))<<"\n"*/
        <<"kf = "<< fdot[k] << " frop = " << frop[k] << " kr = " << rdot[k] << " rrop = " << rrop[k] << " Keq = kf/kr = " << fdot[k]/rdot[k]<< "\n\n";
    }
    Cantera::vector_fp wdot(surf.nTotalSpecies());
    surf.getNetProductionRates(wdot.data());
    std::cout<<"Printing species ROP..."<<std::endl;
    for (size_t kk = 0; kk < surf.nTotalSpecies(); kk++)
        std::cout << surf.kineticsSpeciesName(kk)<< " wdot = " << wdot[surf.kineticsSpeciesIndex(surf.kineticsSpeciesName(kk))] <<std::endl;

}
