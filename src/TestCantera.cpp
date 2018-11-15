//============================================================================
// Name        : TestCantera.cpp
// Author      : Manik Mayur
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <cantera/thermo/BinarySolutionTabulatedThermo.h>
#include <iostream>
#include <string>
#include <sstream>
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/Interface.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/kinetics/importKinetics.h"
#include "cantera/Edge.h"
#include "cantera/reactionpaths.h"


using namespace Cantera;

void thermoInitAnode_demo(std::string inFile);
void thermoInitCathodePhases(std::string inFile);
void thermoInitAnodeROP_demo(std::string inFile);
void thermoInitRefdemo(std::string inFile);
void writeRxnPathDiagram(double time, ReactionPathBuilder& b,
        Kinetics& reaction, std::ostream& logfile, std::ostream& outfile);
void thermoTestSPM(std::string inFile);
void printThermoReactionO2(std::string inFile);
void printThermoReactionTDPA(std::string inFile);
void printThermoReactionLi(std::string inFile);

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
        //printThermoReactionO2("cantera\/Work_LiO2_organic_LiBaLu_CV_TDPA.xml");
        //printThermoReactionTDPA("cantera\/Work_LiO2_organic_LiBaLu_CV_TDPA.xml");
        printThermoReactionLi("cantera\/Work_LiO2_organic_LiBaLu_CV_TDPA.xml");
        //thermoInitElectrolyte_demo(inFile);
        //simple_demo2();
    }
    catch (CanteraError& err){
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
        ThermoPhase* tp = (newPhase(inFile,"anode"));
        BinarySolutionTabulatedThermo* Li_ion = dynamic_cast<BinarySolutionTabulatedThermo*> (tp);
        Li_ion->setState_TP(T,P);

        MultiSpeciesThermo& sp = Li_ion->speciesThermo();
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

void printData(ThermoPhase* tp) {
    double c[4];
    double minTemp, maxTemp, refPressure;
    double T=298.15;
    vector_fp Xo(tp->nSpecies());
    vector_fp G(tp->nSpecies());
    vector_fp H(tp->nSpecies());
    vector_fp S(tp->nSpecies());
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

void thermoInitCathodePhases(std::string inFile) {
    std::string s = "";
    double T=298.15;
    double P=101325;
    double dG, dH, dS, dG0;
    try {
        ThermoPhase* tp = (newPhase(inFile,"elyte"));
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
    } catch (CanteraError& err){
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

    ThermoPhase* tp = (newPhase(inFile,"anode"));
    ThermoPhase* tp1 = (newPhase(inFile,"electrolyte"));
    ThermoPhase* tp2 = (newPhase(inFile,"electron"));

    std::vector<ThermoPhase*> phaseList;

    phaseList.push_back(tp);
    phaseList.push_back(tp1);
    phaseList.push_back(tp2);

    Edge surf(inFile, "edge_anode_electrolyte", phaseList);

    //Interface surf(inFile, "edge_anode_electrolyte", phaseList);
    //Interface refSurf(inFile, "edge_lithium_electrolyte", refPhaseList);
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

    std::vector<ThermoPhase*> refPhaseList;

    ThermoPhase* tp1 = (newPhase(inFile,"electrolyte"));
    ThermoPhase* tp2 = (newPhase(inFile,"electron"));
    ThermoPhase* tp = (newPhase(inFile,"lithium_metal"));

    refPhaseList.push_back(tp);
    refPhaseList.push_back(tp1);
    refPhaseList.push_back(tp2);

    Edge refSurf(inFile, "edge_reference_electrode", refPhaseList);

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

void writeRxnPathDiagram(double time, ReactionPathBuilder& b,
        Kinetics& reaction, std::ostream& logfile, std::ostream& outfile)
{
    // create a new empty diagram
    ReactionPathDiagram d;
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
    d.flow_type = OneWayFlow; //OneWayFlow; // set to either NetFlow or OneWayFlow
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
        ThermoPhase* tp = (newPhase(inFile,"anode"));
        ThermoPhase* tp1 = (newPhase(inFile,"electron"));
        ThermoPhase* tp2 = (newPhase(inFile,"electrolyte"));
        ThermoPhase* tp3 = (newPhase(inFile,"cathode"));
        std::vector<ThermoPhase*> phaseList;
        phaseList.push_back(tp);
        phaseList.push_back(tp1);
        phaseList.push_back(tp2);
        Interface surf(inFile, "edge_anode_electrolyte", phaseList);

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

    } catch (CanteraError& err){
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
        ThermoPhase* tp2 = (newPhase(inFile,"conductor"));
        /*surfName = "Li2O2_precipitation_from_solution";
        ThermoPhase* tp2 = (newPhase(inFile,"Li2O2"));*/
        /*surfName = "LiO2_precipitation";
        ThermoPhase* tp2 = (newPhase(inFile,"conductor"));
        ThermoPhase* tp3 = (newPhase(inFile,"LiO2"));
        ThermoPhase* tp4 = (newPhase(inFile,"Li2O2"));*/
        ThermoPhase* tp1 = (newPhase(inFile,"elyte"));
        std::vector<ThermoPhase*> phaseList;
        phaseList.push_back(tp1);
        //phaseList.push_back(tp3);
        //phaseList.push_back(tp4);
        phaseList.push_back(tp2);
        Interface surf(inFile, surfName, phaseList);
        // Set phase standard state
        for (size_t k = 0; k < phaseList.size(); k++) {
            phaseList[k]->setState_TP(T,P);
        }
        surf.setState_TP(T,P);
        double mmu0;
        std::cout<<"\nPrinting species thermodynamics..."<<std::endl;
        for (size_t k = 0; k < surf.nPhases(); k++) {
            std::cout << "Name: "<<surf.thermo(k).name()<<" #Species: "<< surf.thermo(k).nSpecies()
            		<<" density: "<<surf.thermo(k).density()
					<<" Epot: "<<surf.thermo(k).electricPotential()
					<<std::endl;
            vector_fp mu0(surf.thermo(k).nSpecies());
            vector_fp mu(surf.thermo(k).nSpecies());
            vector_fp h0(surf.thermo(k).nSpecies());
            vector_fp s0(surf.thermo(k).nSpecies());
            vector_fp ac(surf.thermo(k).nSpecies());

            surf.thermo(k).getStandardChemPotentials(mu0.data());
            surf.thermo(k).getChemPotentials(mu.data());
            surf.thermo(k).getActivityConcentrations(ac.data());

            for (size_t kk = 0; kk < surf.thermo(k).nSpecies(); kk++) {
                mmu0 = mu0[kk] + Faraday * surf.thermo(k).electricPotential()*surf.thermo(k).charge(kk);
                //mmu0 -= surf.thermo(0).RT() * surf.thermo(k).logStandardConc(k);
                std::cout <<"Species: "<< surf.thermo(k).speciesName(kk)
                        <<" : c0 = "<<surf.thermo(k).standardConcentration(kk)<<" : c = "<<surf.thermo(k).concentration(kk)
                        <<" : ac = "<< ac[kk]
                        <<" mmu0: "<<mmu0<<" mu0: "<<mu0[kk]<<" mu: "<<mu[kk]<<" h0: "<<h0[kk]*GasConstant*T
                        <<" s0: "<<s0[kk]*GasConstant<<" mu0_calc=h0-T*s0: " <<h0[kk]*GasConstant*T-T*s0[kk]*GasConstant<<std::endl;
            }
            std::cout<<std::endl;
        }
        // Reaction thermodynamics
        std::cout<<"Printing reaction thermodynamics..."<<std::endl;
        vector_fp dG0(surf.nReactions());
        vector_fp dG(surf.nReactions());
        vector_fp dH(surf.nReactions());
        vector_fp dS(surf.nReactions());
        vector_fp kc(surf.nReactions());
        vector_fp fdot(surf.nReactions());
        vector_fp rdot(surf.nReactions());
        vector_fp frop(surf.nReactions());
        vector_fp rrop(surf.nReactions());
        surf.getDeltaSSGibbs(dG0.data());
        surf.getDeltaGibbs(dG.data());
        //surf.getDeltaEnthalpy(dH.data());
        //surf.getDeltaEntropy(dS.data());
        surf.getEquilibriumConstants(kc.data());
        surf.getFwdRateConstants(fdot.data());
        surf.getRevRateConstants(rdot.data());
        surf.getFwdRatesOfProgress(frop.data());
        surf.getRevRatesOfProgress(rrop.data());
        for (size_t k = 0; k < surf.nReactions(); k++) {
            std::cout<<"No: "<<k<<" "<<surf.reactionString(k)<<std::endl;
            std::cout<<"dG(J/kmol)= "<<dG[k]<<" dG0(J/kmol)= "<<dG0[k]<<std::endl;
            std::cout<<"dH(J/kmol)= "<<dH[k]<<" dS(J/(kmol.K))= "<<dS[k]<<std::endl;
            std::cout<<"G_cal(J/kmol)=dH-TdS= " <<dH[k]-T*dS[k]<<std::endl;
            std::cout<<"Keq = " <<kc[k]<<" K_cal=c0P/c0R*exp(-dG/RT)= "<<tp2->standardConcentration(0)/tp1->standardConcentration(0)*(std::exp(-dG[k]/(GasConstant*T)));
            std::cout<<"K_cal=c0P/c0R*exp(-dG0/RT)= "<<tp2->standardConcentration(0)/tp1->standardConcentration(0)*(std::exp(-dG0[k]/(GasConstant*T)))<<std::endl;
            std::cout<<"fdot = "<< fdot[k] << " frop = " << frop[k] << " rdot = " << rdot[k] << " rrop = " << rrop[k] << " Keq = " << fdot[k]/rdot[k]<< std::endl;
        }
        vector_fp wdot(surf.nTotalSpecies());
        surf.getNetProductionRates(wdot.data());
        std::cout<<"Printing species ROP..."<<std::endl;
        for (size_t kk = 0; kk < surf.nTotalSpecies(); kk++)
            std::cout << surf.kineticsSpeciesName(kk)<< " wdot = " << wdot[kk]<<std::endl;
    } catch (CanteraError& err){
        std::cout<<err.what()<< std::endl;
    }
}

void printThermoReactionO2(std::string inFile) {
    std::string surfName;
    double T=298.15;
    double P=101325;
    try {
        surfName = "O_surface";
        ThermoPhase* tp2 = (newPhase(inFile,"gas_cathode"));
        ThermoPhase* tp1 = (newPhase(inFile,"elyte"));
        std::vector<ThermoPhase*> phaseList;
        phaseList.push_back(tp1);
        phaseList.push_back(tp2);
        Interface surf(inFile, surfName, phaseList);
        // Set phase standard state
        for (size_t k = 0; k < phaseList.size(); k++) {
            phaseList[k]->setState_TP(T,P);
        }
        surf.setState_TP(T,P);
        vector_fp c(phaseList[0]->nSpecies());
        vector_fp x(phaseList[0]->nSpecies());
        //tp2->setPressure(P*1.958);

        for (size_t k = 0; k < phaseList[0]->nSpecies(); k++) {
                    std::cout << "Name: "<<phaseList[0]->speciesName(k)<<" x: "<< x[k]<<std::endl;
        }
        double mmu0;
        std::cout<<"\nPrinting species thermodynamics..."<<std::endl;
        for (size_t k = 0; k < surf.nPhases(); k++) {
            std::cout << "Name: "<<surf.thermo(k).name()<<" #Species: "<< surf.thermo(k).nSpecies()<<" density: "<<surf.thermo(k).density()<<std::endl;
            vector_fp mu0(surf.thermo(k).nSpecies());
            vector_fp mu(surf.thermo(k).nSpecies());
            vector_fp h0(surf.thermo(k).nSpecies());
            vector_fp s0(surf.thermo(k).nSpecies());
            vector_fp ac(surf.thermo(k).nSpecies());
            surf.thermo(k).getStandardChemPotentials(mu0.data());
            surf.thermo(k).getChemPotentials(mu.data());
            surf.thermo(k).getEnthalpy_RT(h0.data());
            surf.thermo(k).getEntropy_R(s0.data());
            surf.thermo(k).getActivityConcentrations(ac.data());
            for (size_t kk = 0; kk < surf.thermo(k).nSpecies(); kk++) {
                mmu0 = mu0[kk] + Faraday * surf.thermo(k).electricPotential()*surf.thermo(k).charge(kk);
                mmu0 -= surf.thermo(0).RT() * surf.thermo(k).logStandardConc(k);
                std::cout <<"Species: "<< surf.thermo(k).speciesName(kk)
                        <<" : c0 = "<<surf.thermo(k).standardConcentration(kk)<<" : c = "<<surf.thermo(k).concentration(kk)
                        <<" : ac = "<< ac[kk]
                        <<" mmu0: "<<mmu0<<" mu0: "<<mu0[kk]<<" mu: "<<mu[kk]<<" h0: "<<h0[kk]*GasConstant*T
                        <<" s0: "<<s0[kk]*GasConstant<<" mu0_calc=h0-T*s0: " <<h0[kk]*GasConstant*T-T*s0[kk]*GasConstant<<std::endl;
            }
            std::cout<<std::endl;
        }
        // Reaction thermodynamics
        std::cout<<"Printing reaction thermodynamics..."<<std::endl;
        vector_fp dG0(surf.nReactions());
        vector_fp dG(surf.nReactions());
        vector_fp dH(surf.nReactions());
        vector_fp dS(surf.nReactions());
        vector_fp kc(surf.nReactions());
        vector_fp fdot(surf.nReactions());
        vector_fp rdot(surf.nReactions());
        vector_fp frop(surf.nReactions());
        vector_fp rrop(surf.nReactions());
        surf.getDeltaSSGibbs(dG0.data());
        surf.getDeltaGibbs(dG.data());
        //surf.getDeltaEnthalpy(dH.data());
        //surf.getDeltaEntropy(dS.data());
        surf.getEquilibriumConstants(kc.data());
        surf.getFwdRateConstants(fdot.data());
        surf.getRevRateConstants(rdot.data());
        surf.getFwdRatesOfProgress(frop.data());
        surf.getRevRatesOfProgress(rrop.data());
        for (size_t k = 0; k < surf.nReactions(); k++) {
            std::cout<<"No: "<<k<<" "<<surf.reactionString(k)<<std::endl;
            std::cout<<"dG(J/kmol)= "<<dG[k]<<" dG0(J/kmol)= "<<dG0[k]<<std::endl;
            std::cout<<"dH(J/kmol)= "<<dH[k]<<" dS(J/(kmol.K))= "<<dS[k]<<std::endl;
            std::cout<<"G_cal(J/kmol)=dH-TdS= " <<dH[k]-T*dS[k]<<std::endl;
            std::cout<<"Keq = " <<kc[k]<<" K_cal=c0P/c0R*exp(-dG/RT)= "<<tp2->standardConcentration(0)/tp1->standardConcentration(0)*(std::exp(-dG[k]/(GasConstant*T)));
            std::cout<<" K_cal=c0P/c0R*exp(-dG0/RT)= "<<tp2->standardConcentration(0)/tp1->standardConcentration(0)*(std::exp(-dG0[k]/(GasConstant*T)))<<std::endl;
            std::cout<<std::endl;
        }
        vector_fp wdot(surf.nTotalSpecies());
        surf.getNetProductionRates(wdot.data());
        std::cout<<"Printing species ROP..."<<std::endl;
        for (size_t kk = 0; kk < surf.nTotalSpecies(); kk++)
            std::cout << surf.kineticsSpeciesName(kk)<< " wdot = " << wdot[kk] << " fdot ="<< fdot[kk] << " frop = " << frop[kk] << " rdot = " << rdot[kk] << " rrop = " << rrop[kk] << " Keq = " << fdot[kk]/rdot[kk]<< std::endl;
    } catch (CanteraError& err){
        std::cout<<err.what()<< std::endl;
    }
}

void printThermoReactionLi(std::string inFile) {
    std::string surfName;
    double T=298.15;
    double P=101325;
    try {
        surfName = "Li_surface";
        ThermoPhase* tp1 = (newPhase(inFile,"lithium"));
        ThermoPhase* tp2 = (newPhase(inFile,"elyte"));
        ThermoPhase* tp3 = (newPhase(inFile,"conductor"));
        std::vector<ThermoPhase*> phaseList;
        phaseList.push_back(tp1);
        phaseList.push_back(tp2);
        phaseList.push_back(tp3);
        Interface surf(inFile, surfName, phaseList);
        // Set phase standard state
        for (size_t k = 0; k < phaseList.size(); k++) {
            phaseList[k]->setState_TP(T,P);
        }
        surf.setState_TP(T,P);

        double mmu0;
        std::cout<<"\nPrinting species thermodynamics..."<<std::endl;
        for (size_t k = 0; k < surf.nPhases(); k++) {
            std::cout << "Name: "<<surf.thermo(k).name()<<" #Species: "<< surf.thermo(k).nSpecies()
            		<<" density: "<<surf.thermo(k).density()
					<<" Epot: "<<surf.thermo(k).electricPotential() <<std::endl;
            vector_fp mu0(surf.thermo(k).nSpecies());
            vector_fp mu(surf.thermo(k).nSpecies());
            vector_fp h0(surf.thermo(k).nSpecies());
            vector_fp s0(surf.thermo(k).nSpecies());
            vector_fp ac(surf.thermo(k).nSpecies());
            surf.thermo(k).getStandardChemPotentials(mu0.data());
            surf.thermo(k).getChemPotentials(mu.data());
            surf.thermo(k).getEnthalpy_RT(h0.data());
            surf.thermo(k).getEntropy_R(s0.data());
            surf.thermo(k).getActivityConcentrations(ac.data());
            surf.thermo(1).setElectricPotential(1);

            for (size_t kk = 0; kk < surf.thermo(k).nSpecies(); kk++) {
                mmu0 = mu0[kk] + Faraday * surf.thermo(k).electricPotential()*surf.thermo(k).charge(kk);
                mmu0 -= surf.thermo(0).RT() * surf.thermo(k).logStandardConc(k);
                std::cout <<"Species: "<< surf.thermo(k).speciesName(kk)
                        <<" : c0 = "<<surf.thermo(k).standardConcentration(kk)<<" : c = "<<surf.thermo(k).concentration(kk)
                        <<" : ac = "<< ac[kk]
                        <<" mmu0: "<<mmu0<<" mu0: "<<mu0[kk]<<" mu: "<<mu[kk]<<" h0: "<<h0[kk]*GasConstant*T
                        <<" s0: "<<s0[kk]*GasConstant<<" mu0_calc=h0-T*s0: " <<h0[kk]*GasConstant*T-T*s0[kk]*GasConstant<<std::endl;
            }
            std::cout<<std::endl;
        }
        // Reaction thermodynamics
        std::cout<<"Printing reaction thermodynamics..."<<std::endl;
        vector_fp dG0(surf.nReactions());
        vector_fp dG(surf.nReactions());
        vector_fp dH(surf.nReactions());
        vector_fp dS(surf.nReactions());
        vector_fp kc(surf.nReactions());
        vector_fp fdot(surf.nReactions());
        vector_fp rdot(surf.nReactions());
        vector_fp frop(surf.nReactions());
        vector_fp rrop(surf.nReactions());
        surf.getDeltaSSGibbs(dG0.data());
        surf.getDeltaGibbs(dG.data());
        surf.getEquilibriumConstants(kc.data());
        surf.getFwdRateConstants(fdot.data());
        surf.getRevRateConstants(rdot.data());
        surf.getFwdRatesOfProgress(frop.data());
        surf.getRevRatesOfProgress(rrop.data());
        for (size_t k = 0; k < surf.nReactions(); k++) {
            std::cout<<"No: "<<k<<" "<<surf.reactionString(k)<<" Type:"<<surf.reactionType(k)<<std::endl;
            std::cout<<"dG(J/kmol)= "<<dG[k]<<" dG0(J/kmol)= "<<dG0[k]<<std::endl;
            std::cout<<"dH(J/kmol)= "<<dH[k]<<" dS(J/(kmol.K))= "<<dS[k]<<std::endl;
            std::cout<<"G_cal(J/kmol)=dH-TdS= " <<dH[k]-T*dS[k]<<std::endl;
            std::cout<<"Keq = " <<kc[k]<<" K_cal=exp(-dG0/RT)= "<<std::exp(-dG0[k]/(GasConstant*T));
            std::cout<<" K_cal=c0P/c0R*exp(-dG0/RT)= "<<tp2->standardConcentration(0)/tp1->standardConcentration(0)*(std::exp(-dG0[k]/(GasConstant*T)))<<std::endl;
            std::cout<<"kf = "<< fdot[k] << " frop = " << frop[k] << " kr = " << rdot[k] << " rrop = " << rrop[k] << " Keq = kf/kr = " << fdot[k]/rdot[k]<< std::endl;
            std::cout<<std::endl;
        }
        vector_fp wdot(surf.nTotalSpecies());
        surf.getNetProductionRates(wdot.data());
        std::cout<<"Printing species ROP..."<<std::endl;
        for (size_t kk = 0; kk < surf.nTotalSpecies(); kk++)
            std::cout << surf.kineticsSpeciesName(kk)<< " wdot = " << wdot[kk] <<std::endl;
    } catch (CanteraError& err){
        std::cout<<err.what()<< std::endl;
    }
}
