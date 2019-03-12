/*
 * calc_potCantera.h
 *
 *  Created on: 31.01.2019
 *      Author: Manik
 */

#ifndef CANTERAFUNCTIONS_H_
#define CANTERAFUNCTIONS_H_

#include "algorithms/algorithms.h"
#include "cantera/thermo.h"
#include "cantera/Interface.h"
#include "cantera/Edge.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/Interface.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/kinetics/importKinetics.h"

enum phaseType {phCA,phEL,phAN};

void initCanteraCV();
void calc_ropCanteraCV_2s(doublereal cA, doublereal cB, doublereal t, doublereal *rop, doublereal phis);
//
void initCanteraSPM();
double calc_potCantera(std::string surfName, double xLi, double phiL, double Iapp, double T);
double calc_entropyCantera(std::string surfName, double xLi, double T);
double calc_ilocCantera(double phiS);
size_t getPhaseIdxbyName(std::string phaseName, Cantera::Interface *surface);

#endif /* CANTERAFUNCTIONS_H_ */
