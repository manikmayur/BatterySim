/*
 * calc_itotCantera.h
 *
 *  Created on: 15 May 2018
 *      Author: Manik
 */

#ifndef CALC_ITOTCANTERA_H_
#define CALC_ITOTCANTERA_H_

#include "cantera/thermo.h"
#include "cantera/Interface.h"
#include "cantera/Edge.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/Interface.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/kinetics/importKinetics.h"

void calc_itotCantera2s(doublereal cA, doublereal cB, doublereal t, doublereal &itot, doublereal phis);
void calc_itotCantera3s(doublereal cA, doublereal cB, doublereal cC, doublereal t, doublereal &itot, doublereal phis);
void thermoInitCathodeROP_demo(std::string xmlFile);


#endif /* CALC_ITOTCANTERA_H_ */
