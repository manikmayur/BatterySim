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
#include "parameters_2s.h"

void initCantera();
void calc_ropCantera2S(doublereal cA, doublereal cB, doublereal t, doublereal rop[nSpecies+1], doublereal phis);

#endif /* CALC_ITOTCANTERA_H_ */
