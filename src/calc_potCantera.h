/*
 * calc_potCantera.h
 *
 *  Created on: 31.01.2019
 *      Author: Manik
 */

#ifndef CALC_POTCANTERA_H_
#define CALC_POTCANTERA_H_

#include <algorithms/algorithms.h>
#include "cantera/thermo.h"
#include "cantera/Interface.h"
#include "cantera/Edge.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/Interface.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/kinetics/importKinetics.h"

enum phaseType {CA,EL,AN};

void initCanteraSPM();
double calc_potCantera(double phiL, std::string surface, Cantera::compositionMap speciesMoleFrac);
double calc_ilocCantera(double phiS);
size_t getPhaseIdxbyName(std::string phaseName, Cantera::Interface *surface);

#endif /* CALC_POTCANTERA_H_ */
