/*
 * calcStress.cpp
 *
 *  Created on: 07.02.2019
 *      Author: Manik
 */
#include "algorithms/algorithms.h"
#include "singleParticleModel/parameters_SPM.h"

inline double cP(double ri) {return 0.0;}//interpolate(r,cLi_an(end,:),ri);}
inline double cS(double r) {return 0*r;}
inline double c0(double r) {return 0;}
inline double TP(double r) {return exp(-r)*0;}
inline double TS(double r) {return 0*r;}
inline double T0(double r) {return 0;}
inline double func_cP(double r) {return (cP(r) - c0(r))*pow(r,2);}
inline double func_cS(double r) {return 0*r;}
inline double int_cP(double r) {return integrate(func_cP,0.0,0.0,r);}
inline double int_cS(double r) {return integrate(func_cS,0.0,0.0,r);}
inline double func_TP(double r) {return (TP(r) - T0(r))*pow(r,2); }
inline double func_TS(double r) {return 0.*r;}
inline double int_TP(double r) {return integrate(func_TP,0.0,0.0,r)*0;}
inline double int_TS(double r) {return integrate(func_TS,0.0,0.0,r)*0;}

double aP()
{
	double aP0 = (1+(1+p_nuS)/(2*(1-2*p_nuS))*pow(p_RS/p_RP,3))/(1+(1+p_nuS)/(2*(1-2*p_nuS))*pow(p_RS/p_RP,3)
			-(p_modES/p_modEP)*(1-2*p_nuP)/(1-2*p_nuS)*(1-pow(p_RS/p_RP,3)));
	return aP0*(3*p_alphaP*pow(1/p_RP,3)*int_TP(p_RP)+p_omegaP*pow(1/p_RP,3)*int_cP(p_RP))
		    - (1+p_nuP)/(1-p_nuP)*p_alphaP*pow(1/p_RP,3)*int_TP(p_RP) - (1/3)*(1+p_nuP)/(1-p_nuP)*p_omegaP*pow(1/p_RP,3)*int_cP(p_RP);

}
double calc_sigtP(double r)
{
	double sigtP = p_alphaP*(p_modEP/(1-2*p_nuP))*(1/pow(r,3))*int_TP(r)
			-(1/3)*p_omegaP*(p_modEP/(1-2*p_nuP))*(1/pow(r,3))*int_cP(r)
	     	+ aP()*p_modEP/(1-2*p_nuP) + p_alphaP*(p_modEP/(1-p_nuP))*(TP(r)-T0(r))
			+ (1/3)*p_omegaP*(p_modEP/(1-p_nuP))*(cP(r)-c0(r));
	return sigtP;

}
