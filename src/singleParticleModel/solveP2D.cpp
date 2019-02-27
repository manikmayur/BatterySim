/*
 * solveP2D.cpp
 *
 *  Created on: 18.02.2019
 *      Author: Manik
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <ida/ida.h>                   /* prototypes for IDA fcts., consts.    */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h>        /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h>        /* access to dense SUNLinearSolver      */
#include <sundials/sundials_types.h>   /* definition of type realtype */
#include <sundials/sundials_math.h>      /* contains the macros ABS, SUNSQR, and EXP */
#include "singleParticleModel/parameters_SPM.h"
#include "cantera/canteraFunctions.h"

/* Problem Constants */

#define NOUT  100
#define NX	NXTOTAL             /* NX = number of x mesh points */
#define NSPECIES 6
#define NEQ   NSPECIES*NX
#define ZERO  RCONST(0.0)
#define ONE   RCONST(1.0)
#define TWO   RCONST(2.0)
#define BVAL  RCONST(0.1)
#define DT	RCONST(p_tTotal/p_NT)       /* number of seconds in two hours  */

#define XMIN	ZERO                 /* grid boundaries in x  */
#define XMAX    RCONST(ONE)

/* Type: UserData */

typedef struct {
	sunindextype mm;
	realtype dr;
	realtype dx;
	size_t idxCsAvg, idxCs, idxCe, idxT, idxphiS, idxphiL;
} *UserData;

/* Prototypes of functions called by IDA */

int heatres(realtype tres, N_Vector uu, N_Vector up, N_Vector resval, void *user_data);
#define IJth(vdata,i,j) (vdata[i-1 + (j)*NSPECIES])

/* Prototypes of private functions */

static void printOutputP2D(void *mem, realtype t, N_Vector u, std::vector<FILE*> fp, void *user_data);
static int SetInitialProfile(UserData data, N_Vector uu, N_Vector up,
		N_Vector id, N_Vector res);
static int check_retval(void *returnvalue, const char *funcname, int opt);
static void openCircuitPotential(size_t jx, double cS,double T, double &U_p, double &dudt_p, double &U_n, double &dudt_n);
/* Load problem constants in data */

static void InitUserData(UserData data)
{
	data->dr = (XMAX-XMIN)/(NX-1); // jx = 0, x = XMIN, jx = 1, x = XMAX, x = XMIN + jx*dx
	data->dx = (XMAX-XMIN)/(NX-1); // jx = 0, x = XMIN, jx = 1, x = XMAX, x = XMIN + jx*dx
	size_t i=0;
	data->idxCsAvg = ++i;
	data->idxCs = ++i;
	data->idxCe = ++i;
	data->idxT = ++i;
	data->idxphiS = ++i;
	data->idxphiL = ++i;
}


static int SetInitialProfile(UserData data, N_Vector uu, N_Vector up,
		N_Vector id, N_Vector res)
{
	realtype *udata, *updata, *iddata;
	domain dom;
	udata = N_VGetArrayPointer(uu);
	updata = N_VGetArrayPointer(up);
	iddata = N_VGetArrayPointer(id);

	  /* Initialize id to 1's. */
	  N_VConst(ONE, id);
	/* Initialize uu on all grid points. */
	for (size_t jx=0; jx < NX; jx++)
	{
		dom = getDomain(jx);
		IJth(udata,data->idxCsAvg,jx) = (dom.domType==CA||dom.domType==AN)?dom.xLiInit:ZERO;
		IJth(udata,data->idxCs,jx) = (dom.domType==CA||dom.domType==AN)?dom.xLiInit:ZERO;
		IJth(udata,data->idxCe,jx) = (jx<ca.idx0||jx>an.idxL)?ZERO:p_cE;
		IJth(udata,data->idxT,jx) = Tref;
		IJth(udata,data->idxphiS,jx) = ZERO;;
		IJth(udata,data->idxphiL,jx) = ZERO;;
	}
	/* Initialize up vector to 0. */
	 N_VConst(ZERO, up);

	  /* heatres sets res to negative of ODE RHS values at interior points. */
	  //heatres(ZERO, uu, up, res, data);

	  /* Copy -res into up to get correct interior initial up values. */
	 // N_VScale(-ONE, res, up);

	return(0);
}

// Setup and run solver
int runP2D(void)
{
	void *mem;
	UserData data;
	N_Vector uu, up, id, res;
	int retval;
	long int netf, ncfn;
	realtype rtol, atol, t0, t1, tret;
	SUNMatrix A;
	SUNLinearSolver LS;
	size_t iout = 1;
	std::vector<FILE*> files(NSPECIES);

	mem = NULL;
	data = NULL;
	uu = up = id = res = NULL;
	A = NULL;
	LS = NULL;

	/* Create vectors uu, up, res, constraints, id. */
	uu = N_VNew_Serial(NEQ);
	if(check_retval((void *)uu, "N_VNew_Serial", 0)) return(1);
	up = N_VNew_Serial(NEQ);
	if(check_retval((void *)up, "N_VNew_Serial", 0)) return(1);
	id = N_VNew_Serial(NEQ);
	if(check_retval((void *)id, "N_VNew_Serial", 0)) return(1);
	res = N_VNew_Serial(NEQ);
	if(check_retval((void *)res, "N_VNew_Serial", 0)) return(1);

	/* Create and load problem data block. */
	data = (UserData) malloc(sizeof *data);
	if(check_retval((void *)data, "malloc", 2)) return(1);

	InitUserData(data);

	/* Initialize uu, up, id. */
	SetInitialProfile(data, uu, up, id, res);

	/* Set remaining input parameters. */
	t0   = ZERO;
	t1   = RCONST(0.01);
	rtol = RTOL;
	atol = ATOL;

	/* Initialite Cantera */
	initCanteraSPM();
	/* Call IDACreate and IDAMalloc to initialize solution */
	mem = IDACreate();
	if(check_retval((void *)mem, "IDACreate", 0)) return(1);

	retval = IDASetUserData(mem, data);
	if(check_retval(&retval, "IDASetUserData", 1)) return(1);

	// Set which components are algebraic or differential
	retval = IDASetId(mem, id);
	if(check_retval(&retval, "IDASetId", 1)) return(1);

	retval = IDAInit(mem, heatres, t0, uu, up);
	if(check_retval(&retval, "IDAInit", 1)) return(1);

	retval = IDASStolerances(mem, rtol, atol);
	if(check_retval(&retval, "IDASStolerances", 1)) return(1);

	// Create dense SUNMatrix for use in linear solves
	A = SUNDenseMatrix(NEQ, NEQ);
	if(check_retval((void *)A, "SUNDenseMatrix", 0)) return(1);

	/* Create dense SUNLinearSolver object for use by CVode */
	LS = SUNLinSol_Dense(uu, A);
	if(check_retval((void *)LS, "SUNDenseLinearSolver", 0)) return(1);

	/* Call CVDlsSetLinearSolver to attach the matrix and linear solver to CVode */
	retval = IDASetLinearSolver(mem, LS, A);
	if(check_retval(&retval, "IDASetLinearSolver", 1)) return(1);


	/* Print output heading. */
	files[data->idxCsAvg-1] = fopen("outP2D_concSLiAvg.dat","w");
	files[data->idxCs-1] = fopen("outP2D_concSLi.dat","w");
	files[data->idxCe-1] = fopen("outP2D_concElyte.dat","w");
	files[data->idxT-1] = fopen("outP2D_Temperature.dat","w");
	files[data->idxphiS-1] = fopen("outP2D_phiS.dat","w");
	files[data->idxphiL-1] = fopen("outP2D_phiL.dat","w");

	printOutputP2D(mem, t0, uu, files, data);


	/* Loop over output times, call IDASolve, and print results. */

	for (double tout = t1; iout <= p_NT; iout++, tout += DT)
	{
		retval = IDASolve(mem, tout, &tret, uu, up, IDA_NORMAL);
		if(check_retval(&retval, "IDASolve", 1)) return(1);
		printOutputP2D(mem, tret, uu, files, data);
	}

	/* Print remaining counters and free memory. */
	retval = IDAGetNumErrTestFails(mem, &netf);
	check_retval(&retval, "IDAGetNumErrTestFails", 1);
	retval = IDAGetNumNonlinSolvConvFails(mem, &ncfn);
	check_retval(&retval, "IDAGetNumNonlinSolvConvFails", 1);
	printf("\n netf = %ld,   ncfn = %ld \n", netf, ncfn);
	fclose(files[data->idxCsAvg-1]);
	fclose(files[data->idxCs-1]);
	fclose(files[data->idxCe-1]);
	fclose(files[data->idxT-1]);
	fclose(files[data->idxphiS-1]);
	fclose(files[data->idxphiL-1]);
	IDAFree(&mem);
	SUNLinSolFree(LS);
	SUNMatDestroy(A);
	N_VDestroy(uu);
	N_VDestroy(up);
	N_VDestroy(res);
	free(data);

	return(0);
}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY IDA
 *--------------------------------------------------------------------
 */

/*
 * heatres: heat equation system residual function
 * This uses 5-point central differencing on the interior points, and
 * includes algebraic equations for the boundary values.
 * So for each interior point, the residual component has the form
 *    res_i = u'_i - (central difference)_i
 * while for each boundary point, it is res_i = u_i.
 */

int heatres(realtype tres, N_Vector uu, N_Vector up, N_Vector resval,
		void *user_data)
{
	realtype *udata, *updata, *resv;
	realtype Cs, CsAvg, Ce, Celt, Cert, T, Tlt, Trt;
	realtype phiS, phiSlt, phiSrt, diff_phiS;
	realtype phiL, phiLlt, phiLrt, diff_phiL1, diff_phiL2, g;
	realtype diff_T, diff_Ce, sCe, sCs, sCsAvg, sphiS, sphiL;
	realtype qIrr, qRev, qOut, qTot, qOhm;
	realtype Uca, Uan, Ucell, Ueq, dUdTca, dUdTan, iloc=0.0;
	domain dom;
	UserData data;

	Uca=Uan=Ucell=Ueq=dUdTca=dUdTan=ZERO;
	qIrr=qRev=qOut=qTot=qOhm=ZERO;
	sCe=sCs=sCsAvg=sphiS=sphiL=ZERO;

	udata = N_VGetArrayPointer(uu);
	updata = N_VGetArrayPointer(up);
	resv = N_VGetArrayPointer(resval);

	data = (UserData)user_data;

	CsAvg = IJth(udata,data->idxCsAvg,NX-1);
	Cs = IJth(udata,data->idxCs,NX-1);
	Ce = IJth(udata,data->idxCe,NX-1);
	T = IJth(udata,data->idxT,NX-1);
	//phiS = IJth(udata,data->idxphiS,NX-1);
	//phiL = IJth(udata,data->idxphiS,NX-1);
	//Cs = 0.5;
	/*
	Uca = calc_potCantera(p_nameCathodeSurf, 0.5, p_Iapp*p_Rel(Tref), p_Iapp, Tref);
	Uan = calc_potCantera(p_nameAnodeSurf, 0.2, 0.0, p_Iapp, Tref);
	Ucell = Uca - Uan;
	Uca = calc_potCantera(p_nameCathodeSurf, 0.5, 0.0, 0.0, Tref);
	Uan = calc_potCantera(p_nameAnodeSurf, 0.2, 0.0, 0.0, Tref);
	Ueq = Uca - Uan;
	dUdTca = calc_entropyCantera(p_nameCathodeSurf, 0.5, Tref)/F;
	dUdTan = calc_entropyCantera(p_nameAnodeSurf, 0.2, Tref)/F;*/

	// Loop over interior points; set res = up - (central difference).
	// Loop over all grid points.

	for (size_t jx=0; jx < NX; jx++)
	{
		// Extract c1 and c2, and set kinetic rate terms
		dom = getDomain(jx);
		//
		CsAvg = IJth(udata,data->idxCsAvg,jx);
		Cs = IJth(udata,data->idxCs,jx);
		Ce = IJth(udata,data->idxCe,jx);
		T = IJth(udata,data->idxT,jx);
		phiS = IJth(udata,data->idxphiS,jx);
		phiL = IJth(udata,data->idxphiL,jx);
		//
		if (dom.domType==CA||dom.domType==AN)
			openCircuitPotential(jx, Cs,T, Uca, dUdTca, Uan, dUdTan);
		//
		Celt = (jx == ca.idx0) ? Ce : IJth(udata,data->idxCe,jx-1);
		Cert = (jx == an.idxL) ? Ce : IJth(udata,data->idxCe,jx+1);
		//
		Tlt = ONE/(ONE+(p_h*al.dx/(TWO*al.kappaS)))*(p_h*al.dx/al.kappaS*p_Tamb
			+(ONE-p_h*al.dx/(TWO*al.kappaS))*IJth(udata,data->idxT,al.idx0));
		Tlt = (jx == al.idx0) ? Tlt : IJth(udata,data->idxT,jx-1);
		Trt = ONE/(ONE+(p_h*cu.dx/(TWO*cu.kappaS)))*(p_h*cu.dx/cu.kappaS*p_Tamb
			+(ONE-p_h*cu.dx/(TWO*cu.kappaS))*IJth(udata,data->idxT,cu.idxL));
		Trt = (jx == cu.idxL) ? Trt : IJth(udata,data->idxT,jx+1);
		//
		phiSlt = (jx == ca.idx0) ? phiS+dom.dx*p_Iapp/dom.sigmaS : IJth(udata,data->idxphiS,jx-1);
		phiSlt = (jx == an.idx0) ? phiS : IJth(udata,data->idxphiS,jx-1);
		phiSrt = (jx == ca.idxL) ? phiS : IJth(udata,data->idxphiS,jx+1);
		phiSrt = (jx == an.idxL) ? phiS+dom.dx*p_Iapp/dom.sigmaS : IJth(udata,data->idxphiS,jx+1);
		//
		phiLlt = (jx == ca.idx0) ? phiL : IJth(udata,data->idxphiL,jx-1);
		phiLrt = (jx == an.idxL) ? -phiL : IJth(udata,data->idxphiL,jx+1);
		//
		diff_Ce = diffL(jx)/dx2(jx)*(Cert-Ce)-diffL(jx-1)/dx2(jx-1)*(Ce-Celt);
		diff_T = kappaS(jx)/dx2(jx)*(Trt-T)-kappaS(jx-1)/dx2(jx-1)*(T-Tlt);
		diff_phiS = (dom.domType==AL||dom.domType==EL||dom.domType==CU)?
				ZERO:sigmaS(jx)/dx2(jx)*(phiSrt-phiS)-sigmaS(jx-1)/dx2(jx-1)*(phiS-phiSlt);
		diff_phiL1 = (dom.domType==AL||dom.domType==CU)?
				ZERO:sigmaL(jx)/dx2(jx)*(phiLrt-phiL)-sigmaL(jx-1)/dx2(jx-1)*(phiL-phiLlt);
		g = TWO*(1-p_tp)*R/F;
		diff_phiL2 = (dom.domType==AL||dom.domType==CU)?
				ZERO:sigmaL(jx)*(T+Trt)*g/(TWO*dx2(jx))*log(Cert/Ce)-sigmaL(jx-1)*(T+Tlt)*g/(TWO*dx2(jx-1))*log(Ce/Celt);
		// Set source terms
		Ueq = Uca - Uan;
		qIrr = p_Iapp*(Ucell-Ueq);
		qRev = p_Iapp*T*(dUdTca-dUdTan);
		qOut = p_h*(T-p_Tamb);

		//iloc = TWO*dom.sigmaL*dom.cLiMax*std::sqrt(Ce*(1-Cs)*Cs)*std::sinh(0.5*R/(F*T)*(phiS-phiL-Ueq));
		sCsAvg = (dom.domType==AL||dom.domType==EL||dom.domType==CU)?
				ZERO:-3.0*iloc/dom.rP;
		sCs = (dom.domType==AL||dom.domType==EL||dom.domType==CU)?
				ZERO:-dom.rP*iloc/(dom.diffS*5);
		sCe = (dom.domType==AL||dom.domType==EL||dom.domType==CU)?
				ZERO:dom.aLi*(1-p_tp)*iloc;
		sphiS = (dom.domType==AL||dom.domType==EL||dom.domType==CU)?
				ZERO:dom.aLi*F*iloc;
		sphiL = (dom.domType==AL||dom.domType==EL||dom.domType==CU)?
				ZERO:-dom.aLi*F*iloc;

		//

		IJth(resv,data->idxCsAvg, jx) = IJth(updata,data->idxCsAvg,jx) - sCsAvg;
		IJth(resv,data->idxCs, jx) = IJth(updata,data->idxCs,jx);
		IJth(resv,data->idxCe,jx) = dom.por*IJth(updata,data->idxCe,jx) - diff_Ce - sCe;
		IJth(resv,data->idxT,jx) = dom.rho*dom.cP*IJth(updata,data->idxT,jx) - diff_T;
		IJth(resv,data->idxphiS,jx) = IJth(updata,data->idxphiS,jx);
		IJth(resv,data->idxphiL,jx) = IJth(updata,data->idxphiL,jx);
	}

	return(0);

}

/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/*
 * SetInitialProfile: routine to initialize u, up, and id vectors.
 */

/* Free data memory */

static void FreeUserData(UserData data)
{
	free(data);
}

/*
 * Print Output
 */

static void printOutputP2D(void *mem, realtype t, N_Vector uu, std::vector<FILE*>fp,
		void *user_data)
{
	int retval;
	realtype umax, hused;
	realtype *udata;
	double x = ZERO;
	long int nst, nni, nje, nre, nreLS;
	int kused;
	UserData data;

	data = (UserData)user_data;
	umax = N_VMaxNorm(uu);
	udata = N_VGetArrayPointer(uu);

	retval = IDAGetLastOrder(mem, &kused);
	check_retval(&retval, "IDAGetLastOrder", 1);
	retval = IDAGetNumSteps(mem, &nst);
	check_retval(&retval, "IDAGetNumSteps", 1);
	retval = IDAGetNumNonlinSolvIters(mem, &nni);
	check_retval(&retval, "IDAGetNumNonlinSolvIters", 1);
	retval = IDAGetNumResEvals(mem, &nre);
	check_retval(&retval, "IDAGetNumResEvals", 1);
	retval = IDAGetLastStep(mem, &hused);
	check_retval(&retval, "IDAGetLastStep", 1);
	retval = IDAGetNumJacEvals(mem, &nje);
	check_retval(&retval, "IDAGetNumJacEvals", 1);
	retval = IDAGetNumLinResEvals(mem, &nreLS);
	check_retval(&retval, "IDAGetNumLinResEvals", 1);

	for (size_t jx=0; jx<NX; jx++)
	{
		domain dom = getDomain(jx);
		x = (jx == 0) ? (x+dom.dx/2):(x+dom.dx);
		fprintf(fp[data->idxCsAvg-1],"%.2e %I64u %12.7e %12.3e\n",t, jx, x, IJth(udata,data->idxCsAvg,jx));
		fprintf(fp[data->idxCs-1],"%.2e %I64u %12.7e %12.3e\n",t, jx, x, IJth(udata,data->idxCs,jx));
		fprintf(fp[data->idxCe-1],"%.2e %I64u %12.7e %12.3e\n",t, jx, x, IJth(udata,data->idxCe,jx));
		fprintf(fp[data->idxT-1],"%.2e %I64u %12.7e %12.3e\n",t, jx, x, IJth(udata,data->idxT,jx));
		fprintf(fp[data->idxphiS-1],"%.2e %I64u %12.7e %12.3e\n",t, jx, x, IJth(udata,data->idxphiS,jx));
		fprintf(fp[data->idxphiL-1],"%.2e %I64u %12.7e %12.3e\n",t, jx, x, IJth(udata,data->idxphiL,jx));
	}
}

/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns an integer value so check if
 *            retval < 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer
 */

static int check_retval(void *returnvalue, const char *funcname, int opt)
{
	int *retval;

	/* Check if SUNDIALS function returned NULL pointer - no memory allocated */
	if (opt == 0 && returnvalue == NULL) {
		fprintf(stderr,
				"\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
				funcname);
		return(1);
	} else if (opt == 1) {
		/* Check if retval < 0 */
		retval = (int *) returnvalue;
		if (*retval < 0) {
			fprintf(stderr,
					"\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
					funcname, *retval);
			return(1);
		}
	} else if (opt == 2 && returnvalue == NULL) {
		/* Check if function returned NULL pointer - no memory allocated */
		fprintf(stderr,
				"\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
				funcname);
		return(1);
	}

	return(0);
}

static void openCircuitPotential(size_t jx, double cS,double T, double &U_p, double &dudt_p, double &U_n, double &dudt_n)
{
	// Calculate the open circuit voltage of the battery in the positive
	// electrode
	domain dom = getDomain(jx);
	double theta_p  = cS/dom.cLiMax;

	// Compute the variation of OCV with respect to temperature variations [V/K]
	dudt_p = -0.001*(0.199521039-0.928373822*theta_p + 1.364550689000003*std::pow(theta_p,2)-0.6115448939999998*std::pow(theta_p,3));
	dudt_p = dudt_p/(1-5.661479886999997*theta_p +11.47636191*std::pow(theta_p,2)-9.82431213599998*std::pow(theta_p,3)+3.048755063*std::pow(theta_p,4));

	// Define the OCV for the positive electrode
	U_p = -4.656+88.669*std::pow(theta_p,2) - 401.119*std::pow(theta_p,4) + 342.909*std::pow(theta_p,6) - 462.471*std::pow(theta_p,8) + 433.434*std::pow(theta_p,10);
	U_p = U_p/(-1+18.933*std::pow(theta_p,2)-79.532*std::pow(theta_p,4)+37.311*std::pow(theta_p,6)-73.083*std::pow(theta_p,8)+95.96*std::pow(theta_p,10));
	U_p = U_p + (T-Tref)*dudt_p;


	// Calculate the open circuit voltage of the battery in the negative
	// electrode
	double theta_n = cS/dom.cLiMax;

	// Compute the variation of OCV with respect to temperature variations [V/K]
	dudt_n = 0.001*(0.005269056 +3.299265709*theta_n-91.79325798*std::pow(theta_n,2)+1004.911008*std::pow(theta_n,3)-5812.278127*std::pow(theta_n,4)
		+ 19329.7549*std::pow(theta_n,5) - 37147.8947*std::pow(theta_n,6) + 38379.18127*std::pow(theta_n,7)-16515.05308*std::pow(theta_n,8));
	dudt_n = dudt_n/(1-48.09287227*theta_n+1017.234804*std::pow(theta_n,2)-10481.80419*std::pow(theta_n,3)+59431.3*std::pow(theta_n,4)-195881.6488*std::pow(theta_n,5)
	+ 374577.3152*std::pow(theta_n,6) - 385821.1607*std::pow(theta_n,7) + 165705.8597*std::pow(theta_n,8));

	// Define the OCV for the negative electrode
	U_n   = 0.7222 + 0.1387*theta_n + 0.029*std::pow(theta_n,0.5) - 0.0172/theta_n + 0.0019/std::pow(theta_n,1.5) + 0.2808*std::exp(0.9-15*theta_n)-0.7984*std::exp(0.4465*theta_n - 0.4108);
	U_n = U_n + (T-Tref)*dudt_n;
}
