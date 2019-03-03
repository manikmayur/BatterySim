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
static void openCircuitPotential(size_t jx, double cS,double T, double &U, double &dUdT);
/* Load problem constants in data */

static void InitUserData(UserData data)
{
	data->dr = (XMAX-XMIN)/(NX-1); // jx = 0, x = XMIN, jx = 1, x = XMAX, x = XMIN + jx*dx
	data->dx = (XMAX-XMIN)/(NX-1); // jx = 0, x = XMIN, jx = 1, x = XMAX, x = XMIN + jx*dx
	size_t i=0;
	data->idxCsAvg = ++i;
	data->idxCe = ++i;
	data->idxT = ++i;
	data->idxphiS = ++i;
	data->idxCs = ++i;
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
		IJth(udata,data->idxCsAvg,jx) = (dom.domType==CA||dom.domType==AN)?dom.cLiInit:ZERO;
		IJth(udata,data->idxCs,jx) = (dom.domType==CA||dom.domType==AN)?dom.cLiInit:ZERO;
		IJth(udata,data->idxCe,jx) = (jx<ca.idx0||jx>an.idxL)?ZERO:p_cE;
		IJth(udata,data->idxT,jx) = Tref;
		IJth(udata,data->idxphiS,jx) = ZERO;
		IJth(udata,data->idxphiL,jx) = ZERO;
		//
		IJth(iddata,data->idxCs,jx) = ZERO;
		IJth(iddata,data->idxphiS,jx) = ZERO;
		IJth(iddata,data->idxphiL,jx) = ZERO;
		IJth(updata,data->idxCs,jx) = ZERO;
		IJth(updata,data->idxphiS,jx) = ZERO;
		IJth(updata,data->idxphiL,jx) = ZERO;
	}
	/* Initialize up vector to 0. */
	 N_VConst(ZERO, up);

	  /* heatres sets res to negative of ODE RHS values at interior points. */
	 heatres(ZERO, uu, up, res, data);

	  /* Copy -res into up to get correct interior initial up values. */
	 N_VScale(-ONE, res, up);

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
	std::vector<FILE*> files(NSPECIES+1);

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

	/* Call IDACalcIC (with default options) to correct the initial values. */

	retval = IDACalcIC(mem, IDA_YA_YDP_INIT, t1);
	if(check_retval(&retval, "IDACalcIC", 1)) return(1);

	/* Print output heading. */
	files[data->idxCsAvg-1] = fopen("outP2D_concSLiAvg.dat","w");
	files[data->idxCs-1] = fopen("outP2D_concSLi.dat","w");
	files[data->idxCe-1] = fopen("outP2D_concElyte.dat","w");
	files[data->idxT-1] = fopen("outP2D_Temperature.dat","w");
	files[data->idxphiS-1] = fopen("outP2D_phiS.dat","w");
	files[data->idxphiL-1] = fopen("outP2D_phiL.dat","w");
	files[data->idxphiL-1] = fopen("outP2D_phiL.dat","w");
	files[data->idxphiL] = fopen("outP2D_Vcell.dat","w");

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
	fclose(files[data->idxphiL]);
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
	realtype diff_T, diff_Ce, sT, sCe, sCs, sCsAvg, sphiS, sphiL;
	realtype qIrr, qRev, qOut, qTot, qOhm;
	double dx2lt, dx2rt, betal, betar, difflt, diffrt;
	realtype Ucell, Ueq, dUdT, jloc=0.0;
	domain dom;
	UserData data;

	Ucell=Ueq=dUdT=ZERO;
	qIrr=qRev=qOut=qTot=qOhm=ZERO;
	sT=sCe=sCs=sCsAvg=sphiS=sphiL=ZERO;

	udata = N_VGetArrayPointer(uu);
	updata = N_VGetArrayPointer(up);
	resv = N_VGetArrayPointer(resval);

	data = (UserData)user_data;

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
	g = TWO*(1-p_tp)*R/F;

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
		openCircuitPotential(jx, Cs, T, Ueq, dUdT);
		jloc = TWO*rateConst(jx,T)*std::sqrt(Ce*(dom.cLiMax-Cs)*Cs)*std::sinh(0.5*R/(F*T)*(phiS-phiL-Ueq));
		jloc = (dom.domType==AL||dom.domType==EL||dom.domType==CU)?
				ZERO:jloc;
		//
		dx2lt = (jx == dom.idx0 && jx != al.idx0)?
				dx(jx)*(dx(jx-1)+dx(jx))/TWO:SUNSQR(dx(jx));
		dx2rt = (jx == dom.idxL && jx != cu.idxL)?
				dx(jx)*(dx(jx)+dx(jx+1))/TWO:SUNSQR(dx(jx));
		betal = (jx == dom.idx0 && jx != al.idx0)?
				dx(jx-1)/(dx(jx)+dx(jx-1)):0.5;
		betar = (jx == dom.idxL && jx != cu.idxL)?
				dx(jx)/(dx(jx)+dx(jx+1)):0.5;
		// Ce
		Celt = (jx == ca.idx0) ? Ce : IJth(udata,data->idxCe,jx-1);
		Cert = (jx == an.idxL) ? Ce : IJth(udata,data->idxCe,jx+1);
		// T
		Tlt = ONE/(ONE+(p_h*al.dx/(TWO*al.kappaS)))*(p_h*al.dx/al.kappaS*p_Tamb
			+(ONE-p_h*al.dx/(TWO*al.kappaS))*IJth(udata,data->idxT,al.idx0));
		Tlt = (jx == al.idx0) ? Tlt : IJth(udata,data->idxT,jx-1);
		Trt = ONE/(ONE+(p_h*cu.dx/(TWO*cu.kappaS)))*(p_h*cu.dx/cu.kappaS*p_Tamb
			+(ONE-p_h*cu.dx/(TWO*cu.kappaS))*IJth(udata,data->idxT,cu.idxL));
		Trt = (jx == cu.idxL) ? Trt : IJth(udata,data->idxT,jx+1);

		// CsAvg
		sCsAvg = (dom.domType==CA||dom.domType==AN)?
				-3.0*jloc/dom.rP:ZERO;
		IJth(resv,data->idxCsAvg, jx) = IJth(updata,data->idxCsAvg,jx) - sCsAvg;

		// Cs
		sCs = (dom.domType==CA||dom.domType==AN)?
				-dom.rP*jloc/(diffS(jx,T)*5):ZERO;
		IJth(resv,data->idxCs,jx) = IJth(udata,data->idxCs,jx) - IJth(udata,data->idxCsAvg,jx) - sCs;

		// T
		difflt = (jx == dom.idx0 && jx != al.idx0)?
				kappaS(jx)*kappaS(jx-1)/(betal*kappaS(jx)+(1-betal)*kappaS(jx-1)):kappaS(jx);
		diffrt = (jx == dom.idxL && jx != cu.idxL)?
				kappaS(jx)*kappaS(jx+1)/(betar*kappaS(jx+1)+(1-betar)*kappaS(jx)):kappaS(jx);
		diff_T = diffrt/dx2rt*(Trt-T)-difflt/dx2lt*(T-Tlt);
		qIrr = F*dom.aLi*jloc*(phiS-phiL-Ueq);
		qRev = F*dom.aLi*jloc*T*dUdT;
		qOut = p_h*(T-p_Tamb);
		qTot = qIrr+qRev;
		sT = (dom.domType==AL||dom.domType==CU)?
				SUNSQR(p_Iapp)/sigmaS(jx):qTot;
		IJth(resv,data->idxT,jx) = IJth(updata,data->idxT,jx) - 1/(dom.rho*dom.cP)*(diff_T + sT);

		// Ce
		difflt = (jx == dom.idx0 && jx > ca.idx0)?
				diffL(dom.por,Ce,T)*diffL((getDomain(jx-1)).por,Celt,Tlt)/(betal*diffL(dom.por,Ce,T)+(1-betal)*diffL((getDomain(jx-1)).por,Celt,Tlt)):diffL(dom.por,Ce,T);
		diffrt = (jx == dom.idxL && jx < an.idxL)?
				diffL(dom.por,Ce,T)*diffL((getDomain(jx+1)).por,Cert,Trt)/(betar*diffL((getDomain(jx+1)).por,Cert,Trt)+(1-betar)*diffL(dom.por,Ce,T)):diffL(dom.por,Ce,T);
		diff_Ce = diffrt/dx2rt*(Cert-Ce)-difflt/dx2lt*(Ce-Celt);
		sCe = dom.aLi*(1-p_tp)*jloc;
		double rpor = (jx>=ca.idx0 && jx<=an.idxL)?1/dom.por:ZERO;
		IJth(resv,data->idxCe,jx) = IJth(updata,data->idxCe,jx) - rpor*(diff_Ce + sCe);

		// phiS
		if (jx == al.idx0)
			phiSlt = phiS-dom.dx*p_Iapp/dom.sigmaS;
		else if (jx == an.idx0)
			phiSlt = phiS;
		else if (dom.domType==AL||dom.domType==CA||dom.domType==AN||dom.domType==CU)
			phiSlt = IJth(udata,data->idxphiS,jx-1);
		else
			phiSlt = ZERO;
		if (jx == ca.idxL)
			phiSrt = phiS;
		else if (jx == cu.idxL)
			phiSrt = phiS+dom.dx*p_Iapp/dom.sigmaS;
		else if (dom.domType==AL||dom.domType==CA||dom.domType==AN||dom.domType==CU)
			phiSrt = IJth(udata,data->idxphiS,jx+1);
		else
			phiSrt = ZERO;
		difflt = (jx == ca.idx0 || jx == cu.idx0)?
				sigmaS(jx)*sigmaS(jx-1)/(betal*sigmaS(jx)+(1-betal)*sigmaS(jx-1)):sigmaS(jx);
		diffrt = (jx == al.idxL|| jx == an.idxL)?
				sigmaS(jx)*sigmaS(jx+1)/(betar*sigmaS(jx+1)+(1-betar)*sigmaS(jx)):sigmaS(jx);
		diff_phiS = diffrt/dx2rt*(phiSrt-phiS)-difflt/dx2lt*(phiS-phiSlt);
		sphiS = dom.aLi*F*jloc;
		//IJth(resv,data->idxphiS,jx) = IJth(updata,data->idxphiS,jx) - rCdl*(diff_phiS + sphiS);
		double revD = (dom.domType==AL||dom.domType==CA||dom.domType==AN||dom.domType==CU)?
					1/(diffrt/dx2rt+difflt/dx2lt):ZERO;
		IJth(resv,data->idxphiS,jx) = IJth(udata,data->idxphiS,jx)
								-revD*(diffrt/dx2rt*phiSrt + difflt/dx2lt*phiSlt*0.999999999984 - sphiS);
		// phiL
		phiLlt = (jx == ca.idx0) ? phiL : IJth(udata,data->idxphiL,jx-1);
		phiLrt = (jx == an.idxL) ? -phiL : IJth(udata,data->idxphiL,jx+1);
		difflt = (jx == dom.idx0 && jx > ca.idx0)?
				sigmaL(dom.por,Ce,T)*sigmaL((getDomain(jx-1)).por,Celt,Tlt)/(betal*sigmaL(dom.por,Ce,T)+(1-betal)*sigmaL((getDomain(jx-1)).por,Celt,Tlt)):sigmaL(dom.por,Ce,T);
		diffrt = (jx == dom.idxL && jx < an.idxL)?
				sigmaL(dom.por,Ce,T)*sigmaL((getDomain(jx+1)).por,Cert,Trt)/(betar*sigmaL((getDomain(jx+1)).por,Cert,Trt)+(1-betar)*sigmaL(dom.por,Ce,T)):sigmaL(dom.por,Ce,T);
		diff_phiL1 = (jx >= ca.idx0 && jx <= an.idxL)?
				diffrt/dx2rt*(phiLrt-phiL)-difflt/dx2lt*(phiL-phiLlt):ZERO;
		diff_phiL2 = (jx >= ca.idx0 && jx <= an.idxL)?
				diffrt/dx2rt*g*(T+Trt)/TWO*log(Cert/Ce)-difflt/dx2lt*g*(T+Tlt)/TWO*log(Ce/Celt):ZERO;
		sphiL = -dom.aLi*F*jloc;
		//IJth(resv,data->idxphiL,jx) = IJth(updata,data->idxphiL,jx) - rCdl*(diff_phiL1 + diff_phiL2 + sphiL);
		revD = (dom.domType==CA||dom.domType==EL||dom.domType==AN)?
							1/(diffrt/dx2rt+difflt/dx2lt):ZERO;
		IJth(resv,data->idxphiL,jx) = IJth(udata,data->idxphiL,jx)
				- revD*(diffrt/dx2rt*phiLrt + difflt/dx2lt*phiLlt - diff_phiL2 - sphiL);
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

/*
 * Print Output
 */

static void printOutputP2D(void *mem, realtype t, N_Vector uu, std::vector<FILE*>fp,
		void *user_data)
{
	realtype *udata;
	double x = ZERO;
	UserData data;

	data = (UserData)user_data;
	udata = N_VGetArrayPointer(uu);

	for (size_t jx=0; jx<NX; jx++)
	{
		domain dom = getDomain(jx);
		x = (jx == 0) ? (x+dom.dx/2):(x+dom.dx);
		fprintf(fp[data->idxCsAvg-1],"%.5e %I64u %12.7e %12.7e\n",t, jx, x, IJth(udata,data->idxCsAvg,jx));
		fprintf(fp[data->idxCs-1],"%.5e %I64u %12.7e %12.7e\n",t, jx, x, IJth(udata,data->idxCs,jx));
		fprintf(fp[data->idxCe-1],"%.5e %I64u %12.7e %12.7e\n",t, jx, x, IJth(udata,data->idxCe,jx));
		fprintf(fp[data->idxT-1],"%.5e %I64u %12.7e %12.7e\n",t, jx, x, IJth(udata,data->idxT,jx));
		fprintf(fp[data->idxphiS-1],"%.5e %I64u %12.7e %12.7e\n",t, jx, x, IJth(udata,data->idxphiS,jx));
		fprintf(fp[data->idxphiL-1],"%.5e %I64u %12.7e %12.7e\n",t, jx, x, IJth(udata,data->idxphiL,jx));
	}
	fprintf(fp[data->idxphiL],"%.2e %12.7e %12.3e\n",t,IJth(udata,data->idxphiS,al.idx0),IJth(udata,data->idxphiS,cu.idxL));
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

static void openCircuitPotential(size_t jx, double cS,double T, double &U, double &dUdT)
{
	// Calculate the open circuit voltage of the battery in the positive
	// electrode
	domain dom = getDomain(jx);
	double theta  = cS/dom.cLiMax;

	if (dom.domType == CA)
	{
		// Compute the variation of OCV with respect to temperature variations [V/K]
		dUdT = -0.001*(0.199521039-0.928373822*theta + 1.364550689000003*std::pow(theta,2)-0.6115448939999998*std::pow(theta,3));
		dUdT = dUdT/(1-5.661479886999997*theta +11.47636191*std::pow(theta,2)-9.82431213599998*std::pow(theta,3)+3.048755063*std::pow(theta,4));

		// Define the OCV for the positive electrode
		U = -4.656+88.669*std::pow(theta,2) - 401.119*std::pow(theta,4) + 342.909*std::pow(theta,6) - 462.471*std::pow(theta,8) + 433.434*std::pow(theta,10);
		U = U/(-1+18.933*std::pow(theta,2)-79.532*std::pow(theta,4)+37.311*std::pow(theta,6)-73.083*std::pow(theta,8)+95.96*std::pow(theta,10));
		U = U + (T-Tref)*dUdT;
	}
	else if (dom.domType == AN)
	{
		// Compute the variation of OCV with respect to temperature variations [V/K]
		dUdT = 0.001*(0.005269056 +3.299265709*theta-91.79325798*std::pow(theta,2)+1004.911008*std::pow(theta,3)-5812.278127*std::pow(theta,4)
		+ 19329.7549*std::pow(theta,5) - 37147.8947*std::pow(theta,6) + 38379.18127*std::pow(theta,7)-16515.05308*std::pow(theta,8));
		dUdT = dUdT/(1-48.09287227*theta+1017.234804*std::pow(theta,2)-10481.80419*std::pow(theta,3)+59431.3*std::pow(theta,4)-195881.6488*std::pow(theta,5)
		+ 374577.3152*std::pow(theta,6) - 385821.1607*std::pow(theta,7) + 165705.8597*std::pow(theta,8));

		// Define the OCV for the negative electrode
		U = 0.7222 + 0.1387*theta + 0.029*std::pow(theta,0.5) - 0.0172/theta + 0.0019/std::pow(theta,1.5) + 0.2808*std::exp(0.9-15*theta)-0.7984*std::exp(0.4465*theta - 0.4108);
		U = U + (T-Tref)*dUdT;
	}
	else
	{
		U = ZERO;
		dUdT = ZERO;
	}
}
