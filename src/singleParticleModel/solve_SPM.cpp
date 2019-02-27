/* -----------------------------------------------------------------
 * Programmer(s): Manik Mayur
 *
 * This example solves a Single Particle Model for intercalation in
 * Lithium-ion batteries
 *
 * Example problem:
 *
 * An ODE system is generated from the following 2-species diurnal
 * kinetics advection-diffusion PDE system in 1 space dimensions:
 *
 * dX(i)/dt = (1/r^2)*d/dr(r^2*D*(d/dr)X(i))      for i = 1,2,   where
 * D is a concentration dependent
 * The problem is posed on the line
 *   0 <= x <= L,
 * with homogeneous Neumann boundary conditions, and for time t in
 *   0 <= t <= 2*tp sec.
 * The PDE system is treated by central differences on a uniform
 * 100 x 100 mesh, with simple polynomial initial profiles.
 * The problem is solved with CVODE, with the BDF/GMRES,
 * BDF/Bi-CGStab, and BDF/TFQMR methods (i.e. using the SUNSPGMR,
 * SUNSPBCGS and SUNSPTFQMR linear solvers) and the block-diagonal
 * part of the Newton matrix as a left preconditioner. A copy of
 * the block-diagonal part of the Jacobian is saved and
 * conditionally reused within the Precond routine.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>

#include <cvode/cvode.h>                 /* main integrator header file       */
#include <cvode/cvode_spils.h>           /* access to CVSpils interface       */
#include <sunlinsol/sunlinsol_spgmr.h>   /* access to SPGMR SUNLinearSolver   */
#include <sunlinsol/sunlinsol_spbcgs.h>  /* access to SPBCGS SUNLinearSolver  */
#include <sunlinsol/sunlinsol_sptfqmr.h> /* access to SPTFQMR SUNLinearSolver */
#include <nvector/nvector_serial.h>      /* serial N_Vector types, fct. and macros */
#include <sundials/sundials_dense.h>     /* use generic DENSE solver in preconditioning */
#include <sundials/sundials_types.h>     /* definition of realtype */
#include <sundials/sundials_math.h>      /* contains the macros ABS, SUNSQR, and EXP */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sunmatrix/sunmatrix_band.h>  /* access to band SUNMatrix             */
#include <sunlinsol/sunlinsol_band.h>  /* access to band SUNLinearSolver       */
#include <cvode/cvode_direct.h>        /* access to CVDls interface            */
#include "cantera/canteraFunctions.h"
#include "singleParticleModel/parameters_SPM.h"
#include "solvep2D.h"

/* Problem Constants */

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)
#define TWO  RCONST(2.0)

#define DT		RCONST(p_tTotal/p_NT)       /* number of seconds in two hours  */

#define XMIN	ZERO                 /* grid boundaries in x  */
#define XMAX    RCONST(ONE)

#define MR		100             /* MR = number of x mesh points */

/* CVodeInit Constants */

//#define NEQ     (nVariables*MR)  /* NEQ = number of equations */

/* Linear Solver Loop Constants */

#define USE_SPGMR   0
#define USE_SPBCG   1
#define USE_SPTFQMR 2

/* User-defined vector and matrix accessor macros: IJth, IJth */

/* IJth is defined in order to isolate the translation from the
   mathematical 3-dimensional structure of the dependent variable vector
   to the underlying 1-dimensional storage. IJth is defined in order to
   write code which indexes into dense matrices with a (row,column)
   pair, where 1 <= row, column <= nVariables.

   IJth(vdata,i,j,k) references the element in the vdata array for
   species i at mesh point (j,k), where 1 <= i <= nVariables,
   0 <= j <= MR-1, 0 <= k <= MY-1. The vdata array is obtained via
   the call vdata = N_VGetArrayPointer(v), where v is an N_Vector.
   For each mesh point (j,k), the elements for species i and i+1 are
   contiguous within vdata.
 */
#define IJth(vdata,i,j) (vdata[i-1 + (j)*nVariables])

/* Type : UserData
   contains preconditioner blocks, pivot arrays, and problem constants */

typedef struct
{
	double dr,hdcoA,hdcoB;
} *UserData;

size_t nVariables;                 /* number of variables */
/* Private Helper Functions */

static UserData AllocUserData(void);
static void InitUserData(UserData data);
static void FreeUserData(UserData data);
static void SetInitProfiles(N_Vector u);
int runSolver();
static void PrintOutput(void *cvode_mem, N_Vector u, realtype t, FILE *fp);
static void PrintFinalStats(void *cvode_mem, int linsolver);
static int check_flag(void *flagvalue, const char *funcname, int opt);

/* Functions Called by the Solver */

int (*f)(double, N_Vector, N_Vector, void *); /* solver right hand function template */
static int fSPM(double t, N_Vector u, N_Vector udot, void *user_data);
static int fSPMT(double t, N_Vector u, N_Vector udot, void *user_data);

/* Main Program */

int main(void)
{
	//int val = runSolver();
	int val = runP2D();
}

int runSolver()
{
	double t;
	N_Vector u;
	UserData data;
	//SUNMatrix A;
	SUNLinearSolver LS;
	SUNMatrix A;
	void *cvode_mem;
	int flag;
	size_t NEQ;  /* NEQ = number of equations */
	std::string outFile;
	size_t iout = 1;
	FILE * pFile;

	//A = NULL;
	u = NULL;
	data = NULL;
	LS = NULL;
	cvode_mem = NULL;
	try
	{
		switch(p_model)
		{
			case SPM:
				f = &fSPM;
				nVariables = 2;
				pFile = fopen("outSPM.dat","w");
				fprintf(pFile, "# t xLi_ca xLi_an V_ca V_an\n");
				break;
			case SPMT:
				f = &fSPMT;
				nVariables = 3;
				pFile = fopen("outSPMT.dat","w");
				fprintf(pFile, "# t xLi_ca xLi_an T V_ca V_an\n");
				break;
			default:
				std::cout<<"Unknown model"<<std::endl;
		}
		NEQ = nVariables*MR;

		/* Allocate memory, and set problem data, initial values, tolerances */
		u = N_VNew_Serial(NEQ);
		if(check_flag((void *)u, "N_VNew_Serial", 0)) return(1);
		data = AllocUserData();
		if(check_flag((void *)data, "AllocUserData", 2)) return(1);
		InitUserData(data);
		SetInitProfiles(u);
		/* Call CVodeCreate to create the solver memory and specify the
		 * Backward Differentiation Formula and the use of a Newton iteration */
		cvode_mem = CVodeCreate(CV_BDF);
		if(check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);

		/* Set the pointer to user-defined data */
		flag = CVodeSetUserData(cvode_mem, data);
		if(check_flag(&flag, "CVodeSetUserData", 1)) return(1);

		/* Call CVodeInit to initialize the integrator memory and specify the
		 * user's right hand side function in u'=f(t,u), the inital time p_tInit, and
		 * the initial dependent variable vector u. */
		flag = CVodeInit(cvode_mem, f, p_tInit, u);
		if(check_flag(&flag, "CVodeInit", 1)) return(1);

		/* Call CVodeSStolerances to specify the scalar relative tolerance
		 * and scalar absolute tolerances */
		flag = CVodeSStolerances(cvode_mem, RTOL, ATOL);
		if (check_flag(&flag, "CVodeSStolerances", 1)) return(1);

		/* Call SUNSPGMR to specify the linear solver SPGMR with
         left preconditioning and the default maximum Krylov dimension */
		//LS = SUNSPGMR(u, PREC_NONE, 0);
		//if(check_flag((void *)LS, "SUNSPGMR", 0)) return(1);

		/* Create banded SUNMatrix for use in linear solves */
		A = SUNBandMatrix(NEQ, MR, MR);
		if(check_flag((void *)A, "SUNBandMatrix", 0)) return(1);

		/* Create banded SUNLinearSolver object */
		LS = SUNLinSol_Band(u, A);
		if(check_flag((void *)LS, "SUNLinSol_Band", 0)) return(1);

		flag = CVodeSetLinearSolver(cvode_mem, LS, A);
		//flag = CVSpilsSetLinearSolver(cvode_mem, LS);
		if(check_flag(&flag, "CVSpilsSetLinearSolver", 1)) return 1;

		/* Initialite Cantera */
		initCanteraSPM();
		//double test = IJth(u,3,1);
		//std::cout<<nVariables<<" "<<test<<std::endl;
		// 2nd
		/* Create dense SUNMatrix for use in linear solves */
		//A = SUNDenseMatrix(NEQ, NEQ);
		//if(check_flag((void *)A, "SUNDenseMatrix", 0)) return(1);

		/* Create dense SUNLinearSolver object for use by CVode */
		//LS = SUNDenseLinearSolver(u, A);
		//if(check_flag((void *)LS, "SUNDenseLinearSolver", 0)) return(1);

		/* Call CVDlsSetLinearSolver to attach the matrix and linear solver to CVode */
		//flag = CVDlsSetLinearSolver(cvode_mem, LS, A);
		//if(check_flag(&flag, "CVDlsSetLinearSolver", 1)) return(1);

		//2nd

		/* In loop over output points, call CVode, print results, test for error */
		printf(" \nLIB Single Particle Model\n\n");

		for (double tout = DT; iout <= p_NT; iout++, tout += DT)
		{
			flag = CVode(cvode_mem, tout, u, &t, CV_NORMAL);
			PrintOutput(cvode_mem, u, t, pFile);
			if(check_flag(&flag, "CVode", 1)) break;
		}

		PrintFinalStats(cvode_mem, 0);

		/* Free memory */
		fclose(pFile);
		//Cantera::appdelete();
		N_VDestroy(u);
		FreeUserData(data);
		CVodeFree(&cvode_mem);
		SUNLinSolFree(LS);
		return(0);
	}
	catch(std::exception &err)
	{
		printf(err.what());
	}
}

/* Private helper functions */

/* Allocate memory for data structure of type UserData */

static UserData AllocUserData(void)
{
	UserData data;
	data = (UserData) malloc(sizeof *data);
	return(data);
}

/* Load problem constants in data */

static void InitUserData(UserData data)
{
	data->dr = (XMAX-XMIN)/(MR-1); // jx = 0, x = XMIN, jx = 1, x = XMAX, x = XMIN + jx*dx
	data->hdcoA = p_DLi_ca(Tref)/SUNSQR(data->dr);
	data->hdcoB = p_DLi_an(Tref)/SUNSQR(data->dr);
}

/* Free data memory */

static void FreeUserData(UserData data)
{
	free(data);
}

/* Set initial conditions in u */

static void SetInitProfiles(N_Vector u)
{
	realtype *udata;

	/* Set pointer to data array in vector u. */
	udata = N_VGetArrayPointer(u);

	/* Load initial profiles of cA and cB into u vector */
	switch (p_model)
	{
	case SPM:
		for (int jx=0; jx < MR; jx++)
		{
			IJth(udata,1,jx) = p_xLimin_ca;
			IJth(udata,2,jx) = p_xLimax_an;
		}
		break;
	case SPMT:
		for (int jx=0; jx < MR; jx++)
		{
			IJth(udata,1,jx) = p_xLimin_ca;
			IJth(udata,2,jx) = p_xLimax_an;
			IJth(udata,3,jx) = p_Tamb;
		}
		break;
	default:
		std::cout<<"Unknown model"<<std::endl;
	}
}

/* Print current t, step count, order, stepsize, and sampled c1,c2 values */
static void PrintOutput(void *cvode_mem, N_Vector u, realtype t, FILE *fp)
{
	long int nst;
	int qu, flag;
	realtype hu, *udata, potAN, potCA;
	int mx1 = MR - 1;

	udata = N_VGetArrayPointer(u);

	flag = CVodeGetNumSteps(cvode_mem, &nst);
	check_flag(&flag, "CVodeGetNumSteps", 1);
	flag = CVodeGetLastOrder(cvode_mem, &qu);
	check_flag(&flag, "CVodeGetLastOrder", 1);
	flag = CVodeGetLastStep(cvode_mem, &hu);
	check_flag(&flag, "CVodeGetLastStep", 1);
	potCA = calc_potCantera(p_nameCathodeSurf, IJth(udata,1,mx1), p_Iapp*p_Rel(Tref), p_Iapp);
	potAN = calc_potCantera(p_nameAnodeSurf, IJth(udata,2,mx1), 0.0, p_Iapp);
	switch (p_model)
	{
	case SPM:
		fprintf(fp, "%.2e %12.3e %12.3e %12.3e %12.3e\n",t, IJth(udata,1,mx1), IJth(udata,2,mx1), potCA+3.0, potAN+3.0);
		break;
	case SPMT:
		fprintf(fp, "%.2e %12.3e %12.3e %12.3e %12.3e %12.3e\n",t, IJth(udata,1,mx1), IJth(udata,2,mx1), IJth(udata,3,mx1), potCA+3.0, potAN+3.0);
		break;
	default:
		std::cout<<"Unknown model"<<std::endl;
		break;
	}
}

/* Get and print final statistics */
static void PrintFinalStats(void *cvode_mem, int linsolver)
{
	long int lenrw, leniw ;
	long int lenrwLS, leniwLS;
	long int nst, nfe, nsetups, nni, ncfn, netf;
	long int nli, npe, nps, ncfl, nfeLS;
	int flag;

	flag = CVodeGetWorkSpace(cvode_mem, &lenrw, &leniw);
	check_flag(&flag, "CVodeGetWorkSpace", 1);
	flag = CVodeGetNumSteps(cvode_mem, &nst);
	check_flag(&flag, "CVodeGetNumSteps", 1);
	flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
	check_flag(&flag, "CVodeGetNumRhsEvals", 1);
	flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
	check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
	flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
	check_flag(&flag, "CVodeGetNumErrTestFails", 1);
	flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
	check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
	flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
	check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);

	flag = CVSpilsGetWorkSpace(cvode_mem, &lenrwLS, &leniwLS);
	check_flag(&flag, "CVSpilsGetWorkSpace", 1);
	flag = CVSpilsGetNumLinIters(cvode_mem, &nli);
	check_flag(&flag, "CVSpilsGetNumLinIters", 1);
	flag = CVSpilsGetNumPrecEvals(cvode_mem, &npe);
	check_flag(&flag, "CVSpilsGetNumPrecEvals", 1);
	flag = CVSpilsGetNumPrecSolves(cvode_mem, &nps);
	check_flag(&flag, "CVSpilsGetNumPrecSolves", 1);
	flag = CVSpilsGetNumConvFails(cvode_mem, &ncfl);
	check_flag(&flag, "CVSpilsGetNumConvFails", 1);
	flag = CVSpilsGetNumRhsEvals(cvode_mem, &nfeLS);
	check_flag(&flag, "CVSpilsGetNumRhsEvals", 1);

	printf("\nFinal Statistics.. \n\n");
	printf("lenrw   = %5ld     leniw   = %5ld\n"  , lenrw, leniw);
	printf("lenrwLS = %5ld     leniwLS = %5ld\n"  , lenrwLS, leniwLS);
	printf("nst     = %5ld\n"                     , nst);
	printf("nfe     = %5ld     nfeLS   = %5ld\n"  , nfe, nfeLS);
	printf("nni     = %5ld     nli     = %5ld\n"  , nni, nli);
	printf("nsetups = %5ld     netf    = %5ld\n"  , nsetups, netf);
	printf("npe     = %5ld     nps     = %5ld\n"  , npe, nps);
	printf("ncfn    = %5ld     ncfl    = %5ld\n\n", ncfn, ncfl);

	if (linsolver < 2)
		printf("======================================================================\n\n");
}

/* Check function return value...
     opt == 0 means SUNDIALS function allocates memory so check if
              returned NULL pointer
     opt == 1 means SUNDIALS function returns a flag so check if
              flag >= 0
     opt == 2 means function allocates memory so check if returned
              NULL pointer */

static int check_flag(void *flagvalue, const char *funcname, int opt)
{
	int *errflag;

	/* Check if SUNDIALS function returned NULL pointer - no memory allocated */
	if (opt == 0 && flagvalue == NULL) {
		fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
				funcname);
		return(1); }

	/* Check if flag < 0 */
	else if (opt == 1) {
		errflag = (int *) flagvalue;
		if (*errflag < 0) {
			fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
					funcname, *errflag);
			return(1); }}

	/* Check if function returned NULL pointer - no memory allocated */
	else if (opt == 2 && flagvalue == NULL) {
		fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
				funcname);
		return(1); }

	return(0);
}

/* Functions called by the solver */

/* f routine. Compute RHS function f(t,u). */
static int fSPM(realtype t, N_Vector u, N_Vector udot, void *user_data)
{
	realtype cA, cB, cAlt, cBlt, cArt, cBrt;
	realtype hordA, hordB;
	realtype hordcoA, hordcoB;
	realtype *udata, *dudata;
	int jx;
	UserData data;

	data   = (UserData) user_data;
	udata  = N_VGetArrayPointer(u);
	dudata = N_VGetArrayPointer(udot);

	/* Loop over all grid points. */

	for (jx=0; jx < MR; jx++) {

		/* Extract c1 and c2, and set kinetic rate terms. */
		cA = IJth(udata,1,jx);
		cB = IJth(udata,2,jx);
		//if (cA<0 || cB<0) return(1);

		/* Set horizontal diffusion terms. */
		/*
		 * At r = 1; DdX/dr = Ri = -2dr*(rCa/(DCa*cMaxCa))*iApp/(n*F*sCa)
		 * c_j+1 - c_j-1 = -2dr*(rCa/(DCa*cMaxCa))*iApp/(n*F*sCa)
		 * At r = 0; DdX/dr = 0
		 * IJth(vdata,i,j,k) references the element in the vdata array for
		 * species i at mesh point (j,k), where 1 <= i <= nVariables,
		 * 0 <= j <= MR-1, 0 <= k <= MY-1.
		 */

		cAlt = (jx == 0) ? IJth(udata,1,1) : IJth(udata,1,jx-1);
		cBlt = (jx == 0) ? IJth(udata,2,1) : IJth(udata,2,jx-1);
		cArt = (jx == MR-1) ? IJth(udata,1,MR-2) - TWO*data->dr*(p_rP_ca/(p_DLi_ca(Tref)*p_csMax_ca))*(p_Iapp/(F*p_S_ca)) : IJth(udata,1,jx+1);
		cBrt = (jx == MR-1) ? IJth(udata,2,MR-2) + TWO*data->dr*(p_rP_an/(p_DLi_an(Tref)*p_csMax_an))*(p_Iapp/(F*p_S_an)) : IJth(udata,2,jx+1);
		//DA/SUNSQR(data->dr) * ( IJth(udata,1,jx+iright) - 2*IJth(udata,1,jx) + IJth(udata,1,jx+ileft))
		hordA = p_DLi_ca(Tref)/(TWO*jx*SUNSQR(data->dr))*((jx+2)*cArt - TWO*jx*cA + (jx-2)*cAlt)
				+ p_DLi_ca(Tref)/(TWO*SUNSQR(data->dr))*(cArt - cA)
				- p_DLi_ca(Tref)/(TWO*SUNSQR(data->dr))*(cA - cAlt);
		hordcoA = (jx == 0) ? 6*p_DLi_ca(Tref)/SUNSQR(data->dr)*(IJth(udata,1,1)-IJth(udata,1,0)) : hordA;
		//DB/SUNSQR(data->dr) * ( IJth(udata,1,jx+iright) - 2*IJth(udata,1,jx) + IJth(udata,1,jx+ileft))
		hordB = p_DLi_an(Tref)/(TWO*jx*SUNSQR(data->dr))*((jx+2)*cBrt - TWO*jx*cB + (jx-2)*cBlt)
				+ p_DLi_an(Tref)/(TWO*SUNSQR(data->dr))*(cBrt - cB)
				- p_DLi_an(Tref)/(TWO*SUNSQR(data->dr))*(cB - cBlt);
		hordcoB = (jx == 0) ? 6*p_DLi_an(Tref)/SUNSQR(data->dr)*(IJth(udata,1,2)-IJth(udata,1,0)) : hordB;
		/* Load all terms into udot. */
		/*Equation: dci/dt = (Di*dci/dx)
		 * dci/dt = Di*(c1rt - TWO*c1 + c1lt)/dx2
		 * */

		IJth(dudata, 1, jx) = hordcoA/SUNSQR(p_rP_ca);
		IJth(dudata, 2, jx) = hordcoB/SUNSQR(p_rP_an);
	}

	return(0);
}

/* f routine. Compute RHS function f(t,u). */
static int fSPMT(realtype t, N_Vector u, N_Vector udot, void *user_data)
{
	realtype cA, cB, cAlt, cBlt, cArt, cBrt;
	realtype hordA, hordB;
	realtype hordcoA, hordcoB;
	realtype *udata, *dudata;
	realtype qIrr, qRev, qOut, T;
	realtype Uca, Uan, Ucell, Ueq, dUdTca, dUdTan;
	int jx;
	UserData data;

	data   = (UserData) user_data;
	udata  = N_VGetArrayPointer(u);
	dudata = N_VGetArrayPointer(udot);

	cA = IJth(udata,1,MR-1);
	cB = IJth(udata,2,MR-1);
	T = IJth(udata,3,MR-1);
	Uca = calc_potCantera(p_nameCathodeSurf, cA, p_Iapp*p_Rel(T), p_Iapp, T);
	Uan = calc_potCantera(p_nameAnodeSurf, cB, 0.0, p_Iapp, T);
	Ucell = Uca - Uan;
	Uca = calc_potCantera(p_nameCathodeSurf, cA, 0.0, 0.0, T);
	Uan = calc_potCantera(p_nameAnodeSurf, cB, 0.0, 0.0, T);
	Ueq = Uca - Uan;
	dUdTca = calc_entropyCantera(p_nameCathodeSurf, cA, T)/F;
	dUdTan = calc_entropyCantera(p_nameAnodeSurf, cB, T)/F;

	/* Loop over all grid points. */

	for (jx=0; jx < MR; jx++) {

		/* Extract c1 and c2, and set kinetic rate terms. */
		cA = IJth(udata,1,jx);
		cB = IJth(udata,2,jx);
		T = IJth(udata,3,jx);

		//if (cA<0 || cB<0) return(1);

		/* Set horizontal diffusion terms. */
		/*
		 * At r = 1; DdX/dr = Ri = -2dr*(rCa/(DCa*cMaxCa))*iApp/(n*F*sCa)
		 * c_j+1 - c_j-1 = -2dr*(rCa/(DCa*cMaxCa))*iApp/(n*F*sCa)
		 * At r = 0; DdX/dr = 0
		 * IJth(vdata,i,j,k) references the element in the vdata array for
		 * species i at mesh point (j,k), where 1 <= i <= nVariables,
		 * 0 <= j <= MR-1, 0 <= k <= MY-1.
		 */

		cAlt = (jx == 0) ? IJth(udata,1,1) : IJth(udata,1,jx-1);
		cBlt = (jx == 0) ? IJth(udata,2,1) : IJth(udata,2,jx-1);
		cArt = (jx == MR-1) ? IJth(udata,1,MR-2) - TWO*data->dr*(p_rP_ca/(p_DLi_ca(Tref)*p_csMax_ca))*(p_Iapp/(F*p_S_ca)) : IJth(udata,1,jx+1);
		cBrt = (jx == MR-1) ? IJth(udata,2,MR-2) + TWO*data->dr*(p_rP_an/(p_DLi_an(Tref)*p_csMax_an))*(p_Iapp/(F*p_S_an)) : IJth(udata,2,jx+1);
		//DA/SUNSQR(data->dr) * ( IJth(udata,1,jx+iright) - 2*IJth(udata,1,jx) + IJth(udata,1,jx+ileft))
		hordA = p_DLi_ca(Tref)/(TWO*jx*SUNSQR(data->dr))*((jx+2)*cArt - TWO*jx*cA + (jx-2)*cAlt)
				+ p_DLi_ca(Tref)/(TWO*SUNSQR(data->dr))*(cArt - cA)
				- p_DLi_ca(Tref)/(TWO*SUNSQR(data->dr))*(cA - cAlt);
		hordcoA = (jx == 0) ? 6*p_DLi_ca(Tref)/SUNSQR(data->dr)*(IJth(udata,1,1)-IJth(udata,1,0)) : hordA;
		//DB/SUNSQR(data->dr) * ( IJth(udata,1,jx+iright) - 2*IJth(udata,1,jx) + IJth(udata,1,jx+ileft))
		hordB = p_DLi_an(Tref)/(TWO*jx*SUNSQR(data->dr))*((jx+2)*cBrt - TWO*jx*cB + (jx-2)*cBlt)
				+ p_DLi_an(Tref)/(TWO*SUNSQR(data->dr))*(cBrt - cB)
				- p_DLi_an(Tref)/(TWO*SUNSQR(data->dr))*(cB - cBlt);
		hordcoB = (jx == 0) ? 6*p_DLi_an(Tref)/SUNSQR(data->dr)*(IJth(udata,1,2)-IJth(udata,1,0)) : hordB;
		/* Load all terms into udot. */
		/*Equation: dci/dt = (Di*dci/dx)
		 * dci/dt = Di*(c1rt - TWO*c1 + c1lt)/dx2
		 */
		qIrr = p_Iapp*(Ucell-Ueq);
		qRev = p_Iapp*T*(dUdTca-dUdTan);
		qOut = p_hA*(T-p_Tamb);

		IJth(dudata, 1, jx) = hordcoA/SUNSQR(p_rP_ca);
		IJth(dudata, 2, jx) = hordcoB/SUNSQR(p_rP_an);
		IJth(dudata, 3, jx) = 1/(p_rhoCell*p_volCell*p_cpCell)*(qIrr+qRev-qOut);
	}
	return(0);
}
