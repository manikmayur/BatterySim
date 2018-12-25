/* -----------------------------------------------------------------
 * Programmer(s): Manik Mayur
 *
 * This example solves a diffusion only CV problem
 *
 * Example problem:
 *
 * An ODE system is generated from the following 2-species diurnal
 * kinetics advection-diffusion PDE system in 1 space dimensions:
 *
 * dc(i)/dt = Kh*(d/dx)^2 c(i) + Ri(c1,c2,t)      for i = 1,2,   where
 *   R1(c1,c2,t) = 0 ,
 *   R2(c1,c2,t) = 0 ,
 * Kh is a constant
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

#include "calc_itotCantera.h"

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
#include <cvode/cvode_direct.h>        /* access to CVDls interface            */
#include "parameters_3s.h"


/* Problem Constants */

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)
#define TWO  RCONST(2.0)

#define NUM_SPECIES  3            /* number of species         */
#define T0		ZERO              /* initial time */
#define NOUT	100               /* number of output times */
#define DT		RCONST(2*tp/NOUT) /* number of seconds in two hours  */
#define XMIN	ZERO              /* grid boundaries in x  */
#define XMAX    RCONST(L)
#define XMID    RCONST(L/2)       /* grid midpoints in x,y */
#define MX		100               /* MX = number of x mesh points */
#define NSMX	NUM_SPECIES*MX    /* NSMX = NUM_SPECIES*MX */

/* CVodeInit Constants */

#define RTOL    RCONST(1.0e-12)   /* scalar relative tolerance */
#define FLOOR   RCONST(1.0)       /* value of C1 or C2 at which tolerances */
/* change from relative to absolute      */
#define ATOL    (RTOL*FLOOR)      /* scalar absolute tolerance */
#define NEQ     (NUM_SPECIES*MX)  /* NEQ = number of equations */

/* Linear Solver Loop Constants */

#define USE_SPGMR   0
#define USE_SPBCG   1
#define USE_SPTFQMR 2

/* User-defined vector and matrix accessor macros: IJKth, IJth */

/* IJKth is defined in order to isolate the translation from the
   mathematical 3-dimensional structure of the dependent variable vector
   to the underlying 1-dimensional storage. IJth is defined in order to
   write code which indexes into dense matrices with a (row,column)
   pair, where 1 <= row, column <= NUM_SPECIES.

   IJKth(vdata,i,j,k) references the element in the vdata array for
   species i at mesh point (j,k), where 1 <= i <= NUM_SPECIES,
   0 <= j <= MX-1, 0 <= k <= MY-1. The vdata array is obtained via
   the call vdata = N_VGetArrayPointer(v), where v is an N_Vector.
   For each mesh point (j,k), the elements for species i and i+1 are
   contiguous within vdata.
 */

#define IJKth(vdata,i,j) (vdata[i-1 + (j)*NUM_SPECIES])

/* Type : UserData
   contains preconditioner blocks, pivot arrays, and problem constants */

typedef struct {
	realtype dx, hdco1, hdco2, hdco3;
} *UserData;

/* Private Helper Functions */

static UserData AllocUserData(void);
static void InitUserData(UserData data);
static void FreeUserData(UserData data);
static void SetInitialProfiles(N_Vector u, realtype dx);
static void PrintOutput(void *cvode_mem, N_Vector u, realtype t, FILE *fp);
static void PrintFinalStats(void *cvode_mem, int linsolver);
static int check_flag(void *flagvalue, const char *funcname, int opt);

/* Functions Called by the Solver */

static int f(realtype t, N_Vector u, N_Vector udot, void *user_data);
static int calc_itot(realtype cA, realtype cB, realtype cC, realtype t, realtype &itot, realtype &phis);

/* Main Program */

int main(void)
{
	realtype abstol, reltol, t, tout;
	N_Vector u;
	UserData data;
	SUNMatrix A;
	SUNLinearSolver LS;
	void *cvode_mem;
	int iout, flag;
	realtype itot=0, eta=0;

	A = NULL;
	u = NULL;
	data = NULL;
	LS = NULL;
	cvode_mem = NULL;
	try {
		/* Allocate memory, and set problem data, initial values, tolerances */
		u = N_VNew_Serial(NEQ);
		if(check_flag((void *)u, "N_VNew_Serial", 0)) return(1);
		data = AllocUserData();
		if(check_flag((void *)data, "AllocUserData", 2)) return(1);
		InitUserData(data);
		SetInitialProfiles(u, data->dx);
		abstol=ATOL;
		reltol=RTOL;

		/* Call CVodeCreate to create the solver memory and specify the
		 * Backward Differentiation Formula and the use of a Newton iteration */
		cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
		if(check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);

		/* Set the pointer to user-defined data */
		flag = CVodeSetUserData(cvode_mem, data);
		if(check_flag(&flag, "CVodeSetUserData", 1)) return(1);

		/* Call CVodeInit to initialize the integrator memory and specify the
		 * user's right hand side function in u'=f(t,u), the inital time T0, and
		 * the initial dependent variable vector u. */
		flag = CVodeInit(cvode_mem, f, T0, u);
		if(check_flag(&flag, "CVodeInit", 1)) return(1);

		/* Call CVodeSStolerances to specify the scalar relative tolerance
		 * and scalar absolute tolerances */
		flag = CVodeSStolerances(cvode_mem, reltol, abstol);
		if (check_flag(&flag, "CVodeSStolerances", 1)) return(1);

		/* Call SUNSPGMR to specify the linear solver SPGMR with
         left preconditioning and the default maximum Krylov dimension */
		LS = SUNSPGMR(u, PREC_NONE, 0);
		if(check_flag((void *)LS, "SUNSPGMR", 0)) return(1);

		flag = CVSpilsSetLinearSolver(cvode_mem, LS);
		if(check_flag(&flag, "CVSpilsSetLinearSolver", 1)) return 1;

		/* In loop over output points, call CVode, print results, test for error */
		printf(" \n3-species diffusion problem\n\n");
		FILE * pFile;
		pFile = fopen ("CV_3Species.dat","w");
		for (iout=1, tout = DT; iout <= NOUT; iout++, tout += DT) {
			flag = CVode(cvode_mem, tout, u, &t, CV_NORMAL);
			PrintOutput(cvode_mem, u, t, pFile);
			if(check_flag(&flag, "CVode", 1)) break;
		}

		PrintFinalStats(cvode_mem, 0);

		/* Free memory */
		fclose (pFile);

		N_VDestroy(u);
		FreeUserData(data);
		CVodeFree(&cvode_mem);
		SUNLinSolFree(LS);

		return(0);
	}
	catch(std::exception &err){
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
	data->dx = (XMAX-XMIN)/(MX-1); // jx = 0, x = XMIN, jx = 1, x = XMAX, x = XMIN + jx*dx
	data->hdco1 = DA/SUNSQR(data->dx);
	data->hdco2 = DB/SUNSQR(data->dx);
	data->hdco3 = DC/SUNSQR(data->dx);
}

/* Free data memory */

static void FreeUserData(UserData data)
{
	free(data);
}

/* Set initial conditions in u */

static void SetInitialProfiles(N_Vector u, realtype dx)
{
	int jx;
	realtype *udata;
	/* Set pointer to data array in vector u. */
	udata = N_VGetArrayPointer(u);

	/* Load initial profiles of c1 and c2 into u vector */

	for (jx=0; jx < MX; jx++) {
		IJKth(udata,1,jx) = c_bulk;
		IJKth(udata,2,jx) = c_bulk_p1;
		IJKth(udata,2,jx) = c_bulk_p2;
	}
}

/* Print current t, step count, order, stepsize, and sampled c1,c2 values */

static void PrintOutput(void *cvode_mem, N_Vector u, realtype t, FILE *fp)
{
	long int nst;
	int qu, flag;
	realtype hu, *udata, phis=0.0, itot=0.0;
	int mxh = MX/2 - 1, mx1 = MX - 1;

	udata = N_VGetArrayPointer(u);

	flag = CVodeGetNumSteps(cvode_mem, &nst);
	check_flag(&flag, "CVodeGetNumSteps", 1);
	flag = CVodeGetLastOrder(cvode_mem, &qu);
	check_flag(&flag, "CVodeGetLastOrder", 1);
	flag = CVodeGetLastStep(cvode_mem, &hu);
	check_flag(&flag, "CVodeGetLastStep", 1);

	printf("t = %.2e   no. steps = %ld   order = %d   stepsize = %.2e\n",
			t, nst, qu, hu);
	printf("c1 (bot.left/middle/top rt.) = %12.3e  %12.3e  %12.3e\n",
			IJKth(udata,1,0), IJKth(udata,1,mxh), IJKth(udata,1,mx1));
	printf("c2 (bot.left/middle/top rt.) = %12.3e  %12.3e  %12.3e\n\n",
			IJKth(udata,2,0), IJKth(udata,2,mxh), IJKth(udata,2,mx1));
	calc_itot(IJKth(udata,1,0), IJKth(udata,2,0), t, itot, phis);
	//printf ("%.2e %12.3e %12.3e %12.3e %12.3e\n",t, phis, IJKth(udata,1,0), IJKth(udata,2,0), itot);
	fprintf (fp, "%.2e %12.3e %12.3e %12.3e %12.3e\n",t, phis, IJKth(udata,1,0), IJKth(udata,2,0), IJKth(udata,3,0), itot);

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
static int f(realtype t, N_Vector u, N_Vector udot, void *user_data)
{
	realtype c1, c2, c3, c1lt, c2lt, c3lt, c1l, c2l, c3l;
	realtype c1rt, c2rt, c3rt, hord1, hord2, hord3;
	realtype hordco1, hordco2, hordco3;
	realtype *udata, *dudata, itot, phis;
	int jx, ileft, iright;
	UserData data;

	data   = (UserData) user_data;
	udata  = N_VGetArrayPointer(u);
	dudata = N_VGetArrayPointer(udot);

	/* Make local copies of problem variables, for efficiency. */
	hordco1  = data->hdco1;
	hordco2  = data->hdco2;
	hordco3  = data->hdco3;

	/* Loop over all grid points. */

	for (jx=0; jx < MX; jx++) {

		/* Extract c1 and c2, and set kinetic rate terms. */

		c1 = IJKth(udata,1,jx);
		c2 = IJKth(udata,2,jx);
		c3 = IJKth(udata,3,jx);

		c1l = IJKth(udata,1,0);
		c2l = IJKth(udata,2,0);
		c3l = IJKth(udata,3,0);

		//calc_itot(c1l, c2l, t, itot, phis);
		calc_itotCantera(c1l, c2l, c3l, t, itot, phis);

		/* Set horizontal diffusion terms. */
		/*
		 * At x = 0; Ddc/dx = Ri = -vi*itot/(n*F)
		 * c_j+1 - c_j = -(dx/Di)*vi*itot/(n*F)
		 * At x = L: u(L) = [cA0;cB0]
		 * IJKth(vdata,i,j,k) references the element in the vdata array for
		 * species i at mesh point (j,k), where 1 <= i <= NUM_SPECIES,
		 * 0 <= j <= MX-1, 0 <= k <= MY-1.
		 */

		ileft = (jx == 0) ? 1 : -1;
		iright = (jx == MX-1) ? -1 : 1;
		c1lt = (jx == 0) ? IJKth(udata,1,0) - vA*itot/(DA*F)*data->dx : IJKth(udata,1,jx-1);
		c2lt = (jx == 0) ? IJKth(udata,2,0) - vB*itot/(DB*F)*data->dx : IJKth(udata,2,jx-1);
		c3lt = (jx == 0) ? IJKth(udata,3,0) - vC*itot/(DB*F)*data->dx : IJKth(udata,3,jx-1);
		c1rt = (jx == MX-1) ? c_bulk : IJKth(udata,1,jx+1);
		c2rt = (jx == MX-1) ? c_bulk_p1 : IJKth(udata,2,jx+1);
		c3rt = (jx == MX-1) ? c_bulk_p2 : IJKth(udata,3,jx+1);
		//DA/SUNSQR(data->dx) * ( IJKth(udata,1,jx+iright) - 2*IJKth(udata,1,jx) + IJKth(udata,1,jx+ileft))
		hord1 = hordco1*(c1rt - TWO*c1 + c1lt);
		hord2 = hordco2*(c2rt - TWO*c2 + c2lt);
		hord3 = hordco3*(c3rt - TWO*c3 + c3lt);

		/* Load all terms into udot. */
		/*Equation: dci/dt = Di*d^2ci/dx^2
		 * dci/dt = Di*(c1rt - TWO*c1 + c1lt)/dx2
		 * */

		IJKth(dudata, 1, jx) = hord1;
		IJKth(dudata, 2, jx) = hord2;
		IJKth(dudata, 3, jx) = hord3;
	}
	return(0);
}
