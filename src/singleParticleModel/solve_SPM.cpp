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

#include <ida/ida.h>                   /* prototypes for IDA fcts., consts.    */
#include <nvector/nvector_serial.h>      /* serial N_Vector types, fct. and macros */
#include <sundials/sundials_dense.h>     /* use generic DENSE solver in preconditioning */
#include <sundials/sundials_types.h>     /* definition of realtype */
#include <sundials/sundials_math.h>      /* contains the macros ABS, SUNSQR, and EXP */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
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

static void InitUserData(UserData data);
static void FreeUserData(UserData data);
static void SetInitialProfile(UserData data, N_Vector uu, N_Vector up, N_Vector id, N_Vector res);
int runSolver();
static void PrintOutput(void *mem, N_Vector u, realtype t, FILE *fp);
static int check_flag(void *flagvalue, const char *funcname, int opt);

/* Functions Called by the Solver */

int (*f)(double, N_Vector, N_Vector, N_Vector, void *); /* solver right hand function template */
static int fSPM(double t, N_Vector u, N_Vector udot, N_Vector resval,
        void *user_data);
//static int fSPMT(double t, N_Vector u, N_Vector udot, void *user_data);
static int fSPMT(realtype t, N_Vector u, N_Vector udot, N_Vector resval,
        void *user_data);

/* Main Program */

int main(void)
{
    int val = runSolver();
    //int val = runP2D();
}

int runSolver()
{
    void *mem;
    UserData data;
    N_Vector uu, up, id, res;
    int retval;
    realtype rtol, atol, t0, t1, tret;
    SUNMatrix A;
    SUNLinearSolver LS;
    size_t iout = 1;

    mem = NULL;
    data = NULL;
    uu = up = id = res = NULL;
    A = NULL;
    LS = NULL;

    size_t NEQ;  /* NEQ = number of equations */
    FILE *pFile;

    try
    {
        switch(p_model)
        {
            case SPM:
                f = &fSPM;
                nVariables = 4;
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
                throw std::runtime_error("Unknown model");
        }
        NEQ = nVariables*MR;

        /* Allocate memory, and set problem data, initial values, tolerances */
        uu = N_VNew_Serial(NEQ);
        if(check_flag((void *)uu, "N_VNew_Serial", 0)) return(1);

        up = N_VNew_Serial(NEQ);
        if(check_flag((void *)up, "N_VNew_Serial", 0)) return(1);

        id = N_VNew_Serial(NEQ);
        if(check_flag((void *)id, "N_VNew_Serial", 0)) return(1);

        res = N_VNew_Serial(NEQ);
        if(check_flag((void *)res, "N_VNew_Serial", 0)) return(1);

        /* Create and load problem data block. */
        data = (UserData) malloc(sizeof *data);
        if(check_flag((void *)data, "malloc", 2)) return(1);

        InitUserData(data);

        /* Initialite Cantera */
        initCanteraSPM();

        /* Initialize uu, up, id. */
        SetInitialProfile(data, uu, up, id, res);

        /* Set remaining input parameters. */
        t0   = ZERO;
        t1   = RCONST(0.01);
        rtol = RTOL;
        atol = ATOL;

        /* Call IDACreate and IDAMalloc to initialize solution */
        mem = IDACreate();
        if(check_flag((void *)mem, "IDACreate", 0)) return(1);

        retval = IDASetUserData(mem, data);
        if(check_flag(&retval, "IDASetUserData", 1)) return(1);

        // Set which components are algebraic or differential
        retval = IDASetId(mem, id);
        if(check_flag(&retval, "IDASetId", 1)) return(1);

        switch(p_model)
        {
            case SPM:
                retval = IDAInit(mem, fSPM, t0, uu, up);
                if(check_flag(&retval, "IDAInit", 1)) return(1);
                break;
            case SPMT:
                retval = IDAInit(mem, fSPMT, t0, uu, up);
                if(check_flag(&retval, "IDAInit", 1)) return(1);
                break;
            default:
                throw std::runtime_error("Unknown model");
        }

        retval = IDASStolerances(mem, rtol, atol);
        if(check_flag(&retval, "IDASStolerances", 1)) return(1);

        // Create dense SUNMatrix for use in linear solves
        A = SUNDenseMatrix(NEQ, NEQ);
        if(check_flag((void *)A, "SUNDenseMatrix", 0)) return(1);

        /* Create dense SUNLinearSolver object for use by CVode */
        LS = SUNLinSol_Dense(uu, A);
        if(check_flag((void *)LS, "SUNDenseLinearSolver", 0)) return(1);

        /* Call CVDlsSetLinearSolver to attach the matrix and linear solver to CVode */
        retval = IDASetLinearSolver(mem, LS, A);
        if(check_flag(&retval, "IDASetLinearSolver", 1)) return(1);

        /* Call IDACalcIC (with default options) to correct the initial values. */

        retval = IDACalcIC(mem, IDA_YA_YDP_INIT, t1);
        if(check_flag(&retval, "IDACalcIC", 1)) return(1);

        /* In loop over output points, call CVode, print results, test for error */
        printf(" \nLIB Single Particle Model\n\n");

        for (double tout = DT; iout <= p_NT; iout++, tout += DT)
        {
            retval = IDASolve(mem, tout, &tret, uu, up, IDA_NORMAL);
            PrintOutput(mem, uu, tret, pFile);
            if(check_flag(&retval, "IDASolve", 1)) break;
        }

        /* Free memory */
        fclose(pFile);
        FreeUserData(data);
        SUNLinSolFree(LS);
        IDAFree(&mem);
        SUNMatDestroy(A);
        N_VDestroy(uu);
        N_VDestroy(up);
        N_VDestroy(res);
        //free(data);
        return(0);
    }
    catch(std::exception &err)
    {
        printf(err.what());
    }
    catch(...)
    {
        std::cout<<"Unknown runtime exception\n";
    }
}

/* Private helper functions */

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

static void SetInitialProfile(UserData data, N_Vector uu, N_Vector up,
        N_Vector id, N_Vector res)
{
    realtype *udata, *iddata;

    udata = N_VGetArrayPointer(uu);
    iddata = N_VGetArrayPointer(id);

    /* Initialize id to 1's. */
    N_VConst(ONE, id);

    /* Initialize up vector to 0. */
    N_VConst(ZERO, up);

    /* Load initial profiles of cA and cB into u vector */
    switch (p_model)
    {
        case SPM:
            for (int jx=0; jx < MR; jx++)
            {
                IJth(udata,1,jx) = p_xLimin_ca;
                IJth(udata,2,jx) = p_xLimax_an;
                IJth(udata,3,jx) = ZERO;
                IJth(iddata,3,jx) = ZERO;
                IJth(udata,4,jx) = ZERO;
                IJth(iddata,4,jx) = ZERO;
            }
            /* fSPM sets res to negative of ODE RHS values at interior points. */
            fSPM(ZERO, uu, up, res, data);
            break;
        case SPMT:
            for (int jx=0; jx < MR; jx++)
            {
                IJth(udata,1,jx) = p_xLimin_ca;
                IJth(udata,2,jx) = p_xLimax_an;
                IJth(udata,3,jx) = p_Tamb;
            }
            /* fSPMT sets res to negative of ODE RHS values at interior points. */
            //fSPMT(ZERO, uu, up, res, data);
            break;
        default:
            std::cout<<"Unknown model"<<std::endl;
    }

    /* Copy -res into up to get correct interior initial up values. */
    N_VScale(-ONE, res, up);
}

/* Print current t, step count, order, stepsize, and sampled c1,c2 values */
static void PrintOutput(void *mem, N_Vector u, realtype t, FILE *fp)
{
    realtype *udata, potAN, potCA;
    int mx1 = MR - 1;

    udata = N_VGetArrayPointer(u);
    switch (p_model)
    {
        case SPM:
            fprintf(fp, "%.2e %12.3e %12.3e %12.3e %12.3e\n",t, IJth(udata,1,mx1), IJth(udata,2,mx1), IJth(udata,3,mx1)+3.0, IJth(udata,4,mx1)+3.0);
            break;
        case SPMT:
            potCA = calc_potCantera(p_nameCathodeSurf, IJth(udata,1,mx1), p_Iapp*p_Rel(IJth(udata,3,mx1)), p_Iapp, IJth(udata,3,mx1));
            potAN = calc_potCantera(p_nameAnodeSurf, IJth(udata,2,mx1), 0.0, p_Iapp, IJth(udata,3,mx1));
            fprintf(fp, "%.2e %12.3e %12.3e %12.3e %12.3e %12.3e\n",t, IJth(udata,1,mx1), IJth(udata,2,mx1), IJth(udata,3,mx1), potCA+3.0, potAN+3.0);
            break;
        default:
            std::cout<<"Unknown model"<<std::endl;
            break;
    }
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
static int fSPM(realtype t, N_Vector u, N_Vector udot, N_Vector resval,
        void *user_data)
{
    realtype cA, cB, cAlt, cBlt, cArt, cBrt;
    realtype hordA, hordB;
    realtype hordcoA, hordcoB;
    realtype *udata, *dudata, *res;
    int jx;
    UserData data;

    data   = (UserData) user_data;
    udata  = N_VGetArrayPointer(u);
    dudata = N_VGetArrayPointer(udot);
    res = N_VGetArrayPointer(resval);

    /* Loop over all grid points. */
    double potCA = calc_potCantera(p_nameCathodeSurf, IJth(udata,1,MR-1), p_Iapp*p_Rel(Tref), p_Iapp, Tref);
    double potAN = calc_potCantera(p_nameAnodeSurf, IJth(udata,2,MR-1), 0, p_Iapp, Tref);

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
        IJth(res, 1, jx) = IJth(dudata, 1, jx) - hordcoA/SUNSQR(p_rP_ca);
        IJth(res, 2, jx) = IJth(dudata, 2, jx) - hordcoB/SUNSQR(p_rP_an);
        IJth(res, 3, jx) = IJth(udata, 3, jx) - potCA;
        IJth(res, 4, jx) = IJth(udata, 4, jx) - potAN;
    }

    return(0);
}

/* f routine. Compute RHS function f(t,u). */
static int fSPMT(realtype t, N_Vector u, N_Vector udot, N_Vector resval,
        void *user_data)
{
    realtype cA, cB, cAlt, cBlt, cArt, cBrt;
    realtype hordA, hordB;
    realtype hordcoA, hordcoB;
    realtype *udata, *dudata, *res;
    realtype qIrr, qRev, qOut, T;
    realtype Uca, Uan, Ucell, Ueq, dUdTca, dUdTan;
    int jx;
    UserData data;

    data   = (UserData) user_data;
    udata  = N_VGetArrayPointer(u);
    dudata = N_VGetArrayPointer(udot);
    res = N_VGetArrayPointer(resval);

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

        IJth(res, 1, jx) = IJth(dudata, 1, jx) - hordcoA/SUNSQR(p_rP_ca);
        IJth(res, 2, jx) = IJth(dudata, 2, jx) - hordcoB/SUNSQR(p_rP_an);
        IJth(res, 3, jx) = IJth(dudata, 3, jx) - 1/(p_rhoCell*p_volCell*p_cpCell)*(qIrr+qRev-qOut);
    }
    return(0);
}
