#ifndef CBF_DRIVER_H
#define CBF_DRIVER_H

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_dense.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>
#include <time.h>

#include "CBF_parameters.h"
#include "CBF_rightHandSide.h"
// #include "ursino.h"

#define Ith(v, i) NV_Ith_S(v, i - 1) /* Ith numbers components 1..NEQ    */
#define IJth(A, i, j) DENSE_ELEM(A, i - 1, j - 1) /* IJth numbers rows,cols 1..NEQ    */
#define RTOL RCONST(1.0e-3) /* scalar relative tolerance            */
#define ATOL RCONST(1.0e-6) /* scalar absolute tolerance components */
#define MAXSTEPS 500000

#define NEQ 57 /* number of equations, ODEs */
#define N_PARAMETER 172 // number of parameters
#define DELTAT 0.01 // time step, seconds. You cannot have time step bigger than a few milliseconds as there are fast processes going on.
#define debg 1 // this is always 1.

// a check to choose if we want to output all or few state variables
// 1 means all, 0 means only important pressures.
#define numBeats 20 // this is total number of beats in a simulation run

int CBF_driver(void* data);

static int check_retval(void* returnvalue, const char* funcname, int opt)
{
    int* retval;

    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    if (opt == 0 && returnvalue == NULL) {
        fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
        return (1);
    }

    /* Check if retval < 0 */
    else if (opt == 1) {
        retval = (int*)returnvalue;
        if (*retval < 0) {
            fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
                funcname, *retval);
            return (1);
        }
    }

    /* Check if function returned NULL pointer - no memory allocated */
    else if (opt == 2 && returnvalue == NULL) {
        fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
        return (1);
    }

    return (0);
}

#endif