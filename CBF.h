#ifndef CBF_H
#define CBF_H

#include <cvodes/cvodes_dense.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
// #include <time.h>

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
#define SIM_TIME 100.0 // this is total number of beats in a simulation run

#define isAFib 0 // Atrial Fibrillation is switched off

#endif // CBF_H