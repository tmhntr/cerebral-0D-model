#ifndef CBF_RHS_H
#define CBF_RHS_H

#include <math.h>

#include "CBF.h"
#include "CBF_parameters.h"

int RHS(realtype t, N_Vector y, N_Vector ydot, void* user_data);

#endif // CBF_RHS_H