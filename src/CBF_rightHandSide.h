#ifndef CBF_RIGHTHANDSIDE_H
#define CBF_RIGHTHANDSIDE_H
#include <math.h>
#include <nvector/nvector_serial.h>
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_types.h>

#include "CBF_driver.h"
#include "CBF_parameters.h"

static int CBF_rightHandSide(realtype t, N_Vector y, N_Vector ydot, void* user_data);

#endif