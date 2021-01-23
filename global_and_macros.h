#pragma once

//Author: 
//Date of last modification : 25/10/2020


//For debugging over console
//#define ON_SCREEN


#define MAX_DOUBLE 1e100
#define MAX_INT	100000000
#define DELTA 1e-2


//Fixed CPLEX parameters
#define PARALLEL 1	//CPX_PARAM_PARALLEL
#define THREADS 1	//CPX_PARAM_THREADS
#define EPAGAP 1e-6	//CPA_PARAM_EPAGAP
#define EPINT 1e-5	//CPX_PARAM_EPINT
#define DELTA_NADIR	1 // In the calculation of the hypervolume, the reference point is Nadir + DELTA_NADIR
#define HV_FONSECA			//Use algorithm of Fonseca & Lopez-Ibanez & Paquete to calculate the hypervolume. If disabled this macro, the program uses WFG algorithm.



#ifndef HV_FONSECA
//#define HV_WFG
#endif // !


#ifdef HV_FONSECA
//#include "hv.c"
#endif

#ifdef HV_WFG
//#include "wfg.h"
#endif




//GLOBAL VARIABLES
#ifdef ON_SCREEN
int active_boxes = 0;
int box_id = 0;
#endif
bool warning = false;

