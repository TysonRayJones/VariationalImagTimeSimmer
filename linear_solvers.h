#ifndef LINEAR_SOLVERS_H_
#define LINEAR_SOLVERS_H_

#include "param_evolver.h"

/** when approx'ing param change by TSVD, truncate SVs smaller than below*max */
extern const double DEFAULT_SVD_TOLERANCE;

/** 
 * when approx'ing the param change by Tikhonov regularisation, this decides
 * how many different Tikhonov Parameters will be tested in the search for the
 * optimal (by L-curve testing) 
 */
extern const double TIKHONOV_PARAM_SEARCH_SIZE;

/** 
 * minimum value allowed of the Tikhonov regularisation param (weighting of min param constraint), 
 * to ensure that the optimal value doesn't cause too large a change in params
 */
extern const double TIKHONOV_REG_MIN_PARAM;

/** 
 * maximum value allowed of the Tikhonov regularisation param (weighting of min param constraint), 
 * to ensure that the change in params remains accurate (keeping <E> above the true minimum)
 */
extern const double TIKHONOV_REG_MAX_PARAM;


/**
 * provided methods for numerically solving for the change in params (pass to evolveParams)
 * Individiaul descriptions are in param_evolver.c
 */
int approxParamsByLUDecomp(EvolverMemory *mem);
int approxParamsByLeastSquares(EvolverMemory *mem);
int approxParamsByRemovingVar(EvolverMemory *mem);
int approxParamsByTSVD(EvolverMemory *mem);
int approxParamsByTikhonov(EvolverMemory *mem);

#endif // LINEAR_SOLVERS_H_