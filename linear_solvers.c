/** @file 
 * A collection of numerical methods for solving A dparam = C
 */

#include "linear_solvers.h" 

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit.h>
 
#include "param_evolver.h"

/** when approx'ing param change by TSVD, truncate SVs smaller than below*max */
const double DEFAULT_SVD_TOLERANCE = 0.01;

/** 
 * when approx'ing the param change by Tikhonov regularisation, this decides
 * how many different Tikhonov Parameters will be tested in the search for the
 * optimal (by L-curve testing) 
 */
const double TIKHONOV_PARAM_SEARCH_SIZE = 3; // must be >= 3

/** 
 * minimum value allowed of the Tikhonov regularisation param (weighting of min param constraint), 
 * to ensure that the optimal value doesn't cause too large a change in params
 */
const double TIKHONOV_REG_MIN_PARAM = 0.0001;

/** 
 * maximum value allowed of the Tikhonov regularisation param (weighting of min param constraint), 
 * to ensure that the change in params remains accurate (keeping <E> above the true minimum)
 */
const double TIKHONOV_REG_MAX_PARAM = 0.01;

/**
 * Attempts to perform direct solving by LU decomposition, which is
 * unstable and liable to fail.
 */
int approxParamsByLUDecomp(EvolverMemory *mem) {
	
	// Simon says: APPAERNTLY THIS OFTEN PICKS ANY SOLUTION: LU comp is stored in Perm
	// Simon gets different answers for resolving
	int swaps, singular;
	gsl_linalg_LU_decomp(mem->matrA, mem->permA, &swaps);
	singular = gsl_linalg_LU_solve(mem->matrA, mem->permA, mem->vecC, mem->paramChange);
	if (singular)
		printf("Direct LU decomposition failed! Aborting...\n");
		
	return singular;
}


/**
 * Householder (orthogonal equations) solve the least-squares approximation
 */
int approxParamsByLeastSquares(EvolverMemory *mem) {

	// compute A^T A and A^T C
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, mem->matrA, mem->matrA, 0.0, mem->matrATA);
	gsl_blas_dgemv(CblasTrans, 1.0, mem->matrA, mem->vecC, 0.0, mem->vecATC);
	
	// solve A^T A paramChange = A^T C, returning success flag
	return gsl_linalg_HH_solve(mem->matrATA, mem->vecATC, mem->paramChange);
}


/**
 * LU solve the same constraints with the final var set to 0
 * (Sam and Xiao's current method)
 */
int approxParamsByRemovingVar(EvolverMemory *mem) {
		
	// Asub = shave off final col of A
	for (int i=0; i < mem->numParams; i++)
		for (int j=0; j < mem->numParams-1; j++)
			gsl_matrix_set(mem->matrASub, i, j, gsl_matrix_get(mem->matrA, i, j));
	
	// Csub = remove final elem of C
	for (int i=0; i < mem->numParams-1; i++)
		gsl_vector_set(mem->vecCSub, i, gsl_vector_get(mem->vecC, i));
	
	// solve Asub paramChangeSub = Csub
	int swaps;
	gsl_linalg_LU_decomp(mem->matrASub, mem->permASub, &swaps);
	int singular = gsl_linalg_LU_solve(mem->matrASub, mem->permASub, mem->vecCSub, mem->paramChangeSub);
	
	// copy paramChangeSub to paramChange
	gsl_vector_set(mem->paramChange, mem->numParams-1, 0);
	for (int i=0; i < mem->numParams - 1; i++)
		gsl_vector_set(mem->paramChange, i, gsl_vector_get(mem->paramChangeSub, i));
	
	// flag whether this LU solve failed
	return singular;
}


/**
 * Solve for paramChange by truncated SVD
 */
int approxParamsByTSVD(EvolverMemory *mem) { 
	
	double residSum;
	size_t singValsKept;
	return gsl_multifit_linear_tsvd(
		mem->matrA, mem->vecC, mem->svdTolerance, mem->paramChange, 
		mem->svdCovar, &residSum, &singValsKept, mem->svdSpace);
}


/**
 * Solve for paramChange by Tikhonov regularisation:
 * min ||C - A paramChange||^2 + tikhonovParam^2 || paramChange ||^2
 * The tikhonov param is estimated using the L-curve criterion, using
 * as many samples as TIKHONOV_PARAM_SEARCH_SIZE (in mem obj lengths)
 * and restricted to be >= TIKHONOV_REG_MIN_PARAM to ensure params
 * don't vary wildly
 * 
 * https://www.gnu.org/software/gsl/doc/html/lls.html
 * http://www2.compute.dtu.dk/~pcha/DIP/chap5.pdf
 * http://people.bath.ac.uk/mamamf/talks/ilas.pdf
 */
int approxParamsByTikhonov(EvolverMemory *mem) {
	
	// whether we fail to sample tikhonovParams, fail to optimise it,
	// or fail to Tikhonov solve
	int failure;
	
	// compute the SVD in the workspace (needed for L-curve and solve)
	failure = gsl_multifit_linear_svd(mem->matrA, mem->svdSpace);
	if (failure) {
		printf("Computing SVD of matrA failed in Titkhonov! Aborting...\n");
		return failure;
	}
	
	// sample the system under different regularisation params (build L-curve)
	failure = gsl_multifit_linear_lcurve(
		mem->vecC, mem->tikhonovParamSamples, 
		mem->tikhonovParamRho, mem->tikhonovParamEta, mem->svdSpace
	);
	if (failure) {
		printf("Populating the Tikhonov regularisation param space failed! Aborting...\n");
		return failure;
	}
	
	// choose the best regularisation param (hardcode 0.02 performs ok)
	size_t tikhonovParamIndex;
	double tikhonovParam;
	failure = gsl_multifit_linear_lcorner(
		mem->tikhonovParamRho, mem->tikhonovParamEta, &tikhonovParamIndex
	);
	if (failure) {
		printf("Choosing the optimal Tikhonov regularisation param (by L-curve corner) failed! Using min...\n");
		tikhonovParam = TIKHONOV_REG_MIN_PARAM;
	} else {
		tikhonovParam = gsl_vector_get(mem->tikhonovParamSamples, tikhonovParamIndex); 
	}
	
	// restrict the regularisation param from being too big or small
	if (tikhonovParam < TIKHONOV_REG_MIN_PARAM)
		tikhonovParam = TIKHONOV_REG_MIN_PARAM;
	if (tikhonovParam > TIKHONOV_REG_MAX_PARAM)
		tikhonovParam = TIKHONOV_REG_MAX_PARAM;
		
	// the error ||C - A paramChange|| can be monitored
	double residualNorm, paramChangeNorm;
	
	// perform Tikhonov regularisation (where L = identity)
	failure = gsl_multifit_linear_solve(
		tikhonovParam, mem->matrA, mem->vecC, mem->paramChange,
		&residualNorm, &paramChangeNorm, mem->svdSpace
	);
	if (failure) {
		printf("Solving Titkhonov under the best regularisation param failed! Aborting...\n");
	}
	
	return failure;
}