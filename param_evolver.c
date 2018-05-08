/** @file 
 * Wick evolves a given set of parameters under a given Hamiltonian
 */
 
#include "param_evolver.h"

#define _USE_MATH_DEFINES

#include <QuEST.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit.h>

#include "hamiltonian_builder.h"

/** when approx'ing param change by TSVD, truncate SVs smaller than below*max */
double DEFAULT_SVD_TOLERANCE = 0.01;

/** 
 * when approx'ing the param change by Tikhonov regularisation, this decides
 * how many different Tikhonov Parameters will be tested in the search for the
 * optimal (by L-curve testing) 
 */
double TIKHONOV_PARAM_SEARCH_SIZE = 3; // must be >= 3

/** 
 * minimum value allowed of the Tikhonov regularisation param (weighting of min param constraint), 
 * to ensure that the optimal value doesn't cause too large a change in params
 */
double TIKHONOV_REG_MIN_PARAM = 0.0001;

/** size of the change in parameter when approxing wavefunction derivatives */
double DERIV_STEP_SIZE = 1E-8; 

/** finite-dif first deriv coefficients of psi(x+nh) for n > 0, or -1*(that for n < 0) */
double FINITE_DIFFERENCE_COEFFS[4][4] = {
	{1/2.0},
	{2/3.0, -1/12.0},
	{3/4.0, -3/20.0, 1/60.0},
	{4/5.0, -1/5.0, 4/105.0, -1/280.0}
};
// https://en.wikipedia.org/wiki/Finite_difference_coefficient

int MAX_NUM_SAVED_STATES = 100;

double EXCITATION_OF_SAVED_STATES = 10;

/**
 * Returns the number of parameters needed to repeat the default ansatz circuit
 * numBlocks times, resulting in an equal number of gates (each with an individual parameter)
 * on each qubit
 */
int chooseDefaultNumParams(MultiQubit qubits, int numBlocks) {
	
	return 3*numBlocks*qubits.numQubits;
}


void setAnsatzInitState(EvolverMemory *mem, MultiQubit qubits) {
	
	for (long long int i=0LL; i < qubits.numAmps; i++)
		mem->initState[i] = getRealAmpEl(qubits, i) + I*getImagAmpEl(qubits,i);
}


void initState(EvolverMemory *mem, MultiQubit qubits) {
	
	for (long long int i=0LL; i < qubits.numAmps; i++) {
		qubits.stateVec.real[i] = creal(mem->initState[i]);
		qubits.stateVec.imag[i] = cimag(mem->initState[i]);		
	}
}






/**
 * Applies a default parameterised ansatz circuit to the zero state, modifying
 * the wavefunction in qubits according to the values in params. This can be passed
 * to evolveParams in lieu of a custom ansatz circuit
 */
void defaultAnsatzCircuit(EvolverMemory *mem, MultiQubit qubits, double* params, int numParams) {
	
	initState(mem, qubits);
	
	int paramInd = 0;
	
	// loop the following circuit until all params are featured
	while (paramInd < numParams) {
		
		// Rx
		for (int qb=0; qb < qubits.numQubits && paramInd < numParams; qb++)
			rotateX(qubits, qb, params[paramInd++]);
		
		// Rz
		for (int qb=0; qb < qubits.numQubits && paramInd < numParams; qb++)
			rotateZ(qubits, qb, params[paramInd++]);
		
		// C(Ry)
		for (int qb=0; qb < qubits.numQubits - 1 && paramInd < numParams; qb++)
			controlledRotateY(qubits, qb, qb+1, params[paramInd++]);
		
		if (paramInd < numParams)
			controlledRotateY(qubits, qubits.numQubits-1, 0, params[paramInd++]);
	}
}


/**
 * Gives the derivative (as a 2^N length array of complex doubles) of the
 * wavefunction produced by ansatzircuit w.r.t a given parameter in the 
 * region of the passed parameters, approximated using the accuracy-order 
 * central finite-difference formula. Deriv must be pre allocated, and is
 * modified
 */
void findDeriv(
	EvolverMemory *mem, 
	double complex *deriv, long long int length, 
	void (*ansatzCircuit)(EvolverMemory *mem, MultiQubit, double*, int),
	MultiQubit qubits, double* params, int numParams, int paramInd, int accuracy
) {
	// clear deriv
	for (long long int i=0LL; i < length; i++)
		deriv[i] = 0;
	
	// approx deriv with finite difference
	double* coeffs = FINITE_DIFFERENCE_COEFFS[accuracy - 1];
	double origParam = params[paramInd];
	
	// repeatly add c*psi(p+ndp) - c*psi(p-ndp) to deriv
	for (int step=1; step <= accuracy; step++) {
		for (int sign = -1; sign <= 1; sign+=2) {
			params[paramInd] = origParam + sign*step*DERIV_STEP_SIZE;
			ansatzCircuit(mem, qubits, params, numParams);
			for (long long int i=0LL; i < length; i++)
				deriv[i] += sign * coeffs[step-1] * (getRealAmpEl(qubits, i) + I*getImagAmpEl(qubits, i));
		}
	}
	
	// divide by the step size
	for (long long int i=0LL; i < length; i++) 
		deriv[i] /= DERIV_STEP_SIZE;
	
	// reset the original param value
	params[paramInd] = origParam;
}


/**
 * Returns the real component of the inner product between two complex vectors 
 * (of the given length). Note Re(vec1 . vec2) = Re(vec2 . vec1)
 */
double realInnerProduct(double complex* vec1, double complex* vec2, long long int length) {
	
	double prod = 0;
	for (long long int i=0LL; i < length; i++)
		prod += creal(vec1[i])*creal(vec2[i]) + cimag(vec1[i])*cimag(vec2[i]);
		
	return prod;
}


double complex innerProduct(double complex* vec1, double complex* vec2, long long int length) {
	
	double complex prod = 0;
	for (long long int i=0LL; i < length; i++)
		prod += conj(vec1[i])*vec2[i];
	
	return prod;
}
double complex innerProductOnQubits(double complex* vec1, MultiQubit qubits, long long int length) {
	
	double complex prod = 0;
	for (long long int i=0LL; i < length; i++)
		prod += conj(vec1[i])*(getRealAmpEl(qubits,i) + I*getImagAmpEl(qubits,i));
	
	return prod;
}


/**
 * Populates the A and C matrices of derivatives inside of mem, using derivatives of 
 * the wavefunction which is produced by ansatzCircuit with the given parameters,
 * approximated by finite difference of accuracy in {1, 2, 3, 4}
 * mem.matrA: A_ij = 2*Re{ d <psi(p)| dp_i . d |psi(p)> dp_j }
 * mem.vecC: C_i = -2*Re{ d <psi(p)| dp_i . H . |psi(p)> }
 * The derivs and hamilState data structures in mem are also modified
 */
void computeDerivMatrices(
	EvolverMemory *mem, void (*ansatzCircuit)(EvolverMemory *mem, MultiQubit, double*, int),
	MultiQubit qubits, double* params, Hamiltonian hamil, int accuracy) 
{	
	// collect wavef derivs w.r.t each parameter
	for (int i=0; i < mem->numParams; i++)
		findDeriv(mem, mem->derivs[i], mem->stateSize, ansatzCircuit, qubits, params, mem->numParams, i, accuracy);
	
	// populate matrix A with inner product of derivs
	for (int i=0; i < mem->numParams; i++)
		for (int j=0; j < mem->numParams; j++)
			gsl_matrix_set(mem->matrA, i, j, realInnerProduct(mem->derivs[i], mem->derivs[j], mem->stateSize));
	
	// populate vector C with inner product of derivs and H on psi
	ansatzCircuit(mem, qubits, params, mem->numParams);
	applyHamil(mem->hamilState, qubits, hamil);
	for (int i=0; i < mem->numParams; i++)
		gsl_vector_set(mem->vecC, i, -realInnerProduct(mem->derivs[i], mem->hamilState, mem->stateSize));
}


/**
 * Updates the deriv matrices to reflect excitations in the Hamiltonian for each of the
 * saved states (those recorded by exciteStateInHamiltonian)
 */
void exciteSavedStatesInDerivMatrices(EvolverMemory *mem, MultiQubit qubits) {
	
	for (int s=0; s < mem->numSavedStates; s++) {
		double complex saveProjOnQubits = innerProductOnQubits(mem->savedStates[s], qubits, qubits.numAmps);
	
		for (int i=0; i < mem->numParams; i++) {
			double complex saveProjOnDeriv = innerProduct(mem->derivs[i], mem->savedStates[s], qubits.numAmps);
			
			double currC = gsl_vector_get(mem->vecC, i);
			double newC = currC - creal(saveProjOnDeriv * saveProjOnQubits) *  EXCITATION_OF_SAVED_STATES;
			gsl_vector_set(mem->vecC, i, newC);
		}
	}
}


/**
 * adds noise to A and C matrices in mem. 
 * Each element experiences a chance in value of [-1,1]*fractionalVar * value
 */
void addNoiseToDerivMatrices(EvolverMemory *mem, double fractionalVar) {
	
	double oldval, newval, randd;
			
	for (int i=0; i < mem->numParams; i++) {
		
		oldval = gsl_vector_get(mem->vecC, i);
		randd = rand() / (double) RAND_MAX;
		newval = oldval + fractionalVar * fabs(oldval) * (2*randd - 1);
		gsl_vector_set(mem->vecC, i, newval);
		
		for (int j=0; j < mem->numParams; j++) {
			
			oldval = gsl_matrix_get(mem->matrA, i, j);
			randd = rand() / (double) RAND_MAX;
			newval = oldval + fractionalVar * fabs(oldval) * (2*randd - 1);
			gsl_matrix_set(mem->matrA, i, j, newval);
		}
	}
}


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
	
	// restrict the regularisation param from being too small (to keep ||paramChange|| small)
	if (tikhonovParam < TIKHONOV_REG_MIN_PARAM)
		tikhonovParam = TIKHONOV_REG_MIN_PARAM;
		
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


/**
 * Given a list of parameters, a parameterised wavefunction (through ansatzCircuit) and
 * a diagonal Hamiltonian, modifies the parameters by a single time-step under imaginary time 
 * evolution, using Euler's method. defaultAnsatzCircuit may be passed in lieu of a custom one.
 * Param evolution is done by repeatedly simulating a parameterised circuit and using finite-
 * difference approximations of derivatives to populate and here solve a family of
 * linear equations of the parameters. The passed function inversionMethod is called to approximate 
 * a solution to A paramChange = C, which must modify  mem->paramChange, and return 0 for success or
 * 1 for failure. In lieu of a custom method, approxParamsByLeastSquares, ...ByRemovingVar, ...ByTVSD,
 * ...ByTikhonov may be passed.
 * Percent noise in the A and C elements can be passed to emulate their inference from experiment.
 * evolveParams should update the parameters so that the parameterised wavefunction moves closer to the 
 * ground-state of the given Hamiltonian, though doesn't necessarily long-term converge to sol.
 * A return of evolveOutcome:SUCCESS indicates  numerical updating of the params worked,
 * and FAILED indicates the inversionMethod failed.
 * evolveParams finally updates the wavefunction in qubits under the new parameters, so that you can
 * monitor wavefunction properties over evolution.
 * mem contains memory for matrices and arrays which are modified
  * @param mem						data structures created with prepareEvolverMemory()
  * @param ansatzCircuit			parameterised circuit to apply to qubits which generates wavefunc(params)
  * @param inversionMethod function to solve/update paramChange when direct LU solution fails
  * @param qubits					QuEST qubits instance
  * @param params					list of current values of the parameters, to be updated
  * @param hamil					Hamiltonian built by hamiltonian_builder
  * @param timeStepSize			size of the step-size in imag-time
  * @param wrapParams				1 to keep params in [0, 2pi) by wrap-around, 0 to let them grow
  * @param derivAccuracy			accuracy of finite-difference approx to param derivs in {1, 2, 3, 4}
  * @param matrNoise				noise (in [0, 1]) to add to A and C matrices before solving. each elem += +- noise*val
  * @return SUCCESS 				indicates numerical updating (by inversionMethod) of the params worked 
  * @return FAILED					indicates inversionMethod failed
  */
evolveOutcome evolveParams(
	EvolverMemory *mem, 
	void (*ansatzCircuit)(EvolverMemory *mem, MultiQubit, double*, int), 
	int (*inversionMethod)(EvolverMemory*),
	MultiQubit qubits, double* params, Hamiltonian hamil, double timeStepSize, 
	int wrapParams, int derivAccuracy, double matrNoise) 
{
	evolveOutcome outcome = SUCCESS;
	
	// compute matrices A and C
	computeDerivMatrices(mem, ansatzCircuit, qubits, params, hamil, derivAccuracy);
	
	// morph C to excite Hamiltonain states
	exciteSavedStatesInDerivMatrices(mem, qubits);
	
	// add a little noise to A and C
	addNoiseToDerivMatrices(mem, matrNoise);
	
	// solve A paramChange = C
	int singular = inversionMethod(mem);
	if (singular)
		return FAILED;

	// update params
	for (int i=0; i < mem->numParams; i++)
		params[i] += timeStepSize*gsl_vector_get(mem->paramChange, i);
	
	// wrap-around params in [0, 2pi] to avoid overflow
	if (wrapParams) {
		for (int i=0; i < mem->numParams; i++)
			params[i] = fmod(params[i], 2*M_PI);
	}
	
	// update the wavefunction with new params
	ansatzCircuit(mem, qubits, params, mem->numParams);
	
	// indicate params updated successfully
	return outcome;
}


/**
 * Behaves similarly to evolveParams, but using gradient descent (disregards A matrix)
 * and cannot numerically fail (besides repeated use not converging to a solution).
 */
void evolveParamsByGradientDescent(
	EvolverMemory *mem, 
	void (*ansatzCircuit)(EvolverMemory *mem, MultiQubit, double*, int), 
	MultiQubit qubits, double* params, Hamiltonian hamil, double timeStepSize, int wrapParams,
	int derivAccuracy) 
{
	// compute matrices A and C
	computeDerivMatrices(mem, ansatzCircuit, qubits, params, hamil, derivAccuracy);
	
	// update params
	for (int i=0; i < mem->numParams; i++)
		params[i] += timeStepSize*gsl_vector_get(mem->vecC, i);
	
	// wrap-around params in [0, 2pi] to avoid overflow
	if (wrapParams) {
		for (int i=0; i < mem->numParams; i++)
			params[i] = fmod(params[i], 2*M_PI);
	}
	
	ansatzCircuit(mem, qubits, params, mem->numParams);
}


/**
 * Allocates memory for the data structures needed by the evolveParams function,
 * for a given number of parameters. The memory should be kept for the lifetime of the
 * evolver simulation, and afterward freed with freeEvolverMemory
 */
EvolverMemory prepareEvolverMemory(MultiQubit qubits, int numParams) {
	
	// disable errors, so we can detect/handle when A is singular
	gsl_set_error_handler_off();
	
	// allocate arrays
	EvolverMemory memory;
	memory.numParams = numParams;
	memory.stateSize = pow(2LL, qubits.numQubits);
	memory.initState = malloc(memory.stateSize * sizeof *memory.initState);
	initStateZero(&qubits);
	setAnsatzInitState(&memory, qubits);
	memory.hamilState = malloc(memory.stateSize * sizeof *memory.hamilState);
	memory.derivs = malloc(memory.numParams * sizeof *memory.derivs);
	for (int i=0; i < memory.numParams; i++)
		memory.derivs[i] = malloc(memory.stateSize * sizeof **memory.derivs);
		
	// allocate saved state objects
	memory.numSavedStates = 0;
	memory.savedStates = malloc(MAX_NUM_SAVED_STATES * sizeof *memory.savedStates);
	
	// allocate LU-decomp objects
	memory.matrA = gsl_matrix_alloc(memory.numParams, memory.numParams);
	memory.permA = gsl_permutation_alloc(memory.numParams);
	memory.vecC = gsl_vector_alloc(memory.numParams);
	memory.paramChange = gsl_vector_alloc(memory.numParams);
	
	// allocate least-squares approx objects
	memory.matrATA = gsl_matrix_alloc(memory.numParams, memory.numParams);
	memory.vecATC = gsl_vector_alloc(memory.numParams);
	
	// allocate column-removal objects
	memory.matrASub = gsl_matrix_alloc(memory.numParams, memory.numParams - 1);
	memory.vecCSub = gsl_vector_alloc(memory.numParams - 1);
	memory.permASub = gsl_permutation_alloc(memory.numParams);
	memory.paramChangeSub = gsl_vector_alloc(memory.numParams - 1);
	
	// allocate TSVD objects
	memory.svdTolerance = DEFAULT_SVD_TOLERANCE;
	memory.svdSpace = gsl_multifit_linear_alloc(memory.numParams, memory.numParams);
	memory.svdCovar = gsl_matrix_alloc(memory.numParams, memory.numParams);
	
	// allocate Tikhonov regularisation objects (uses SVD space)
	//memory.tikhonovParam = DEFAULT_TIKHONOV_PARAM;
	memory.tikhonovParamSamples = gsl_vector_alloc(TIKHONOV_PARAM_SEARCH_SIZE);
	memory.tikhonovParamRho = gsl_vector_alloc(TIKHONOV_PARAM_SEARCH_SIZE);
	memory.tikhonovParamEta = gsl_vector_alloc(TIKHONOV_PARAM_SEARCH_SIZE);
	memory.tikhonovVecL = gsl_vector_alloc(memory.numParams);
	for (int i=0; i < memory.numParams; i++)
		gsl_vector_set(memory.tikhonovVecL, i, 1);
	
	return memory;
}


void exciteStateInHamiltonian(EvolverMemory *memory, MultiQubit qubits) {
	
	// abort if we've run out of space
	// TODO: don't abort, just re-allocate a bigger-sized array
	if (memory->numSavedStates == MAX_NUM_SAVED_STATES) {
		printf("WARNING! The maximum number of states (%d) have been excited in the Hamiltonian.\n"
			   "This and future calls to exciteStateInHamiltonian do nothing!\n"
			   "Remedy by increasing MAX_NUM_SAVED_STATES in param_evolver.c\n", MAX_NUM_SAVED_STATES);
		return;
	}
	
	// copy the state in qubits
	double complex *state = malloc(qubits.numAmps * sizeof *state);
	for (long long int i=0LL; i < qubits.numAmps; i++)
		state[i] = getRealAmpEl(qubits, i) + I*getImagAmpEl(qubits, i);
	
	// save it 
	memory->savedStates[(memory->numSavedStates)++] = state;
}


void clearExcitedStates(EvolverMemory *memory) {
	
	for (int i=0; i < memory->numSavedStates; i++)
		free(memory->savedStates[i]);
	memory->numSavedStates = 0;
}


/**
 * Frees the memory allocated with prepareEvolverMemory. 
 * This should be done when there are no more calls to evolveParams
 */
void freeEvolverMemory(EvolverMemory *memory) {
	
	// free state info
	for (int i=0; i < memory->numParams; i++)
		free(memory->derivs[i]);
	
	// free arrays
	free(memory->derivs);
	free(memory->hamilState);
	free(memory->initState);
	
	// free saved states
	for (int i=0; i < memory->numSavedStates; i++)
		free(memory->savedStates[i]);
	free(memory->savedStates);
	
	// free LU-decomp structures
	gsl_matrix_free(memory->matrA);
	gsl_vector_free(memory->vecC);
	gsl_vector_free(memory->paramChange);
	gsl_permutation_free(memory->permA);
	
	// free least-squares structures
	gsl_matrix_free(memory->matrATA);
	gsl_vector_free(memory->vecATC);
	
	// free column-removal structures
	gsl_matrix_free(memory->matrASub);
	gsl_vector_free(memory->vecCSub);
	gsl_permutation_free(memory->permASub);
	gsl_vector_free(memory->paramChangeSub);
	
	// free TSVD structures
	gsl_multifit_linear_free(memory->svdSpace);
	gsl_matrix_free(memory->svdCovar);
	
	// free Tikhonov structures
	gsl_vector_free(memory->tikhonovParamSamples);
	gsl_vector_free(memory->tikhonovParamRho);
	gsl_vector_free(memory->tikhonovParamEta);
	gsl_vector_free(memory->tikhonovVecL);
}




