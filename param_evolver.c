/** @file 
 * Wick evolves a given set of parameters under a given Hamiltonian
 */
 
#include "param_evolver.h"

#define _USE_MATH_DEFINES

#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit.h>

#include "QuEST/qubits.h"

/** greatest change in parameter allowed in a single evolve step before flagged */
double MAX_PARAM_CHANGE = 1.572; // approx pi/2
// NOT USED ATM 

/** when approx'ing ill-posed param change by TSVD, truncate SVs smaller than below*max */
double DEFAULT_SVD_TOLERANCE = 0.001;

/** size of the change in parameter when approxing wavefunction derivatives */
double DERIV_STEP_SIZE = 1E-5; 

/** finite-dif first deriv coefficients of psi(x+nh) for n > 0, or -1*(that for n < 0) */
double FINITE_DIFFERENCE_COEFFS[4][4] = {
	{1/2.0},
	{2/3.0, -1/12.0},
	{3/4.0, -3/20.0, 1/60.0},
	{4/5.0, -1/5.0, 4/105.0, -1/280.0}
};
// https://en.wikipedia.org/wiki/Finite_difference_coefficient


/**
 * Returns the number of parameters needed to repeat the default ansatz circuit
 * numBlocks times, resulting in an equal number of gates (each with an individual parameter)
 * on each qubit
 */
int chooseDefaultNumParams(MultiQubit qubits, int numBlocks) {
	
	return 3*numBlocks*qubits.numQubits;
}


/**
 * Applies a default parameterised ansatz circuit to the zero state, modifying
 * the wavefunction in qubits according to the values in params. This can be passed
 * to evolveParams in lieu of a custom ansatz circuit
 */
void defaultAnsatzCircuit(MultiQubit qubits, double* params, int numParams) {
	
	// initStatePlus makes Ax=b very ill-posed! Why?!
	initStateZero(&qubits);
	
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
	
	/**
	 * START MOSTLY IN GROUND AND SOME EXCITED STATES TO SEE WHY THE EIGENSTATES ARE ATTRACTING
	 * 
	 */
	 
	 
	/**
	 *	PLOT SPECTRUM EVOLUTION (to show getting stuck in first excited)
	 */
}


/**
 * Gives the derivative (as a 2^N length array of complex doubles) of the
 * wavefunction produced by ansatzircuit w.r.t a given parameter in the 
 * region of the passed parameters, approximated using the accuracy-order 
 * central finite-difference formula. Deriv must be pre allocated, and is
 * modified
 */
void findDeriv(
	double complex *deriv, long long int length, 
	void (*ansatzCircuit)(MultiQubit, double*, int),
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
			ansatzCircuit(qubits, params, numParams);
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


/**
 * Takes whatever wavefunction qubits currently has, and gives that due to applying
 * the given diagonal Hamiltonian (by modifying hamilState). 
 * Does not modify the wavefunction of qubits.
 */
void applyDiagHamiltonian(
	double complex* hamilState, long long int length, 
	MultiQubit qubits, double* diagHamiltonian
) {
	for (long long int i=0LL; i < length; i++)
		hamilState[i] = diagHamiltonian[i]*(getRealAmpEl(qubits, i) + I*getImagAmpEl(qubits, i));
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
	evolverMemory *mem, void (*ansatzCircuit)(MultiQubit, double*, int),
	MultiQubit qubits, double* params, double* diagHamiltonian, int accuracy) 
{	
	// collect wavef derivs w.r.t each parameter
	for (int i=0; i < mem->numParams; i++)
		findDeriv(mem->derivs[i], mem->stateSize, ansatzCircuit, qubits, params, mem->numParams, i, accuracy);
	
	// populate matrix A with inner product of derivs
	for (int i=0; i < mem->numParams; i++)
		for (int j=0; j < mem->numParams; j++)
			gsl_matrix_set(mem->matrA, i, j, realInnerProduct(mem->derivs[i], mem->derivs[j], mem->stateSize));
	
	// populate vector C with inner product of derivs and H on psi
	ansatzCircuit(qubits, params, mem->numParams);
	applyDiagHamiltonian(mem->hamilState, mem->stateSize, qubits, diagHamiltonian);
	for (int i=0; i < mem->numParams; i++)
		gsl_vector_set(mem->vecC, i, -realInnerProduct(mem->derivs[i], mem->hamilState, mem->stateSize));
}


/**
 * adds noise to A and C matrices in mem. 
 * Each element experiences += +- fractionalVar * value, where the sign is randomly chosen
 */
void addNoiseToDerivMatrices(evolverMemory *mem, double fractionalVar) {
	
	int sign;
	double oldval, newval;
			
	for (int i=0; i < mem->numParams; i++) {
		
		sign = -1 + 2*(rand() < RAND_MAX/2.0);
		oldval = gsl_vector_get(mem->vecC, i);
		newval = (1 + sign*fractionalVar)*oldval;
		gsl_vector_set(mem->vecC, i, newval);
		
		for (int j=0; j < mem->numParams; j++) {
			
			sign = -1 + 2*(rand() < RAND_MAX/2.0);
			oldval = gsl_matrix_get(mem->matrA, i, j);
			newval = (1 + sign*fractionalVar)*oldval;
			gsl_matrix_set(mem->matrA, i, j, newval);
		}
	}
}


/**
 * Householder (orthogonal equations) solve the least-squares approximation
 */
int approxParamsByLeastSquares(evolverMemory *mem) {
	
	printf("Using Least Squares (Householder)\n");

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
int approxParamsByRemovingVar(evolverMemory *mem) {
	
	printf("Using Var Removal\n");
		
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
int approxParamsByTSVD(evolverMemory *mem) { 
	
	printf("Using TVSD\n");
	
	double residSum;
	size_t singValsKept;
	return gsl_multifit_linear_tsvd(
		mem->matrA, mem->vecC, mem->svdTolerance, mem->paramChange, 
		mem->svdCovar, &residSum, &singValsKept, mem->svdSpace);
}


/**
 * Given a list of parameters, a parameterised wavefunction (through ansatzCircuit) and
 * a diagonal Hamiltonian, modifies the parameters by a single time-step under imaginary time 
 * evolution, using Euler's method. defaultAnsatzCircuit may be passed in lieu of a custom one.
 * Param evolution is done by repeatedly simulating a parameterised circuit and using finite-
 * difference approximations of derivatives to populate and here solve a family of
 * linear equations of the parameters. If the numerical solving fails (the matrices are ill-conditioned),
 * the passed illPosedRecoveryMethod is called to approximate a solution, which must modify 
 * mem->paramChange; approxParamsByLeastSquares, ...ByRemovingVar, ...ByTVSD may be passed.
 * This function should update the parameters so that the parameterised wavefunction moves closer to the 
 * ground-state of the given Hamiltonian, though don't necessarily long-term converge to sol.
 * A return of evolveOutcome:SUCCESS indicates direct numerical updating of the params worked,
 * while RECOVERED indicates the 
 * and FAILED indicates the illPosedRecoveryMethod also numerically failed.
 * Updates the wavefunction in qubits under the new parameters.
 * mem contains memory for matrices and arrays which are modified
  * @param mem						data structures created with prepareEvolverMemory()
  * @param ansatzCircuit			parameterised circuit to apply to qubits which generates wavefunc(params)
  * @param illPosedRecoveryMethod function to solve/update paramChange when direct LU solution fails
  * @param qubits					QuEST qubits instance
  * @param params					list of current values of the parameters, to be updated
  * @param diagHamiltonian			the diagonal terms of the Hamiltonian, under which to imag-time evolve
  * @param timeStepSize			size of the step-size in imag-time
  * @param wrapParams				1 to keep params in [0, 2pi) by wrap-around, 0 to let them grow
  * @param derivAccuracy			accuracy of finite-difference approx to param derivs in {1, 2, 3, 4}
  * @param matrNoise				noise (in [0, 1]) to add to A and C matrices before solving. each elem += +- noise*val
  * @return SUCCESS 				indicates direct numerical updating (by LU decomposition) of the params worked 
  * @return RECOVERED				indicaets illPosedRecoveryMethod had to be used (LU failed) but was successful
  * @return FAILED					indicates direct LU solving and illPosedRecoveryMethod failed
  */
evolveOutcome evolveParams(
	evolverMemory *mem, 
	void (*ansatzCircuit)(MultiQubit, double*, int), 
	int (*illPosedRecoveryMethod)(evolverMemory*),
	MultiQubit qubits, double* params, double* diagHamiltonian, double timeStepSize, 
	int wrapParams, int derivAccuracy, double matrNoise) 
{
	evolveOutcome outcome = SUCCESS;
	
	// compute matrices A and C
	computeDerivMatrices(mem, ansatzCircuit, qubits, params, diagHamiltonian, derivAccuracy);
	
	// add a little noise to A and C
	addNoiseToDerivMatrices(mem, matrNoise);
	
	// solve A paramChange = C
	int swaps, singular;
	gsl_linalg_LU_decomp(mem->matrA, mem->permA, &swaps);
	singular = gsl_linalg_LU_solve(mem->matrA, mem->permA, mem->vecC, mem->paramChange);

	// if that failed, try an approximation
	if (singular) {
		outcome = RECOVERED;
		singular = illPosedRecoveryMethod(mem);
	}

	// if that failed, give up
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
	ansatzCircuit(qubits, params, mem->numParams);
	
	// indicate params updated successfully
	return outcome;
}


/**
 * Behaves similarly to evolveParams, but using gradient descent (disregards A matrix)
 * and cannot numerically fail (besides repeated use not converging to a solution).
 */
void evolveParamsByGradientDescent(
	evolverMemory *mem, void (*ansatzCircuit)(MultiQubit, double*, int), 
	MultiQubit qubits, double* params, double* diagHamiltonian, double timeStepSize, int wrapParams,
	int derivAccuracy) 
{
	// compute matrices A and C
	computeDerivMatrices(mem, ansatzCircuit, qubits, params, diagHamiltonian, derivAccuracy);
	
	// update params
	for (int i=0; i < mem->numParams; i++)
		params[i] += timeStepSize*gsl_vector_get(mem->vecC, i);
	
	// wrap-around params in [0, 2pi] to avoid overflow
	if (wrapParams) {
		for (int i=0; i < mem->numParams; i++)
			params[i] = fmod(params[i], 2*M_PI);
	}
	
	ansatzCircuit(qubits, params, mem->numParams);
}


/**
 * Allocates memory for the data structures needed by the evolveParams function,
 * for a given number of parameters. The memory should be kept for the lifetime of the
 * evolver simulation, and afterward freed with freeEvolverMemory
 */
evolverMemory prepareEvolverMemory(MultiQubit qubits, int numParams) {
	
	// disable errors, so we can detect/handle when A is singular
	gsl_set_error_handler_off();
	
	// allocate arrays
	evolverMemory memory;
	memory.numParams = numParams;
	memory.stateSize = pow(2LL, qubits.numQubits);
	memory.hamilState = malloc(memory.stateSize * sizeof *memory.hamilState);
	memory.derivs = malloc(memory.numParams * sizeof *memory.derivs);
	for (int i=0; i < memory.numParams; i++)
		memory.derivs[i] = malloc(memory.stateSize * sizeof **memory.derivs);
	
	// allocate gsl objects
	memory.matrA = gsl_matrix_alloc(memory.numParams, memory.numParams);
	memory.permA = gsl_permutation_alloc(memory.numParams);
	memory.vecC = gsl_vector_alloc(memory.numParams);
	memory.paramChange = gsl_vector_alloc(memory.numParams);
	
	// allocate ill-posed gsl objects
	memory.matrATA = gsl_matrix_alloc(memory.numParams, memory.numParams);
	memory.vecATC = gsl_vector_alloc(memory.numParams);
	
	memory.matrASub = gsl_matrix_alloc(memory.numParams, memory.numParams - 1);
	memory.vecCSub = gsl_vector_alloc(memory.numParams - 1);
	memory.permASub = gsl_permutation_alloc(memory.numParams);
	memory.paramChangeSub = gsl_vector_alloc(memory.numParams - 1);
	
	memory.svdTolerance = DEFAULT_SVD_TOLERANCE;
	memory.svdSpace = gsl_multifit_linear_alloc(memory.numParams, memory.numParams);
	memory.svdCovar = gsl_matrix_alloc(memory.numParams, memory.numParams);
	
	return memory;
}


/**
 * Frees the memory allocated with prepareEvolverMemory. 
 * This should be done when there are no more calls to evolveParams
 */
void freeEvolverMemory(evolverMemory *memory) {
	
	// free state info
	for (int i=0; i < memory->numParams; i++)
		free(memory->derivs[i]);
		
	free(memory->derivs);
	free(memory->hamilState);
	
	// free invertible-A structures
	gsl_matrix_free(memory->matrA);
	gsl_permutation_free(memory->permA);
	gsl_vector_free(memory->vecC);
	gsl_vector_free(memory->paramChange);
	
	// free ill-posed A structures
	gsl_matrix_free(memory->matrATA);
	gsl_vector_free(memory->vecATC);
	
	gsl_matrix_free(memory->matrASub);
	gsl_vector_free(memory->vecCSub);
	gsl_permutation_free(memory->permASub);
	gsl_vector_free(memory->paramChangeSub);
	
	gsl_multifit_linear_free(memory->svdSpace);
	gsl_matrix_free(memory->svdCovar);
}

