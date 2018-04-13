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

#include "../QuEST/qubits.h"


/** size of the change in parameter when approxing wavefunction derivatives */
double DERIV_STEP_SIZE = 1.0/1000000.0;

/** greatest change in parameter allowed in a single evolve step before flagged */
double MAX_PARAM_CHANGE = 1.572; // approx pi/2


/**
 * Returns the number of parameters needed to repeat the default ancilla circuit
 * numBlocks times, resulting in an equal number of gates (each with an individual parameter)
 * on each qubit
 */
int chooseDefaultNumParams(MultiQubit qubits, int numBlocks) {
	
	return 3*numBlocks*qubits.numQubits;
}


/*
void applyAncillaCircuit(MultiQubit qubits, double* params, int numParams) {
	
	double t = params[0];
	double p = params[1];
	
	double norm = sqrt(
		pow(sin(2*t), 2) + 
		pow(sin(2*p), 2) +
		pow(cos(t + 2*p), 2) +
		pow(sin((t + p)/2), 2)
	);
	
	for (int i=0; i < 4; i++)
		qubits.stateVec.imag[i] = 0;
		
	qubits.stateVec.real[0] = sin((t+p)/2.0)/norm;
	qubits.stateVec.real[1] = sin(2*p)/norm;
	qubits.stateVec.real[2] = cos(t + 2*p)/norm;
	qubits.stateVec.real[3] = sin(2*t)/norm;
}
*/

/*
void applyAncillaCircuit(MultiQubit qubits, double* params, int numParams) {
	
	double norm = sqrt(1 + pow(sin(2*params[0]), 2));
	
	for (int i=0; i < 4; i++)
		qubits.stateVec.imag[i] = 0;
		
	qubits.stateVec.real[0] = sin(params[0])/norm;
	qubits.stateVec.real[1] = 0;
	qubits.stateVec.real[2] = cos(params[0])/norm;
	qubits.stateVec.real[3] = sin(2*params[0])/norm;
	
	double normCheck = 0;
	for (int j=0; j < 4; j++)
		normCheck += getProbEl(qubits, j);
}
*/


/**
 * Applies a default parameterised ancilla circuit to the zero state, modifying
 * the wavefunction in qubits according to the values in params. This can be passed
 * to evolveParams in lieu of a custom ancilla circuit
 */
void defaultAncillaCircuit(MultiQubit qubits, double* params, int numParams) {
	
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
}


/**
 * Gives the derivative (as a 2^N length array of complex doubles) of the
 * wavefunction produced by ancillaCircuit w.r.t a given parameter in the 
 * region of the passed parameters, approximated using the first-order 
 * central finite-difference formula. Deriv must be pre allocated, and is
 * modified
 */
void findDeriv(
	double complex *deriv, long long int length, 
	void (*ancillaCircuit)(MultiQubit, double*, int),
	MultiQubit qubits, double* params, int numParams, int paramInd
) {
	// psi(p + h/2)
	params[paramInd] += DERIV_STEP_SIZE/2;
	ancillaCircuit(qubits, params, numParams);
	for (long long int i=0LL; i < length; i++)
		deriv[i] = getRealAmpEl(qubits, i) + I*getImagAmpEl(qubits, i);

	// - psi(p - h/2)
	params[paramInd] -= DERIV_STEP_SIZE;
	ancillaCircuit(qubits, params, numParams);
	for (long long int i=0LL; i < length; i++)
		deriv[i] -= getRealAmpEl(qubits, i) + I*getImagAmpEl(qubits, i);
	
	// / h
	for (long long int i=0LL; i < length; i++)
		deriv[i] /= DERIV_STEP_SIZE;
		
	// reset p
	params[paramInd] += DERIV_STEP_SIZE/2;
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
 * the wavefunction which is produced by ancillaCircuit with the given parameters. 
 * mem.matrA: A_ij = 2*Re{ d <psi(p)| dp_i . d |psi(p)> dp_j }
 * mem.vecC: C_i = -2*Re{ d <psi(p)| dp_i . H . |psi(p)> }
 * The derivs and hamilState data structures in mem are also modified
 */
void computeDerivMatrices(
	evolverMemory *mem, void (*ancillaCircuit)(MultiQubit, double*, int),
	MultiQubit qubits, double* params, double* diagHamiltonian) 
{
	// collect wavef derivs w.r.t each parameter
	for (int i=0; i < mem->numParams; i++)
		findDeriv(mem->derivs[i], mem->stateSize, ancillaCircuit, qubits, params, mem->numParams, i);
	
	// populate matrix A with inner product of derivs
	for (int i=0; i < mem->numParams; i++)
		for (int j=0; j < mem->numParams; j++)
			gsl_matrix_set(mem->matrA, i, j, realInnerProduct(mem->derivs[i], mem->derivs[j], mem->stateSize));
	
	// populate vector C with inner product of derivs and H on psi
	ancillaCircuit(qubits, params, mem->numParams);
	applyDiagHamiltonian(mem->hamilState, mem->stateSize, qubits, diagHamiltonian);
	for (int i=0; i < mem->numParams; i++)
		gsl_vector_set(mem->vecC, i, -realInnerProduct(mem->derivs[i], mem->hamilState, mem->stateSize));
}


/*
void evolveParamsByGD(
	evolverMemory *mem,
	MultiQubit qubits, double* params, double* diagHamiltonian, double timeStepSize) 
{
	// compute matrices A and C
	computeDerivMatrices(mem, qubits, params, diagHamiltonian);

	for (int i=0; i < mem->numParams; i++) 
		params[i] += fmod(timeStepSize*gsl_vector_get(mem->vecC, i), 2*M_PI);
		 
	applyAncillaCircuit(qubits, params, mem->numParams);
}
*/


/**
 * Given a list of parameters, a parameterised wavefunction (through ancillaCircuit) and
 * a diagonal Hamiltonian, modifies the parameters by a single time-step under imaginary time 
 * evolution, using Euler's method. defaultAncillaCircuit may be passed in lieu of a custom one.
 * Param evolution is done by repeatedly simulating a parameterised circuit and using finite-
 * difference approximations of derivatives to populate and here solve a family of
 * linear equations of the parameters. This function should update the parameters so that
 * the parameterised wavefunction moves closer to the ground-state of the given Hamiltonian.
 * The output is 0 if updating the parameters is successful, otherwise a 1 is returned which 
 * indicates timeStepSize was too large and a parameter experienced too great a change,
 * and none of the parameters have been updated.
 * Updates the wavefunction in qubits under the new parameters.
 * mem contains memory for matrices and arrays which are modified
 */
int evolveParams(
	evolverMemory *mem, void (*ancillaCircuit)(MultiQubit, double*, int),
	MultiQubit qubits, double* params, double* diagHamiltonian, double timeStepSize) 
{
	// compute matrices A and C
	computeDerivMatrices(mem, ancillaCircuit, qubits, params, diagHamiltonian);
	
	// // solve A paramChange = C
	int swaps, singular;
	gsl_linalg_LU_decomp(mem->matrA, mem->permA, &swaps);
	singular = gsl_linalg_LU_solve(mem->matrA, mem->permA, mem->vecC, mem->paramChange);
	
	if (singular) {
		
		/*
		 * LEAST SQUARES METHOD: increased energy drastically
		 * 
		printf("A was singular. Using least squares:\n");
		// compute A^T A
		gsl_matrix *matrTranspAA = gsl_matrix_alloc(mem->numParams, mem->numParams);
		gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, mem->matrA, mem->matrA, 0.0, matrTranspAA);
		
		// compute A^T C
		gsl_vector *vecTranspAC = gsl_vector_alloc(mem->numParams);
		gsl_blas_dgemv(CblasTrans, 1.0, mem->matrA, mem->vecC, 0.0, vecTranspAC);
		
		// solve A^T A paramChange = A^T C
		gsl_linalg_LU_decomp(matrTranspAA, mem->permA, &swaps);
		singular = gsl_linalg_LU_solve(matrTranspAA, mem->permA, vecTranspAC, mem->paramChange);
		
		printf("Reattempted singular: %d\n", singular);
		gsl_matrix_free(matrTranspAA);
		gsl_vector_free(vecTranspAC);
		*/
		

		
		/*
		 * LEAST SQUARES METHOD with ORTHOGONAL EQUATIONS (Householder)
		 * works and is stable, but energy can still increase (and params leave 2PI)
		 */
		
		printf("A was singular. Using least squares with householder:\n");
		// compute A^T A
		gsl_matrix *matrTranspAA = gsl_matrix_alloc(mem->numParams, mem->numParams);
		gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, mem->matrA, mem->matrA, 0.0, matrTranspAA);
		
		// compute A^T C
		gsl_vector *vecTranspAC = gsl_vector_alloc(mem->numParams);
		gsl_blas_dgemv(CblasTrans, 1.0, mem->matrA, mem->vecC, 0.0, vecTranspAC);
		
		// solve A^T A paramChange = A^T C
		//gsl_linalg_LU_decomp(matrTranspAA, mem->permA, &swaps);
		//singular = gsl_linalg_LU_solve(matrTranspAA, mem->permA, vecTranspAC, mem->paramChange);
		singular = gsl_linalg_HH_solve(matrTranspAA, vecTranspAC, mem->paramChange);
		
		printf("Reattempted singular: %d\n", singular);
		gsl_matrix_free(matrTranspAA);
		gsl_vector_free(vecTranspAC);
		
		
		
		
		/**
		 * TRY:
		 * when underdetermined, pick a sol such that some very small params go exactly to 0
		*/
		
		
		
		
		/*
		 * SAM/XIAO SHAVING OFF ROW/COL: increases energy mildly, and often needs to be recursed
		 */
		 
		 		// WHY IS PARAM CHANGE EXTREME FOR NON-SINGULAR MATRICES WHEN STARTING IN |+>???
		/*
		printf("A was singular. Shaving off rows/columns blindly.\n");
		
		// Asub = shave off final row/col of A
		gsl_matrix *matrAsub = gsl_matrix_alloc(mem->numParams - 1, mem->numParams - 1);
		for (int i=0; i < mem->numParams-1; i++)
			for (int j=0; j < mem->numParams-1; j++)
				gsl_matrix_set(matrAsub, i, j, gsl_matrix_get(mem->matrA, i, j));
		
		// Csub = remove final elem of C
		gsl_vector *vecCsub = gsl_vector_alloc(mem->numParams -1);
		for (int i=0; i < mem->numParams-1; i++)
			gsl_vector_set(vecCsub, i, gsl_vector_get(mem->vecC, i));
		
		gsl_permutation *permAsub = gsl_permutation_alloc(mem->numParams - 1);
		gsl_linalg_LU_decomp(matrAsub, permAsub, &swaps);
		
		gsl_vector *paramChangeSub = gsl_vector_alloc(mem->numParams - 1);
		singular = gsl_linalg_LU_solve(matrAsub, permAsub, vecCsub, paramChangeSub);
		
		if (singular)
			printf("The sub matrix was still singular! Skipping...\n");
		printf("Pre updated:\n");
		for (int i=0; i < mem->numParams; i++)
			printf("p%d: %lf\n", i, params[i]);
		
		for (int i=0; i < mem->numParams - 1; i++)
			params[i] += timeStepSize * gsl_vector_get(paramChangeSub, i);
			
		printf("Post updated:\n");
		for (int i=0; i < mem->numParams; i++)
			printf("p%d: %lf\n", i, params[i]);
		
		gsl_matrix_free(matrAsub);
		gsl_vector_free(vecCsub);
		gsl_permutation_free(permAsub);
		gsl_vector_free(paramChangeSub);
		*/
		
		
		
		/* HOUSE HOLDER SOLVE ORIGINAL */
		
		//printf("Lin solve failed: using Householder transformation\n");
		//singular = gsl_linalg_HH_solve(mem->matrA, mem->vecC, mem->paramChange);
		
		if (singular)
			printf("ahhhh it too was singular! Skipping...\n");
		else {
			for (int i=0; i < mem->numParams; i++) {
				printf("change in p%d:\t%lf\n", i, timeStepSize*gsl_vector_get(mem->paramChange,i));
				params[i] += timeStepSize*gsl_vector_get(mem->paramChange,i);
			}
		}
		
		
	} else {
		
		// if a param changes too drastically, don't update and report error
		/*
		for (int i=0; i < mem->numParams; i++)
			if (timeStepSize*gsl_vector_get(mem->paramChange, i) > MAX_PARAM_CHANGE)
				return 1;
		*/
		
		// update params
		for (int i=0; i < mem->numParams; i++)
			params[i] += timeStepSize*gsl_vector_get(mem->paramChange, i);
	}
	
	
	// wrap-around params in [0, 2pi] to avoid overflow
	for (int i=0; i < mem->numParams; i++)
		params[i] = fmod(params[i], 2*M_PI);
	
	// update the wavefunction with new params
	ancillaCircuit(qubits, params, mem->numParams);
	
	// indicate params updated successfully
	return 0;
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
	
	return memory;
}


/**
 * Frees the memory allocated with prepareEvolverMemory. 
 * This should be done when there are no more calls to evolveParams
 */
void freeEvolverMemory(evolverMemory *memory) {
	
	for (int i=0; i < memory->numParams; i++)
		free(memory->derivs[i]);
		
	free(memory->derivs);
	free(memory->hamilState);
	
	gsl_matrix_free(memory->matrA);
	gsl_permutation_free(memory->permA);
	gsl_vector_free(memory->vecC);
	gsl_vector_free(memory->paramChange);
}

