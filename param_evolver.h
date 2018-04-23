#ifndef PARAM_EVOLVER_H_
#define PARAM_EVOLVER_H_

#include <complex.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multifit.h>

#include "QuEST/qubits.h"


/**
 * Returns the number of parameters needed to repeat the default ansatz circuit
 * numBlocks times, resulting in an equal number of gates (each with an individual parameter)
 * on each qubit
 */
int chooseDefaultNumParams(MultiQubit qubits, int numBlocks);

/**
 * Applies a default parameterised ansatz circuit to the zero state, modifying
 * the wavefunction in qubits according to the values in params. This can be passed
 * to evolveParams in lieu of a custom ansatz circuit
 */
void defaultAnsatzCircuit(MultiQubit qubits, double* params, int numParams);

/**
 * A container for the memory used by evolveParams. Is to be created once
 * via prepareEvolverMemory and passed to calls to evolveParams,
 * and finally freed by freeEvolverMemory
 */
typedef struct {
	
	// state and hamiltonian info
	long long int stateSize;
	int numParams;
	double complex **derivs;
	double complex *hamilState;
	
	// gsl objects for invertible A
	gsl_matrix *matrA;
	gsl_permutation *permA;
	gsl_vector *vecC, *paramChange;
	
	// gsl objects for least squares, when A not invertible
	gsl_matrix *matrATA;
	gsl_vector *vecATC;
	
	// gsl objects for removing a variable when A not invertible
	gsl_matrix *matrASub;
	gsl_vector *vecCSub;
	gsl_permutation *permASub;
	gsl_vector *paramChangeSub;
	
	// gsl objects for TSVD when A not invertible
	double svdTolerance;
	gsl_multifit_linear_workspace *svdSpace;
	gsl_matrix *svdCovar;
	
	//double tikhonovParam;
	gsl_vector *tikhonovParamSamples;
	gsl_vector *tikhonovParamRho;
	gsl_vector *tikhonovParamEta;
	gsl_vector *tikhonovVecL;
	
} evolverMemory;

/**
 * provided methods for numerically solving for the change in params (pass to evolveParams)
 * Individiaul descriptions are in param_evolver.c
 */
int approxParamsByLUDecomp(evolverMemory *mem);
int approxParamsByLeastSquares(evolverMemory *mem);
int approxParamsByRemovingVar(evolverMemory *mem);
int approxParamsByTSVD(evolverMemory *mem);
int approxParamsByTikhonov(evolverMemory *mem);

/** flags which indicate the success of inversionMethod passed to evolveParams */
typedef enum evolveOutcome {SUCCESS, FAILED} evolveOutcome;

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
  * @param diagHamiltonian			the diagonal terms of the Hamiltonian, under which to imag-time evolve
  * @param timeStepSize			size of the step-size in imag-time
  * @param wrapParams				1 to keep params in [0, 2pi) by wrap-around, 0 to let them grow
  * @param derivAccuracy			accuracy of finite-difference approx to param derivs in {1, 2, 3, 4}
  * @param matrNoise				noise (in [0, 1]) to add to A and C matrices before solving. each elem += +- noise*val
  * @return SUCCESS 				indicates numerical updating (by inversionMethod) of the params worked 
  * @return FAILED					indicates inversionMethod failed
  */
evolveOutcome evolveParams(
	evolverMemory *mem, 
	void (*ansatzCircuit)(MultiQubit, double*, int), 
	int (*inversionMethod)(evolverMemory*),
	MultiQubit qubits, double* params, double* diagHamiltonian, 
	double timeStepSize, int wrapParams, int derivAccuracy, double matrNoise);
	
/**
 * Behaves similarly to evolveParams, but using gradient descent (disregards A matrix)
 * and cannot numerically fail (besides repeated use not converging to a solution).
 */
void evolveParamsByGradientDescent(
	evolverMemory *mem, void (*ansatzCircuit)(MultiQubit, double*, int), 
	MultiQubit qubits, double* params, double* diagHamiltonian, double timeStepSize, int wrapParams,
	int derivAccuracy);

/**
 * Allocates memory for the data structures needed by the evolveParams function,
 * for a given number of blocks. The memory should be kept for the lifetime of the
 * evolver simulation, and afterward freed with freeEvolverMemory
 */
evolverMemory prepareEvolverMemory(MultiQubit qubits, int numBlocks);

/**
 * Frees the memory allocated with prepareEvolverMemory. 
 * This should be done when there are no more calls to evolveParams
 */
void freeEvolverMemory(evolverMemory *memory);


#endif // PARAM_EVOLVER_H_