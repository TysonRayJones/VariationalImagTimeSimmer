#ifndef PARAM_EVOLVER_H_
#define PARAM_EVOLVER_H_

#include <QuEST.h>
#include <complex.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multifit.h>

#include "hamiltonian_builder.h"


/** number of contiguous iterations considered when checking convergence */
extern const int NUM_ITERS_IN_STUCK_CHECK;

/** size of the change in parameter when approxing wavefunction derivatives */
extern const double DERIV_STEP_SIZE;

extern const int MAX_NUM_SAVED_STATES;
extern const double EXCITATION_OF_SAVED_STATES;





/**
 * A container for the memory used by evolveParams. Is to be created once
 * via prepareEvolverMemory and passed to calls to evolveParams,
 * and finally freed by freeEvolverMemory
 */
typedef struct {
	
	// state and hamiltonian info
	long long int stateSize;
	int numParams;
	double complex **firstDerivs;
    double complex **secondDerivs;
    double complex **mixedDerivs;
	double complex *hamilState;
	double complex *initState;
	
	// saved spectrum states (for suguru's idea)
	int numSavedStates;
	double complex **savedStates;
	
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
	
} EvolverMemory;


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
evolveOutcome evolveParamsByImaginaryTime(
	EvolverMemory *mem, 
	void (*ansatzCircuit)(EvolverMemory *mem, MultiQubit, double*, int), 
	int (*inversionMethod)(EvolverMemory*),
	MultiQubit qubits, double* params, Hamiltonian hamil,
	double timeStepSize, int wrapParams, int derivAccuracy, 
    int shotNoiseNumSamplesA, int shotNoiseNumSamplesC, double decoherenceFactor
);
	
/**
 * Behaves similarly to evolveParams, but using gradient descent (disregards A matrix)
 * and cannot numerically fail (besides repeated use not converging to a solution).
 * @return SUCCESS		always
 */
evolveOutcome evolveParamsByGradientDescent(
	EvolverMemory *mem, void (*ansatzCircuit)(EvolverMemory *mem, MultiQubit, double*, int), 
	MultiQubit qubits, double* params, Hamiltonian hamil, 
	double timeStepSize, int wrapParams, int derivAccuracy,
    int shotNoiseNumSamplesA, int shotNoiseNumSamplesC, double decoherenceFactor
);
    
    
/* adds noise */
typedef void (*NoiseModel)(EvolverMemory*, double);    
void addDumbNoiseToDerivMatrices(EvolverMemory *mem, double fractionalVar);
void addShotNoiseToDerivMatrices(EvolverMemory *mem, double numSamples);
    
	
/**
 * returns whether the current simulation has halted, based on the evolution of the parameters
 * evolution is stopped if, for each of the last NUM_ITERS_IN_STUCK_CHECK iterations, 
 * the sum of (absolute) changes of the parameters is less than MAX_PARAM_CHANGE_WHEN_STUCK
 */
int isStuck(double*** paramEvo, int simIteration, int numParams, int iteration, double paramChangeThreshold);
//int isStuck(double* paramEvo, int simIteration, int numParams, int numIterations, int iteration);

/**
 * Allocates memory for the data structures needed by the evolveParams function,
 * for a given number of blocks. The memory should be kept for the lifetime of the
 * evolver simulation, and afterward freed with freeEvolverMemory
 */
EvolverMemory prepareEvolverMemory(MultiQubit qubits, int numBlocks);

/**
 * Copies the state of qubits to the list of 'saved states', which are excited
 * in the Hamiltonian in each variational iteration
 */
void exciteStateInHamiltonian(EvolverMemory *memory, MultiQubit qubits);

void clearExcitedStates(EvolverMemory *memory);

/**
 * Frees the memory allocated with prepareEvolverMemory. 
 * This should be done when there are no more calls to evolveParams
 */
void freeEvolverMemory(EvolverMemory *memory);



/** needed by variational to project state on spectrum */
double complex innerProductOnQubits(double complex* vec1, MultiQubit qubits, long long int length);


#endif // PARAM_EVOLVER_H_