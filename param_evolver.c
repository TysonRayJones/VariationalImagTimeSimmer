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
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multifit.h>

#include "hamiltonian_builder.h"
#include "ansatz_circuits.h"
#include "linear_solvers.h"

/** number of contiguous iterations considered when checking convergence */
const int NUM_ITERS_IN_STUCK_CHECK = 2;

/** size of the change in parameter when approxing wavefunction derivatives */
const double DERIV_STEP_SIZE = 1E-5; // 1E-8; 

const int MAX_NUM_SAVED_STATES = 1000;
const double EXCITATION_OF_SAVED_STATES = 10;



/** finite-dif first deriv coefficients of psi(x+nh) for n > 0, or -1*(that for n < 0) */
double FIRST_DERIV_FINITE_DIFFERENCE_COEFFS[4][4] = {
	{1/2.0},
	{2/3.0, -1/12.0},
	{3/4.0, -3/20.0, 1/60.0},
	{4/5.0, -1/5.0, 4/105.0, -1/280.0}
};
// https://en.wikipedia.org/wiki/Finite_difference_coefficient


/* psi(x+nh) for n >= 0 (same for n < 0) */
double SECOND_DERIV_FINITE_DIFFERENCE_COEFFS[4][5] = {
    {-2.0, 1.0},
    {-5/2.0, 4/3.0, -1/12.0},
    {-49/18.0, 3/2.0, -3/20.0, 1/90.0},
    {-205/72.0, 8/5.0, -1/5.0, 8/315.0, -1/560.0}
};
// https://en.wikipedia.org/wiki/Finite_difference_coefficient


/** dpsi/dp
 * Gives the derivative (as a 2^N length array of complex doubles) of the
 * wavefunction produced by ansatzircuit w.r.t a given parameter in the 
 * region of the passed parameters, approximated using the accuracy-order 
 * central finite-difference formula. Deriv must be pre allocated, and is
 * modified
 */
void findFirstDeriv(
	EvolverMemory *mem, 
	double complex *deriv, long long int length, 
	void (*ansatzCircuit)(EvolverMemory *mem, MultiQubit, double*, int),
	MultiQubit qubits, double* params, int numParams, int paramInd, int accuracy
) {
	// clear deriv
	for (long long int i=0LL; i < length; i++)
		deriv[i] = 0;
	
	// approx deriv with finite difference
	double* coeffs = FIRST_DERIV_FINITE_DIFFERENCE_COEFFS[accuracy - 1];
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

/** d^2 psi/dp^2
 * Gives the 2nd derivative (as a 2^N length array of complex doubles) of the
 * wavefunction produced by ansatzircuit w.r.t a given parameter in the 
 * region of the passed parameters, approximated using the accuracy-order 
 * central finite-difference formula. Deriv must be pre allocated, and is
 * modified
 */
void findSecondDeriv(
	EvolverMemory *mem, 
	double complex *secondDeriv, long long int length, 
	void (*ansatzCircuit)(EvolverMemory *mem, MultiQubit, double*, int),
	MultiQubit qubits, double* params, int numParams, int paramInd, int accuracy
) {

	// clear deriv
	for (long long int i=0LL; i < length; i++)
		secondDeriv[i] = 0;

	// approx deriv with finite difference
	double* coeffs = SECOND_DERIV_FINITE_DIFFERENCE_COEFFS[accuracy - 1];
	double origParam = params[paramInd];
    
    for (int step=0; step <= accuracy; step++) {
		
        params[paramInd] = origParam + step*DERIV_STEP_SIZE;
        ansatzCircuit(mem, qubits, params, numParams);
		for (long long int i=0LL; i < length; i++)
			secondDeriv[i] += coeffs[step] * (getRealAmpEl(qubits, i) + I*getImagAmpEl(qubits, i));
        
        if (step == 0)
            continue;
        
        params[paramInd] = origParam - step*DERIV_STEP_SIZE;
        ansatzCircuit(mem, qubits, params, numParams);
		for (long long int i=0LL; i < length; i++)
			secondDeriv[i] += coeffs[step] * (getRealAmpEl(qubits, i) + I*getImagAmpEl(qubits, i));
        
    }
    
	// divide by the step size squared
	for (long long int i=0LL; i < length; i++) 
		secondDeriv[i] /= (DERIV_STEP_SIZE * DERIV_STEP_SIZE);
	
	// reset the original param value
	params[paramInd] = origParam;
}

/** d^2 psi/dp1 dp2
 */
void findMixedDeriv(
	EvolverMemory *mem, 
	double complex *mixedDeriv, long long int length, 
	void (*ansatzCircuit)(EvolverMemory *mem, MultiQubit, double*, int),
	MultiQubit qubits, double* params, int numParams, 
    int paramInd1, int paramInd2, int accuracy)
{
	// clear deriv
	for (long long int i=0LL; i < length; i++)
		mixedDeriv[i] = 0;
	
	// approx deriv with finite difference
	double* coeffs = FIRST_DERIV_FINITE_DIFFERENCE_COEFFS[accuracy - 1];
	double origParam1 = params[paramInd1];
    double origParam2 = params[paramInd2];
	
	// repeatly add c*psi(p+ndp) - c*psi(p-ndp) to deriv
	for (int step1=1; step1 <= accuracy; step1++) {
		for (int sign1 = -1; sign1 <= 1; sign1+=2) {
			params[paramInd1] = origParam1 + sign1*step1*DERIV_STEP_SIZE;
            
        	for (int step2=1; step2 <= accuracy; step2++) {
        		for (int sign2 = -1; sign2 <= 1; sign2+=2) {
        			params[paramInd2] = origParam2 + sign2*step2*DERIV_STEP_SIZE;
            
			        ansatzCircuit(mem, qubits, params, numParams);
			        for (long long int i=0LL; i < length; i++)
				        mixedDeriv[i] += sign1*coeffs[step1-1] * sign2*coeffs[step2-1] * (getRealAmpEl(qubits, i) + I*getImagAmpEl(qubits, i));
                }
            }
		}
	}
	
	// divide by the step size
	for (long long int i=0LL; i < length; i++) 
		mixedDeriv[i] /= (DERIV_STEP_SIZE * DERIV_STEP_SIZE);
	
	// reset the original param values
	params[paramInd1] = origParam1; 
	params[paramInd2] = origParam2;  
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
void computeFirstDerivs(
	EvolverMemory *mem, void (*ansatzCircuit)(EvolverMemory *mem, MultiQubit, double*, int),
	MultiQubit qubits, double* params, Hamiltonian hamil, int accuracy) 
{	
	// collect wavef derivs w.r.t each parameter
	for (int i=0; i < mem->numParams; i++)
		findFirstDeriv(mem, mem->firstDerivs[i], mem->stateSize, ansatzCircuit, qubits, params, mem->numParams, i, accuracy);
    
    // also compute hamilState (used by vecC and matrHessian)
	ansatzCircuit(mem, qubits, params, mem->numParams);
	applyHamil(mem->hamilState, qubits, hamil);
}

void computeSecondDerivs(
	EvolverMemory *mem, void (*ansatzCircuit)(EvolverMemory *mem, MultiQubit, double*, int),
	MultiQubit qubits, double* params, Hamiltonian hamil, int accuracy)
{
	// collect wavef 2nd derivs w.r.t each parameter (d^2 psi / dp^2)
	for (int i=0; i < mem->numParams; i++)
		findSecondDeriv(mem, mem->secondDerivs[i], mem->stateSize, ansatzCircuit, qubits, params, mem->numParams, i, accuracy);
    
    // collect mixed derivs (d^2 psi /dp1 dp2)
    int ind=0;
    for (int i=0; i < mem->numParams; i++)
        for (int j=i+1; j < i; j++)
            findMixedDeriv(mem, mem->mixedDerivs[ind++], mem->stateSize, ansatzCircuit, qubits, params, mem->numParams, i, j, accuracy);
}

/** to be called after computeFirstDerivs */
void computeVectorA(EvolverMemory *mem) 
{	
	for (int i=0; i < mem->numParams; i++)
		for (int j=0; j < mem->numParams; j++)
			gsl_matrix_set(mem->matrA, i, j, realInnerProduct(mem->firstDerivs[i], mem->firstDerivs[j], mem->stateSize));
}

/** to be called after computeFirstDerivs */
void computeVectorC(
	EvolverMemory *mem, void (*ansatzCircuit)(EvolverMemory *mem, MultiQubit, double*, int),
	MultiQubit qubits, double* params, Hamiltonian hamil) 
{	
	for (int i=0; i < mem->numParams; i++)
		gsl_vector_set(mem->vecC, i, -realInnerProduct(mem->firstDerivs[i], mem->hamilState, mem->stateSize));
}

void computeHessian(EvolverMemory *mem, void (*ansatzCircuit)(EvolverMemory *mem, MultiQubit, double*, int),
	MultiQubit qubits, double* params, Hamiltonian hamil, int accuracy) {
    
    // populated mem.secondDerivs and mem.mixedDerivs
    computeSecondDerivs(mem, ansatzCircuit, qubits, params, hamil, accuracy);
    
    // set diagonals to second derivs
    for (int i=0; i < mem->numParams; i++)
        gsl_matrix_set(mem->matrHessian, i, i, 
            realInnerProduct(mem->secondDerivs[i], mem->hamilState, mem->stateSize));
    
    // set off-diagonals to mixed derivs
    int ind=0;
    for (int i=0; i < mem->numParams; i++) {
        for (int j=i+1; j < i; j++) {
            REAL innerp = realInnerProduct(mem->mixedDerivs[ind++], mem->hamilState, mem->stateSize);
            gsl_matrix_set(mem->matrHessian, i, j, innerp); 
            gsl_matrix_set(mem->matrHessian, j, i, innerp);
        }
    }
}

void computeAdamMoments() {
    
}




/**
 * Updates the deriv matrices to reflect excitations in the Hamiltonian for each of the
 * saved states (those recorded by exciteStateInHamiltonian), by modifying the C vector
 */
void exciteSavedStatesInDerivMatrices(EvolverMemory *mem, MultiQubit qubits) {
	
	for (int s=0; s < mem->numSavedStates; s++) {
		double complex saveProjOnQubits = innerProductOnQubits(mem->savedStates[s], qubits, qubits.numAmps);
	
		for (int i=0; i < mem->numParams; i++) {
			double complex saveProjOnDeriv = innerProduct(mem->firstDerivs[i], mem->savedStates[s], qubits.numAmps);
			
			double currC = gsl_vector_get(mem->vecC, i);
			double newC = currC - creal(saveProjOnDeriv * saveProjOnQubits) * EXCITATION_OF_SAVED_STATES;
			gsl_vector_set(mem->vecC, i, newC);
		}
	}
}





double sampleNormalDistrib(double mean, double var) {
    
    // sample the the standard normal distribution via Box-Muller
    double unif1 = rand() / (double) RAND_MAX;
    double unif2 = rand() / (double) RAND_MAX;
    double norm1 = sqrt(-2.0 * log(unif1)) * cos(2*M_PI*unif2);
    
    // sample a normal dist with new expected value and variance
	return mean + norm1*var;
}

/** Effects decoherence in the circuits and shot noise in the deriv sampling by changing the matrix
 * elements to be random samples of normal distributions. Overestimates controlled-gate shot-noise to be
 * the same as single-qubit-gate variance. Overestimates the shot noise in the C vector by assuming
 * every hamil term has inner product zero (in order to maximise the measurement variance). Effects
 * decoherence by shrinking every expected measurement value closer to 0 (the mixed state).
 */
void addNoiseToVectorC(
    EvolverMemory *mem, double chemHamilCoeffSquaredSum, 
    int shotNoiseNumSamplesC, double decoherFac
) {
	double oldval, mesvar, newval;
			
	for (int i=0; i < mem->numParams; i++) {
		
        // update C
		oldval = gsl_vector_get(mem->vecC, i);
        if (shotNoiseNumSamplesC == 0)
            mesvar = 0;
        else
            mesvar = chemHamilCoeffSquaredSum/4.0/shotNoiseNumSamplesC;
		newval = sampleNormalDistrib(decoherFac * oldval, mesvar);
		gsl_vector_set(mem->vecC, i, newval);
        
    }
}
void addNoiseToMatrixA(
    EvolverMemory *mem, double chemHamilCoeffSquaredSum, 
    int shotNoiseNumSamplesA, double decoherFac
) {    
	double oldval, mesvar, newval;	
	for (int i=0; i < mem->numParams; i++) {
		for (int j=0; j < mem->numParams; j++) {
            
            // current element
            oldval = gsl_matrix_get(mem->matrA, i, j);
            if (shotNoiseNumSamplesA == 0)
                mesvar = 0;
            else
                mesvar = (1/16.0 - decoherFac*decoherFac * oldval*oldval)/(double)shotNoiseNumSamplesA;
            newval = sampleNormalDistrib(decoherFac * oldval, mesvar);
			gsl_matrix_set(mem->matrA, i, j, newval);
		}
	}
}
void addNoiseToHessian(
    EvolverMemory *mem, double chemHamilCoeffSquaredSum, 
    int shotNoiseNumSamplesHess, double decoherFac
) {    
	double oldval, mesvar, newval;	
	for (int i=0; i < mem->numParams; i++) {
		for (int j=0; j < mem->numParams; j++) {
            
            // current element
            oldval = gsl_matrix_get(mem->matrHessian, i, j);
            if (shotNoiseNumSamplesHess == 0)
                mesvar = 0;
            else
                mesvar = (1/16.0 - decoherFac*decoherFac * oldval*oldval)/(double)shotNoiseNumSamplesHess;
            newval = sampleNormalDistrib(decoherFac * oldval, mesvar);
			gsl_matrix_set(mem->matrHessian, i, j, newval);
		}
	}
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
  * @param shotNoiseNumSamplesA     Number of measurements per A element, used to add shot noise
  * @param shotNoiseNumSamplesC     Number of measurements per Hamiltonian term per C element, used for shot noise
  * @param decoherenceFactor        Multiplies with expected values of measurements to effect decoherence
  * @return SUCCESS 				indicates numerical updating (by inversionMethod) of the params worked 
  * @return FAILED					indicates inversionMethod failed
  */
evolveOutcome evolveParamsByImaginaryTime(
	EvolverMemory *mem, 
	void (*ansatzCircuit)(EvolverMemory*, MultiQubit, double*, int), 
	int (*inversionMethod)(EvolverMemory*, gsl_matrix*),
	MultiQubit qubits, double* params, Hamiltonian hamil, double timeStepSize, 
	int wrapParams, int derivAccuracy, 
    int shotNoiseNumSamplesA, int shotNoiseNumSamplesC, double decoherenceFactor) 
{
	evolveOutcome outcome = SUCCESS;
	
	// compute matrices A and C
	computeFirstDerivs(mem, ansatzCircuit, qubits, params, hamil, derivAccuracy);
    computeVectorA(mem);
    computeVectorC(mem, ansatzCircuit, qubits, params, hamil);
	
	// morph C to excite Hamiltonain states
	exciteSavedStatesInDerivMatrices(mem, qubits);
	
	// add a little noise to A and C
    addNoiseToMatrixA(mem, hamil.termCoeffSquaredSum, shotNoiseNumSamplesA, decoherenceFactor);
    addNoiseToVectorC(mem, hamil.termCoeffSquaredSum, shotNoiseNumSamplesC, decoherenceFactor);
    
	// solve A paramChange = C
	int singular = inversionMethod(mem, mem->matrA);
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
evolveOutcome evolveParamsByGradientDescent(
	EvolverMemory *mem, 
	void (*ansatzCircuit)(EvolverMemory *mem, MultiQubit, double*, int), 
	MultiQubit qubits, double* params, Hamiltonian hamil, double timeStepSize, int wrapParams,
	int derivAccuracy, 
    int shotNoiseNumSamplesA, int shotNoiseNumSamplesC, double decoherenceFactor)
{
	// compute vector C
	computeFirstDerivs(mem, ansatzCircuit, qubits, params, hamil, derivAccuracy);
    computeVectorC(mem, ansatzCircuit, qubits, params, hamil);
	
	// morph C to excite Hamiltonain states
	exciteSavedStatesInDerivMatrices(mem, qubits);
	
	// add a little noise to C
    addNoiseToVectorC(mem, hamil.termCoeffSquaredSum, shotNoiseNumSamplesC, decoherenceFactor);
    
	// update params
	for (int i=0; i < mem->numParams; i++)
		params[i] += timeStepSize*gsl_vector_get(mem->vecC, i);
	
	// wrap-around params in [0, 2pi] to avoid overflow
	if (wrapParams) {
		for (int i=0; i < mem->numParams; i++)
			params[i] = fmod(params[i], 2*M_PI);
	}
	
	// update the wavefunction with new params
	ansatzCircuit(mem, qubits, params, mem->numParams);
	
	return SUCCESS;
}

evolveOutcome evolveParamsByHessian(
	EvolverMemory *mem, 
	void (*ansatzCircuit)(EvolverMemory*, MultiQubit, double*, int), 
	int (*inversionMethod)(EvolverMemory*, gsl_matrix*),
	MultiQubit qubits, double* params, Hamiltonian hamil, double timeStepSize, 
	int wrapParams, int derivAccuracy, 
    int shotNoiseNumSamplesHess, int shotNoiseNumSamplesC, double decoherenceFactor) 
{
	evolveOutcome outcome = SUCCESS;
	
	// compute matrices A and C
	computeFirstDerivs(mem, ansatzCircuit, qubits, params, hamil, derivAccuracy);
    computeVectorC(mem, ansatzCircuit, qubits, params, hamil);
    computeHessian(mem, ansatzCircuit, qubits, params, hamil, derivAccuracy);
	
	/* EXCITING STATES IN HESSIAN METHOD WILL REQUIRE MODIFYING HESSIAN MATRIX */
	
	// add a little noise to C and the hessian matrix
    addNoiseToHessian(mem, hamil.termCoeffSquaredSum, shotNoiseNumSamplesHess, decoherenceFactor);
    addNoiseToVectorC(mem, hamil.termCoeffSquaredSum, shotNoiseNumSamplesC, decoherenceFactor);
    
	// solve -Hessian paramChange = C
	int singular = inversionMethod(mem, mem->matrHessian); // multiply by -1???
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
 * returns whether the current simulation has halted, based on the evolution of the parameters
 * evolution is stopped if, for each of the last NUM_ITERS_IN_STUCK_CHECK iterations, 
 * the sum of (absolute) changes of the parameters is less than MAX_PARAM_CHANGE_WHEN_STUCK
 */
int isStuck(double*** paramEvo, int simIteration, int numParams, int iteration, double paramChangeThreshold) {

	// if too few iterations have happened, we're not stuck
	if (iteration < NUM_ITERS_IN_STUCK_CHECK + 1)
		return 0;
		
	// paramEvo[sims][params][iter]
	for (int i=iteration-NUM_ITERS_IN_STUCK_CHECK; i < iteration; i++) {
		
		// calculate |dp1| + |dp2| + .... for this iteration
		double paramVecChange = 0;
		
		for (int n=0; n < numParams; n++) {
			double currVal = paramEvo[simIteration][n][iteration];
			double prevVal = paramEvo[simIteration][n][iteration - 1];
			paramVecChange += fabs(currVal - prevVal);
		}
		
		// not stuck if params change enough in any iteration in the considered window
		if (paramVecChange > paramChangeThreshold)
			return 0;
	}
	
	// paramChange is small in all considered iterations: we're stuck!
	return 1;
}
 
 
// this is for an array paramEvo, not pointers
/*
int isStuck(double* paramEvo, int simIteration, int numParams, int numIterations, int iteration) {
	
	// if too few iterations have happened, we're not stuck
	if (iteration < NUM_ITERS_IN_STUCK_CHECK + 1)
		return 0;
	
	// start of this simulation's memory block in paramEvo
	int currIterInd = simIteration*numParams*numIterations;
	
	// require paramChange is small for every of these iterations
	for (int i=iteration-NUM_ITERS_IN_STUCK_CHECK; i < iteration; i++) {
	
		// calculate |dp1| + |dp2| + .... for this iteration
		double paramVecChange = 0;
		
		for (int n=0; n < numParams; n++) {
			double currVal = paramEvo[currIterInd + i + n*numIterations];
			double prevVal = paramEvo[currIterInd + i + n*numIterations - 1]; 
			paramVecChange += fabs(currVal - prevVal);			
		}
		
		// not stuck if params change enough in any iteration in the considered window
		if (paramVecChange > MAX_PARAM_CHANGE_WHEN_STUCK)
			return 0;
	}
	
	// paramChange is small in all considered iterations: we're stuck!
	return 1;
}
*/


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
	memory.firstDerivs = malloc(memory.numParams * sizeof *memory.firstDerivs);
	memory.secondDerivs = malloc(memory.numParams * sizeof *memory.secondDerivs);
    memory.mixedDerivs = malloc(numParams*(numParams-1) * sizeof *memory.mixedDerivs);
	for (int i=0; i < memory.numParams; i++) {
		memory.firstDerivs[i] = malloc(memory.stateSize * sizeof **memory.firstDerivs);
		memory.secondDerivs[i] = malloc(memory.stateSize * sizeof **memory.secondDerivs);
    }
    for (int i=0; i < numParams*(numParams-1); i++) {
        memory.mixedDerivs[i] = malloc(memory.stateSize * sizeof **memory.mixedDerivs);
    }
		
	// allocate saved state objects
	memory.numSavedStates = 0;
	memory.savedStates = malloc(MAX_NUM_SAVED_STATES * sizeof *memory.savedStates);
	
	// allocate LU-decomp objects
	memory.matrA = gsl_matrix_alloc(memory.numParams, memory.numParams);
	memory.permA = gsl_permutation_alloc(memory.numParams);
	memory.vecC = gsl_vector_alloc(memory.numParams);
	memory.paramChange = gsl_vector_alloc(memory.numParams);
    memory.matrHessian = gsl_matrix_alloc(memory.numParams, memory.numParams);
	
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
	for (int i=0; i < memory->numParams; i++) {
		free(memory->firstDerivs[i]);
        free(memory->secondDerivs[i]);
    }
    for (int i=0; i < (memory->numParams)*(memory->numParams - 1); i++) {
        free(memory->mixedDerivs[i]);
    }
	
	// free arrays
	free(memory->firstDerivs);
    free(memory->secondDerivs);
    free(memory->mixedDerivs);
	free(memory->hamilState);
	free(memory->initState);
	
	// free saved states
	for (int i=0; i < memory->numSavedStates; i++)
		free(memory->savedStates[i]);
	free(memory->savedStates);
	
	// free LU-decomp structures
	gsl_matrix_free(memory->matrA);
    gsl_matrix_free(memory->matrHessian);
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




