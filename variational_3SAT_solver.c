/** @file 
 * Solves 3SAT problems by Wick rotation
 */

#include <stdio.h>
#include <stdlib.h>

#include "hamiltonian_builder.h"	// for building and applying Hamiltonians
#include "param_evolver.h"			// for variational imag-time simulation
#include "true_evolver.h"			// for verifying variational sim
#include "mmaformatter.h"			// for outputting results to mathematica

#include "QuEST/qubits.h"			// for simulating the ansatz circuit


#define OUTPUT_FILE "wickSATdata.txt"



double getStableTimeStep(double *hamil, long long int stateSize) {
	
	double maxEigVal = 0;
	for (long long int i=0LL; i < stateSize; i++)
		if (hamil[i] > maxEigVal)
			maxEigVal = hamil[i];
			
	return 1/maxEigVal;
}


/**
 * Initialises a random wavefunction, uniform on the unit hypersphere
 * http://www.qetlab.com/RandomStateVector (k=0)
 * https://en.wikipedia.org/wiki/Normal_distribution#Generating_values_from_normal_distribution
 */
void initStateRandom(MultiQubit *qubits) {
	
	double norm_sq = 0;
	
	// generate unnormalised, normal real and imag components
	for (long long int i=0LL; i < qubits->numAmps; i++) {
		
		double x_uni = rand()/(double) RAND_MAX;
		double y_uni = rand()/(double) RAND_MAX;
		
		// generate pair of normal numbers via Box-Muller
		double x_norm = sqrt(- 2 * log(x_uni)) * cos(2 * M_PI * y_uni);
		double y_norm = sqrt(- 2 * log(x_uni)) * sin(2 * M_PI * y_uni);
		
		qubits->stateVec.real[i] = x_norm;
		qubits->stateVec.imag[i] = y_norm;
		
		norm_sq += pow(x_norm,2) + pow(y_norm,2);
	}
	
	// normalise components
	double scaling = 1/sqrt(norm_sq);
	for (long long int i=0LL; i < qubits->numAmps; i++) {
		qubits->stateVec.real[i] *= scaling;
		qubits->stateVec.imag[i] *= scaling;
	}
}



void getSpectrum(
	double *hamil, long long int hamilSize, 
	double** spectrum, int** degeneracy, int** stateToSpecMap, long long int* spectrumSize) {
	
	// collect all energies, pops in roomy arrays
	*spectrumSize = 0;
	double energiesContainer[hamilSize];
	int populationContainer[hamilSize];
	*stateToSpecMap = malloc(hamilSize * sizeof **stateToSpecMap);
	
	for (long long int i=0LL; i < hamilSize; i++) {
		double energy = hamil[i];
		
		// find where this energy lives in our arrays
		int located=0;
		for (long long int j=0LL; j < *spectrumSize; j++) {
			if (energiesContainer[j] == energy) {
				populationContainer[j] += 1;
				(*stateToSpecMap)[i] = j;
				located=1;
				break;
			}
		}
		
		// if it's a new energy, add to the arrays
		if (!located) {
			energiesContainer[*spectrumSize] = energy;
			populationContainer[*spectrumSize] = 1;
			(*stateToSpecMap)[i] = *spectrumSize; 
			*spectrumSize += 1;
		}
	}
	
	// copy over arrays to tighter, dynamic ones
	*spectrum = malloc(*spectrumSize * sizeof **spectrum);
	*degeneracy = malloc(*spectrumSize * sizeof **degeneracy);
	for (long long int i=0LL; i < *spectrumSize; i++) {
		(*spectrum)[i] = energiesContainer[i];
		(*degeneracy)[i] = populationContainer[i];
	}
	
}


int main(int narg, char *varg[]) {
	
	
	
	/*
	 * GET CMD ARGS 
	 */
	 	 
	if (narg != 10) {
		printf("ERROR! Call with arguments:\n");
		printf("num_bools\nnum_params\nrseed\nthreshold[0 to 1]\ntimestep[0 for auto]\nmax_iters\nwrap_params\n");
		printf("deriv_accuracy[1 to 4]\nmatrix_noise[0 to 1]\n");
		return 1;
	}
	
	int numBools = atoi(varg[1]);
	int numParams = atoi(varg[2]);
	srand(atoi(varg[3]));
	double threshold; sscanf(varg[4], "%lf", &threshold);
	double timeStep; sscanf(varg[5], "%lf", &timeStep);
	int maxIterations = atoi(varg[6]);
	int wrapParams = atoi(varg[7]);
	int derivAccuracy = atoi(varg[8]);
	double matrNoise; sscanf(varg[9], "%lf", &matrNoise);
	
	if (numBools < 4) {
		printf("ERROR! Minimum num_bools is 4\n");
		return 1;
	}
	if (numParams < 1) {
		printf("ERROR! Minimum num_params is 1\n");
		return 1;
	}
	if (threshold <= 0 || threshold >= 1) {
		printf("ERROR! Threshold must be in (0, 1)\n");
		return 1;
	}
	if (timeStep < 0) {
		printf("ERROR! time_step must be positive or 0 for auto");
		return 1;
	}
	
	
	
	
	/* debug
	 * 
	 */
	 
	

	
	
	
	/*
	 * PREPARE SIMULATION
	 */
	 
	// TEST: 10qubit pauli file
	Hamiltonian hamil = loadPauliHamilFromFile("hamtest.txt");
	//printHamil(hamil);
	
	// monkeypatch
	long long int solState = 0;
	if (timeStep == 0)
		timeStep = 0.01;
	
	// generate a random 3SAT problem
	/*
	int *equ, *sol;
	int numClauses;
	Hamiltonian hamil = getRandom3SATHamil(numBools, &equ, &sol, &numClauses);
	printHamil(hamil);
	
	// get state index of solution (for monitoring progress of QuEST state)
	int solState = 0;
	for (int i=0; i < numBools; i++)
		solState += sol[numBools - i - 1] << i;
	
	print3SATEquSol(equ, sol, numBools, numClauses);
	printf("\nsol ind:\n%d\n\n", solState);
	
	// choose a sufficiently small time-step
	if (timeStep == 0)
		timeStep = getStableTimeStep(hamil.diagHamil, pow(2, numBools));
	printf("Time step: %lf\n", timeStep);
	*/
	
	// prepare QuEST
	QuESTEnv env;
	initQuESTEnv(&env);
	MultiQubit qubits; 
	createMultiQubit(&qubits, numBools, env);
	
	// prepare the param evolver
	EvolverMemory mem = prepareEvolverMemory(qubits, numParams);
	
	// set initial param values
	double params[numParams];
	for (int i=0; i < numParams; i++)
		params[i] = (rand()/(double) RAND_MAX) * 2 * M_PI;
	
	
	// remove energy degeneracy
	/*
	for (long long int i=0LL; i < qubits.numAmps; i++)
		if (hamil[i] != 0.0)
			hamil[i] += 10*(rand()/(double) RAND_MAX);
	*/
	// remove some excited terms
	/*
	int c=0;
	for (long long int i=0LL; i < qubits.numAmps; i++) {
		if (hamil[i] == 1.0) {
			c += 1;
		
			if (c == 13) {
				hamil[i] = 0;
				break;
			}
		}
	}
	*/

	
	
	/*
	 * PREPARE DATA RECORDING
	 */
	
	// prepare records of param values
	double paramEvo[numParams][maxIterations];
	for (int i=0; i < numParams; i++)
		for (int j=0; j < maxIterations; j++)
			paramEvo[i][j] = -666;
	
	// prepare records of sol prob and expected energy
	double solProbEvo[maxIterations];
	double expectedEnergyEvo[maxIterations];
	for (int i=0; i < maxIterations; i++) {
		solProbEvo[i] = -1;
		expectedEnergyEvo[i] = -1;
	}
	
	// analyse spectrum (only valid for diagonal hamiltonians)
	/*
	double* spectrum;
	int* degeneracy;
	int* stateToSpecMap;
	long long int spectrumSize;
	long long int solStateSpecInd = 0;
	getSpectrum(hamil.diagHamil, qubits.numAmps, &spectrum, &degeneracy, &stateToSpecMap, &spectrumSize);
	
	printf("\nSpectrum size:\t%lld\n", spectrumSize);
	for (long long int i=0LL; i < spectrumSize; i++)
		printf("energy:\t%lf,\tdegeneracy:\t%d\n", spectrum[i], degeneracy[i]);
	for (long long int i=0LL; i < spectrumSize; i++)
		if (spectrum[i] == 0)
			solStateSpecInd = i;
	printf("spectrum[%lld] = 0\n\n", solStateSpecInd);
	
	// prepare records of spectrum evolution
	double specProbEvo[spectrumSize][maxIterations];
	for (long long int i=0LL; i < spectrumSize; i++)
		for (int j=0; j < maxIterations; j++)
			specProbEvo[i][j] = -666;
	*/
	
	
	/*
	 * PERFORM SIMULATION
	 */
	
	evolveOutcome outcome;
	int step=0;
	double prob = 0;
	double energy = 0;
	
	// set the state we'll feed into the ansatz
	initStateZero(&qubits);
	setAnsatzInitState(&mem, qubits);
	
	defaultAnsatzCircuit(&mem, qubits, params, numParams);
	
	// keep evolving until we converge or reach max iterations
	while (step < maxIterations && prob < threshold ) {
		
		// update params under parameterised evolution
		outcome = evolveParams(
			&mem, defaultAnsatzCircuit, approxParamsByTikhonov,
			qubits, params, hamil, timeStep, wrapParams, derivAccuracy, matrNoise);
			
		if (outcome == FAILED) {
			printf("Numerical inversion failed! Aborting entire sim!\n");
			return 1;
		}
		
		// update params under exact evolution
		//evolveWavefunction(qubits, hamil, timeStep);
		
		// monitor convergence
		prob = getProbEl(qubits, solState);
		energy = getExpectedEnergy(mem.hamilState, qubits, hamil);
		printf("t%d: \t prob(sol) = %f \t <E> = %f\n", step, prob, energy);
		
		// record param evo data
		solProbEvo[step] = prob;
		expectedEnergyEvo[step] = energy;
		for (int i=0; i < numParams; i++)
			paramEvo[i][step] = params[i];
			
		// record spectrum evo data
		/*
		for (int i=0; i < spectrumSize; i++)
			specProbEvo[i][step] = 0;
		for (long long int i=0LL; i < qubits.numAmps; i++)
			specProbEvo[stateToSpecMap[i]][step] += getProbEl(qubits, i);
		*/
			 
		
		// randomly wiggle params
		/*
		if (step > 0 && step % 50 == 0) {
			printf("Wiggling params!");
			
			for (int i=0; i < numParams; i++)
				params[i] += 0.01 * 2*M_PI*(rand() / (double) RAND_MAX);
		}
		*/
		
		step += 1;
	}
	
	
	
	/*
	 * SAVE RESULTS TO FILE
	 */
	
	// record <E>, prob evolutions, recovered iterations
	FILE* file = openAssocWrite(OUTPUT_FILE);
	writeIntToAssoc(file, "numIterations", step);
	writeIntToAssoc(file, "derivAccuracy", derivAccuracy);
	writeDoubleToAssoc(file, "matrNoise", matrNoise, 5);
	writeIntToAssoc(file, "wrapParams", wrapParams);
	writeDoubleToAssoc(file, "timeStep", timeStep, 10);
	writeIntToAssoc(file, "numBools", numBools);
	writeIntToAssoc(file, "numParams", numParams);
	writeDoubleToAssoc(file, "threshold", threshold, 20);
	writeDoubleArrToAssoc(file, "solProbEvo", solProbEvo, step, 10);
	writeDoubleArrToAssoc(file, "expectedEnergyEvo", expectedEnergyEvo, step, 10);
	
	// record evolution of every param
	char buf[1024];
	for (int p=0; p < numParams; p++) {
		sprintf(buf, "param%dEvo", p);
		writeDoubleArrToAssoc(file, buf, paramEvo[p], step, 10);
	}
	
	// record evolution of the spectrum
	/*
	writeIntToAssoc(file, "spectrumSize", spectrumSize);
	writeDoubleArrToAssoc(file, "spectrum", spectrum, spectrumSize, 1);
	writeIntArrToAssoc(file, "spectrumDegeneracy", degeneracy, spectrumSize);
	writeIntToAssoc(file, "solStateSpecInd", solStateSpecInd);
	for (int s=0; s < spectrumSize; s++) {
		sprintf(buf, "spec%dEvo", s);
		writeDoubleArrToAssoc(file, buf, specProbEvo[s], step, 10);
	}
	*/
	
	closeAssocWrite(file);
	
	
	
	/*
	 * TIDY UP
	 */
	
	// cleanup
	freeEvolverMemory(&mem);
	freeHamil(hamil);
	/*
	free(equ);
	free(sol);
	free(spectrum);
	free(degeneracy);
	free(stateToSpecMap);
	*/
	
	// unload QuEST
	destroyMultiQubit(qubits, env); 
	closeQuESTEnv(env);
	return 0;
}