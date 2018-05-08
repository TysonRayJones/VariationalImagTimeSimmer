/** @file 
 * Solves 3SAT problems by Wick rotation
 */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <QuEST.h>					// for simulating ansatz

#include "hamiltonian_builder.h"	// for building and applying Hamiltonians
#include "param_evolver.h"			// for variational imag-time simulation
#include "true_evolver.h"			// for verifying variational sim
#include "mmaformatter.h"			// for outputting results to mathematica


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
	Hamiltonian hamilObj, long long int hamilSize, 
	double** spectrum, int** degeneracy, int** stateToSpecMap, long long int* spectrumSize) {
		
	// monkey-patch: abort for non-diagonal Hamiltonians
	if (hamilObj.type != DIAGONAL) {
		printf("Spectrum analysis skipped for non-diagonal Hamiltonian.\n");
		*spectrumSize=0;
		return;
	}
		
		
	double *hamil = hamilObj.diagHamil;
	
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
	
	// list of args
	int numBools;
	int numParams;
	int randSeed;
	double timeStep;
	int maxIterations;
	int wrapParams;
	int derivAccuracy;
	double matrNoise;
	int simRepetitions;
	int progressPrintFrequency;
	
	// interactive arg input
	if (narg == 2 && !strcmp(varg[1], "help")) {
		
		printf("\n");
		printf("Enter the number of qubits (or booleans in the 3SAT equation)\n(must be >= 4)\n");
		printf("num_bools: ");
		scanf("%d", &numBools);
		printf("\n");
		
		printf("Enter the number of parameters total in the ansatz circuit\n");
		printf("num_params: ");
		scanf("%d", &numParams);
		printf("\n");
		
		printf("Enter the random seed for generating 3SATs, initial param values, matrix noise " 
			   "and GSL numerics\n");
		printf("rseed: ");
		scanf("%d", &randSeed);
		printf("\n");
		
		printf("Enter the time-step (float), or 0 for auto\n");
		printf("timestep: ");
		scanf("%lf", &timeStep);
		printf("\n");
		
		printf("Enter the maximum number of variational iterations in each simulation");
		printf("max_iters: ");
		scanf("%d", &maxIterations);
		printf("\n");
		
		printf("Enter whether to wrap-around params to keep them in [0, 2PI)\n"
			   "(1 for yes, 0 for no)\n");
		printf("wrap_params: ");
		scanf("%d", &wrapParams);
		printf("\n");
		
		printf("Enter the accuracy of the derivative estimates, as a finite-difference order\n"
			   "(1 for fastest, 4 for most accurate)\n");
		printf("deriv_accuracy: ");
		scanf("%d", &derivAccuracy);
		printf("\n");
		
		printf("Enter the fraction of random noise in the A and C matrices each iteration\n"
			   "(0 for no noise, e.g. 0.3 for maximum +- 30%% random fluctuation in each element)\n");
		printf("matrix_noise: ");
		scanf("%lf", &matrNoise);
		printf("\n");
		
		printf("Enter the number of times to resimulate the given system with different "
			   "random initial parameters\n");
		printf("sim_reps: ");
		scanf("%d", &simRepetitions);
		printf("\n");
		
		printf("Enter how frequently convergence progress should be printed, which slows execution\n"
			   "(0 for never, 1 for every iteration, n for every n iterations)\n");
		printf("print_progress_every: ");
		scanf("%d", &progressPrintFrequency);
		printf("\n");
		
	// invalid number of args input
	} else if (narg != 11) {
		printf("\nERROR! Call with arguments:\n");
		printf(
			"num_bools\n"
			"num_params\n"
			"rseed\n"
			"timestep[0 for auto]\n"
			"max_iters\n"
			"wrap_params\n"
			"deriv_accuracy[1 to 4]\n"
			"matrix_noise[0 to 1]\n"
			"sim_reps\n"
			"print_progress_every\n"
			"\n"
			"Run './Variational3SATSolver help' to enter arguments interactively\n\n");
			
		return 1;
		
	// cmd arg input
	} else {
	
		numBools = atoi(varg[1]);
		numParams = atoi(varg[2]);
		randSeed = atoi(varg[3]); 
		sscanf(varg[4], "%lf", &timeStep);
		maxIterations = atoi(varg[5]);
		wrapParams = atoi(varg[6]);
		derivAccuracy = atoi(varg[7]);
		sscanf(varg[8], "%lf", &matrNoise);
		simRepetitions = atoi(varg[9]);
		progressPrintFrequency = atoi(varg[10]);
	
	}
	
	// invalid arg values
	if (numBools < 4) {
		printf("ERROR! Minimum num_bools is 4\n");
		return 1;
	}
	if (numParams < 1) {
		printf("ERROR! Minimum num_params is 1\n");
		return 1;
	}
	if (timeStep < 0) {
		printf("ERROR! time_step must be positive or 0 for auto\n");
		return 1;
	}
	if (simRepetitions < 1) {
		printf("ERROR! sim_reps must be greater than 0\n");
		return 1;
	}
	if (progressPrintFrequency < 0) {
		printf("ERROR! print_progress_every must be positive (1 for every iteration, 0 for never)\n");
	}
	
	
	
	/*
	 * PREPARE SIMULATION
	 */
	
	srand(randSeed);
	
	// testing chemistry Hamiltonian
	/*
	Hamiltonian hamil = loadPauliHamilFromFile("hamtest6qb.txt");
	printHamil(hamil);
	
	// monkeypatch
	long long int solState = 0;
	if (timeStep == 0)
		timeStep = 0.01;
	*/
	
	
	// generate a random 3SAT problem
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
	if (timeStep == 0 && hamil.type == DIAGONAL)
		timeStep = getStableTimeStep(hamil.diagHamil, pow(2, numBools));
	printf("Time step: %lf\n", timeStep);
	
	
	// prepare QuEST
	QuESTEnv env;
	initQuESTEnv(&env);
	MultiQubit qubits; 
	createMultiQubit(&qubits, numBools, env);
	
	// prepare the param evolver
	EvolverMemory mem = prepareEvolverMemory(qubits, numParams);
	
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
	 
	// pre-allocate initial param values (GSL messes with seeding)
	double params[numParams];
	double initParams[simRepetitions][numParams];
	for (int s=0; s < simRepetitions; s++)
		for (int i=0; i < numParams; i++)
			initParams[s][i] = (rand()/(double) RAND_MAX) * 2 * M_PI;
	
	// prepare records of param values
	double paramEvo[simRepetitions][numParams][maxIterations];
	for (int i=0; i < simRepetitions; i++)
		for (int j=0; j < numParams; j++)
			for (int k=0; k < maxIterations; k++)
				paramEvo[i][j][k] = -666;
	
	// prepare records expected energy...
	double expectedEnergyEvo[simRepetitions][maxIterations];
	for (int s=0; s < simRepetitions; s++)
		for (int i=0; i < maxIterations; i++)
			expectedEnergyEvo[s][i] = -1;
	
	// and solution prob (only relevant for DIAGONAL Hamiltonians)
	double solProbEvo[simRepetitions][maxIterations];
	for (int s=0; s < simRepetitions; s++)
		for (int i=0; i < maxIterations; i++)
			solProbEvo[s][i] = -1;
	
	// analyse spectrum (only valid for diagonal hamiltonians)
	double* spectrum;
	int* degeneracy;
	int* stateToSpecMap;
	long long int spectrumSize;
	long long int solStateSpecInd = 0;
	getSpectrum(
		hamil, qubits.numAmps, 
		&spectrum, &degeneracy, &stateToSpecMap, &spectrumSize
	);
	double specProbEvo[simRepetitions][spectrumSize][maxIterations];
		
	if (hamil.type == DIAGONAL) {
		printf("\nSpectrum size:\t%lld\n", spectrumSize);
		for (long long int i=0LL; i < spectrumSize; i++)
			printf("energy:\t%lf,\tdegeneracy:\t%d\n", spectrum[i], degeneracy[i]);
		for (long long int i=0LL; i < spectrumSize; i++)
			if (spectrum[i] == 0)
				solStateSpecInd = i;
		printf("spectrum[%lld] = 0\n\n", solStateSpecInd);
		
		// prepare records of spectrum evolution
		for (int r=0; r < simRepetitions; r++)
			for (long long int i=0LL; i < spectrumSize; i++)
				for (int j=0; j < maxIterations; j++)
					specProbEvo[r][i][j] = -666;
	}
	
	
	
	/*
	 * PERFORM SIMULATION
	 */
		
	// set the state we'll feed into the ansatz
	initStateZero(&qubits);
	setAnsatzInitState(&mem, qubits);
	
	evolveOutcome outcome;
	double prob, energy;
	
	// re-simulate many times
	for (int rep=0; rep < simRepetitions; rep++) {
	
		// set random initial param values
		for (int i=0; i < numParams; i++)
			params[i] = initParams[rep][i];
		
		// keep evolving until we converge or reach max iterations
		for (int step = 0; step < maxIterations; step++) {
			
			// update params under parameterised evolution
			//if (step % 2 == 0) {
				outcome = evolveParams(
					&mem, defaultAnsatzCircuit, approxParamsByTikhonov,
					qubits, params, hamil, timeStep, wrapParams, derivAccuracy, matrNoise);
			/*} else {
				evolveParamsByGradientDescent(
					&mem, defaultAnsatzCircuit, qubits, params, hamil, timeStep, wrapParams, derivAccuracy
				);
			}*/
			
			if (outcome == FAILED) {
				printf("Numerical inversion failed! Aborting entire sim!\n");
				return 1;
			}
			
			
			if (step % 100 == 0) {
				exciteStateInHamiltonian(&mem, qubits);
				printf("EXCITED THE CURRENT STATE!\n");
				//clearExcitedStates(&mem);
				
				for (int i=0; i < numParams; i++)
					params[i] = initParams[rep][i];
			}
			
			// update params under exact evolution
			//evolveWavefunction(qubits, hamil, timeStep);
			
			// monitor convergence
			prob = getProbEl(qubits, solState);
			energy = getExpectedEnergy(mem.hamilState, qubits, hamil);
			if (progressPrintFrequency != 0 && step % progressPrintFrequency == 0)
				printf("t%d: \t prob(sol) = %f \t <E> = %f\n", step, prob, energy);
			
			// record param evo data
			solProbEvo[rep][step] = prob;
			expectedEnergyEvo[rep][step] = energy;
			for (int i=0; i < numParams; i++)
				paramEvo[rep][i][step] = params[i];
				
			// record spectrum evo data
			if (hamil.type == DIAGONAL) {
				for (int i=0; i < spectrumSize; i++)
					specProbEvo[rep][i][step] = 0;
				for (long long int i=0LL; i < qubits.numAmps; i++)
					specProbEvo[rep][stateToSpecMap[i]][step] += getProbEl(qubits, i);
			}
			
			// randomly wiggle params
			/*
			if (step > 0 && step % 50 == 0) {
				printf("Wiggling params!");
				
				for (int i=0; i < numParams; i++)
					params[i] += 0.01 * 2*M_PI*(rand() / (double) RAND_MAX);
			}
			*/
		}
		
		// clear any previously-saved excitations
		clearExcitedStates(&mem);
	}
	
	
	
	/*
	 * SAVE RESULTS TO FILE
	 */
	
	// record results
	FILE* file = openAssocWrite(OUTPUT_FILE);
	
	// meta-data
	writeIntToAssoc(file, "simRepetitions", simRepetitions);
	writeIntToAssoc(file, "maxIterations", maxIterations);
	writeIntToAssoc(file, "derivAccuracy", derivAccuracy);
	writeDoubleToAssoc(file, "matrNoise", matrNoise, 5);
	writeIntToAssoc(file, "wrapParams", wrapParams);
	writeDoubleToAssoc(file, "timeStep", timeStep, 10);
	writeIntToAssoc(file, "numBools", numBools);
	writeIntToAssoc(file, "numParams", numParams);
	
	// solution, energy, param evolution
	writeNestedDoubleArrToAssoc(file, "solProbEvos", solProbEvo, 2, (int []) {simRepetitions, maxIterations}, maxIterations, 10);
	writeNestedDoubleArrToAssoc(file, "expectedEnergyEvos", expectedEnergyEvo, 2, (int []) {simRepetitions, maxIterations}, maxIterations, 10);
	writeNestedDoubleArrToAssoc(file, "paramEvos", paramEvo, 3, (int []) {simRepetitions, numParams, maxIterations}, maxIterations, 10);
	
	if (hamil.type == PAULI_TERMS) {
		writeStringToAssoc(file, "hamilType", "PAULI_TERMS");
	}
	if (hamil.type == DIAGONAL) {
		
		// 3SAT equ
		writeStringToAssoc(file, "hamilType", "DIAGONAL");
		writeIntArrToAssoc(file, "3SATEqu", equ, numClauses*3);
		writeIntArrToAssoc(file, "3SATSol", sol, numBools);
		
		// spectrum evolution
		writeIntToAssoc(file, "spectrumSize", spectrumSize);
		writeDoubleArrToAssoc(file, "spectrum", spectrum, spectrumSize, 1);
		writeIntArrToAssoc(file, "spectrumDegeneracy", degeneracy, spectrumSize);
		writeIntToAssoc(file, "solStateSpecInd", solStateSpecInd);
		writeNestedDoubleArrToAssoc(file, "specEvos", specProbEvo, 3, (int []) {simRepetitions, spectrumSize, maxIterations}, maxIterations, 10);
	}

	closeAssocWrite(file);
	
	
	
	/*
	 * TIDY UP
	 */
	 
	if (hamil.type == DIAGONAL) {
		free(equ);
		free(sol);
		free(spectrum);
		free(degeneracy);
		free(stateToSpecMap);
	}
	freeHamil(hamil);
	freeEvolverMemory(&mem);
	destroyMultiQubit(qubits, env); 
	closeQuESTEnv(env);
	
	return 0;
}