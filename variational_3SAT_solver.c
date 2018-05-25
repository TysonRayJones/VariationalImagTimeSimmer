/** @file 
 * Solves 3SAT and Chemistry problems by variational imaginary time propogation
 */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <QuEST.h>					// for simulating ansatz

#include "hamiltonian_builder.h"	// for building and applying Hamiltonians
#include "ansatz_circuits.h"		// for selecting ansatz circuits
#include "linear_solvers.h" 		// for selecting numerical solving method
#include "param_evolver.h"			// for variational imag-time simulation
#include "true_evolver.h"			// for verifying variational sim
#include "mmaformatter.h"			// for outputting results to mathematica


/** max energy distance beteween non-degen eigvals for compactifying */
#define DEGEN_EIGVAL_MAX_DIST 0.01

double getMaxProb(MultiQubit qubits) {
	
	double maxProb = 0;
	for (long long int i=0LL; i < qubits.numAmps; i++) {
		double prob = getProbEl(qubits, i);
		if (prob > maxProb)
			maxProb = prob;
	}
	return maxProb;
}



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


/*
 * eigvals must previously be sorted 
 */
 
void countSpectrumDegeneracy(double* eigvals, int numEigvals, double** uniqueEigVals, int** degeneracy, int* numUniqueEigVals) {
	
	double myEigVals[numEigvals];
	int myDegeneracy[numEigvals];
	
	int myInd = 0;
	myEigVals[myInd] = eigvals[0];
	myDegeneracy[myInd] = 1;
	
	// find unique energies
	for (int i=1; i < numEigvals; i++) {
		
		if (fabs(eigvals[i] - myEigVals[myInd]) < DEGEN_EIGVAL_MAX_DIST) {
			myDegeneracy[myInd]++;
		}
		else {
			myInd++;
			myEigVals[myInd] = eigvals[i];
			myDegeneracy[myInd] = 1;
		}
	}
	
	// move arrays to dynamic memory
	*numUniqueEigVals = myInd + 1;
	*uniqueEigVals = malloc(*numUniqueEigVals * sizeof **uniqueEigVals);
	*degeneracy = malloc(*numUniqueEigVals * sizeof **degeneracy);
	for (int i=0; i < *numUniqueEigVals; i++) {
		(*uniqueEigVals)[i] = myEigVals[i];
		(*degeneracy)[i] = myDegeneracy[i];
	}
}



int main(int narg, char *varg[]) {
	
	
	
	/*
	 * GET CMD ARGS 
	 */
	
	// list of args
	int numQubits;
	int numParams;
	int randSeed;
	double timeStep;
	int maxIterations;
	int wrapParams;
	int derivAccuracy;
	double matrNoise;
	int useGradDesc;
	char* ansatzCircString;
	int exciteWhenStuck;
	int simRepetitions;
	int progressPrintFrequency;
	char* outputFilename;
	
	// interactive arg input
	if (narg == 2 && !strcmp(varg[1], "help")) {
		
		printf("\n");
		printf("Enter the number of qubits (or booleans in the 3SAT equation)\n(must be >= 4)\n");
		printf("num_bools: ");
		scanf("%d", &numQubits);
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
		
		printf("Enter the maximum number of variational iterations in each simulation\n");
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
		
		printf("Enter whether to use imaginary-time evolution (0) or gradient descent (1)\n");
		printf("use_gd: ");
		scanf("%d", &useGradDesc);
		printf("\n");
		
		// TODO: get ansatz
		printf("Enter ansatz circuit\n");
		printf("ansatz: ");
		ansatzCircString = "lowDepth";
		printf("\nUSING DEFAULT lowDepth (until I fix this)\n\n");
		
		printf("Enter whether to excite out of stuck states (Suguru's method)\n");
		printf("excite_when_stuck: ");
		scanf("%d", &exciteWhenStuck);
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
		
		// TODO: get output filename
		printf("Enter the output filename\n");
		printf("output_fn: ");
		outputFilename = "VarData.txt";
		printf("\nUSING DEFAULT VarData.txt (until I fix this)\n\n");
		
		
	// invalid number of args input
	} else if (narg != 15) {
		printf("\nERROR! Call with arguments:\n");
		printf(
			"num_bools\n"
			"num_params[0 for auto]\n"
			"rseed\n"
			"timestep[0 for auto]\n"
			"max_iters\n"
			"wrap_params\n"
			"deriv_accuracy[1 to 4]\n"
			"matrix_noise[0 to 1]\n"
			"use_gd[0 or 1]\n"
			"ansatz\n"
			"excite_when_stuck[0 or 1]\n"
			"sim_reps\n"
			"print_progress_every\n"
			"output_fn\n"
			"\n"
			"Run './Variational3SATSolver help' to enter arguments interactively\n\n");
			
		return 1;
		
	// cmd arg input
	} else {
	
		numQubits = atoi(varg[1]);
		numParams = atoi(varg[2]); ////////////////////
		randSeed = atoi(varg[3]); 
		sscanf(varg[4], "%lf", &timeStep);
		maxIterations = atoi(varg[5]);
		wrapParams = atoi(varg[6]);
		derivAccuracy = atoi(varg[7]);
		sscanf(varg[8], "%lf", &matrNoise);
		useGradDesc = atoi(varg[9]);
		ansatzCircString = varg[10];
		exciteWhenStuck = atoi(varg[11]);
		simRepetitions = atoi(varg[12]);
		progressPrintFrequency = atoi(varg[13]);
		outputFilename = varg[14];
	}
	
	// invalid arg values
	if (numQubits < 4) {
		printf("ERROR! Minimum num_bools is 4\n");
		return 1;
	}
	if (numParams < 0) {
		printf("ERROR! Minimum num_params is 1, or 0 for auto\n");
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
		return 1;
	}
	
	AnsatzCircuit ansatz = getAnsatzFromString(ansatzCircString);
	if (ansatz == NULL) {
		printf("ERROR! Unrecognised ansatz circuit '%s'\n", ansatzCircString);
		return 1;
	}

	// apply (some of the) auto-parameters
	if (numParams == 0)
		numParams = getIdealNumParamsInAnsatz(ansatz, numQubits);
		
	
	// TODO: use ansatz circuit
	
	
	// confirm args
	printf("\n");
	printf("numQubits: %d\n", numQubits);
	printf("numParams: %d\n", numParams);
	printf("randSeed: %d\n", randSeed);
	if (timeStep == 0)
		printf("timeStep: TBC\n");
	else
		printf("timeStep: %lf\n", timeStep);
	printf("iterations: %d\n", maxIterations);
	printf("wrapParams: %d\n", wrapParams);
	printf("derivAccuracy: %d\n", derivAccuracy);
	printf("matrNoise: %lf\n", matrNoise);
	printf("useGradDesc: %d\n", useGradDesc);
	printf("ansatzCircuit: %s\n", ansatzCircString);
	printf("exciteWhenStuck: %d\n", exciteWhenStuck);
	printf("simRepettions: %d\n", simRepetitions);
	printf("printEvery: %d\n", progressPrintFrequency);
	printf("outputFilename: %s\n", outputFilename);
	printf("\n");
	
	
	/*
	 * PREPARE SIMULATION
	 */
	
	srand(randSeed);
	
	// testing chemistry Hamiltonian
	Hamiltonian hamil = loadPauliHamilFromFile("chemHamil6qb.txt");
	if (numQubits == 10) {
		
		printf("USING 10 QUBIT HAMIL\n");
		freeHamil(hamil);
		hamil = loadPauliHamilFromFile("chemHamil10qb.txt");
	}
	//printHamil(hamil);

	
	// monkeypatch
	long long int solState = 0;
	if (timeStep == 0)
		timeStep = 0.01;
	
	
	
	// generate a random 3SAT problem
	int *equ = NULL;
	int *sol = NULL;
	int numClauses = -1;
	/*
	Hamiltonian hamil = getRandom3SATHamil(numQubits, &equ, &sol, &numClauses);
	printHamil(hamil);
	
	// get state index of solution (for monitoring progress of QuEST state)
	long long int solState = 0;
	for (int i=0; i < numQubits; i++)
		solState += sol[numQubits - i - 1] << i;
	
	print3SATEquSol(equ, sol, numQubits, numClauses);
	printf("\nsol ind:\n%lld\n\n", solState);
	
	// choose a sufficiently small time-step
	if (timeStep == 0 && hamil.type == DIAGONAL)
		timeStep = getStableTimeStep(hamil.diagHamil, pow(2, numQubits));
	printf("Time step: %lf\n", timeStep);
	*/
	
	
	
	// prepare QuEST
	QuESTEnv env;
	initQuESTEnv(&env);
	MultiQubit qubits; 
	createMultiQubit(&qubits, numQubits, env);
	
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
	
	
	
	
	// HAMIL DIAG TESTING!!!
	
	double* eigvals = NULL;
	double complex** eigvecs = NULL;
	
	int numUniqEigvals = -1;
	double* uniqEigvals = NULL;
	int* eigvalDegen = NULL;

	if (hamil.type == PAULI_TERMS) {
		
		printf("Finding spectrum...\n");
		getPauliHamilEigvals(hamil, qubits, hamil.numAmps, &eigvals, &eigvecs);
		printf("Compacting spectrum...\n");
		countSpectrumDegeneracy(eigvals, hamil.numAmps, &uniqEigvals, &eigvalDegen, &numUniqEigvals);
	}
	
	
	/*
	 * PREPARE DATA RECORDING
	 */
	 
	 
	//double samsInitParams[42] ={3.243044, 4.81527666, 0.62511592, 0.41921514, 4.24673659, 3.43911744, 0.0877967, 6.06145306, 1.08071278, 6.10536731, 4.16570543, 2.37626469, 5.71440099, 2.77030973, 0.46173575, 3.36354362, 4.70542023, 6.11022712, 2.45688627, -0.11424261, 3.69617933, 1.098452, 4.71331006, 0.99540201, 1.01557048, -0.05416214, 1.13185388, 4.36649297, 0.44485088, 1.89563555, 0.847721051, 2.44248901, 2.76156372, 3.24083925, 0.78477016, 5.50861965, 0.21425296, 2.59175576, 2.44989661, 4.19120921, 1.88714856, 3.23731343};
	//double samsFinalParams[42] = {-10.85422203, 87.96330825, 1.16092286E-05, -0.05925503, -3.08235636, 0.00205020046, -18.84501712, 45.39287869, 3.54246012, 51.8369899, 35.95813837, -32.98477717, -1011.59145386, -3.39717838, -18.78717354, 3.22934895, 3.14376364, 8.69201468, 15.70342282, 74.22959909, 11.0544325, -48.36063456, -221.26090968, 45.59825606, 13.99581773, -0.000161766399, 3.14141416, 0.000886289673, 15.71206725, 3.14159531, 54.64032903, -7.85421824, 4.77108247, -1.7355999, -199.49040412, -12.10106662, 3.14141678, 6.2835601, 6.27902231, 2.39810243, 3.14181229, 0.00264922488};
	
	
	// pre-allocate initial param values (GSL messes with seeding)
	double params[numParams];
	double initParams[simRepetitions][numParams];
	for (int s=0; s < simRepetitions; s++)
		for (int i=0; i < numParams; i++)
			initParams[s][i] = M_PI; //(rand()/(double) RAND_MAX) * 2 * M_PI; //samsInitParams[i]; //0.005;    ///////////// CHANGE THIS BACCCCCCKKKKK
	
	// prepare records of param values
	double*** paramEvo = malloc(simRepetitions * sizeof *paramEvo);
	for (int i=0; i < simRepetitions; i++) {
		paramEvo[i] = malloc(numParams * sizeof **paramEvo);
		for (int j=0; j < numParams; j++) {
			paramEvo[i][j] = malloc(maxIterations * sizeof ***paramEvo);
			for (int k=0; k < maxIterations; k++)
				paramEvo[i][j][k] = -666;
		}
	}
	
	/*
	double paramEvo[simRepetitions][numParams][maxIterations];
	for (int i=0; i < simRepetitions; i++)
		for (int j=0; j < numParams; j++)
			for (int k=0; k < maxIterations; k++)
				paramEvo[i][j][k] = -666;
	*/
				
	
	
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
	// sigmaX(qubits, 0);    ///////////////////////////////////////////// NOTTED FOR HARTREE FOCK
	setAnsatzInitState(&mem, qubits);
	
	
	evolveOutcome outcome;
	double prob, energy;
	
	// re-simulate many times
	for (int rep=0; rep < simRepetitions; rep++) {
		
		int numExcitations = 0;
		
		// whethert his simulation has converged (not necessarily ground)
		//int hasBecomeStuck = 0;
	
		// set random initial param values
		for (int i=0; i < numParams; i++)
			params[i] = initParams[rep][i];
		
		// remove any static states from the Hamiltonian excitations
		clearExcitedStates(&mem);
		energy = 1E5;
		
		// keep evolving until we converge or reach max iterations
		for (int step = 0; step < maxIterations; step++) {
			
			// update params under parameterised evolution
			
			if (useGradDesc)
				outcome = evolveParamsByGradientDescent(
					&mem, ansatz, qubits, params, hamil, timeStep, wrapParams, derivAccuracy, matrNoise);
			else
				outcome = evolveParamsByImaginaryTime(
					&mem, ansatz, approxParamsByTikhonov,
					//approxParamsByLUDecomp,
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
			if (progressPrintFrequency != 0 && step % progressPrintFrequency == 0) {
				if (hamil.type == DIAGONAL)
					printf("t%d: \t prob(sol) = %f \t <E> = %f\n", step, prob, energy);
				else
					printf("t%d: \t<E> = %f\n", step, energy);
			}
			
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
			
			// mark if params have halted
			/*
			if ( //haltIterations[rep] == -1 && 
				isStuck((double*) paramEvo, rep, numParams, maxIterations, step) &&
				energy > 0.9) {
					
				haltIterations[rep] = step;
				
				exciteStateInHamiltonian(&mem, qubits);

				for (int i=0; i < numParams; i++)
					params[i] = initParams[rep][i];
			}
			*/
			
			if (exciteWhenStuck && isStuck(paramEvo, rep, numParams, step)) {
				//haltIterations[rep] = step;
				//hasBecomeStuck = 1;
				exciteStateInHamiltonian(&mem, qubits);
				
				printf("Exciting state!\n");
				printf("Resetting params...\n");
				
				for (int n=0; n < numParams; n++)
					params[n] = initParams[rep][n];
				
				numExcitations++;
			}
		}
	}
	
	
	
	/*
	 * SAVE RESULTS TO FILE
	 */
	
	// record results
	FILE* file = openAssocWrite(outputFilename);
	
	// meta-data
	writeIntToAssoc(file, "simRepetitions", simRepetitions);
	writeIntToAssoc(file, "maxIterations", maxIterations);
	writeIntToAssoc(file, "derivAccuracy", derivAccuracy);
	writeDoubleToAssoc(file, "matrNoise", matrNoise, 5);
	writeIntToAssoc(file, "wrapParams", wrapParams);
	writeDoubleToAssoc(file, "timeStep", timeStep, 10);
	writeIntToAssoc(file, "numQubits", numQubits);
	writeIntToAssoc(file, "numParams", numParams);
	writeIntToAssoc(file, "useGradDesc", useGradDesc);
			
	// solution, energy, param evolution, iteration when evolutio halted
	writeNestedDoubleArrToAssoc(file, "solProbEvos", solProbEvo, 2, (int []) {simRepetitions, maxIterations}, maxIterations, 10);
	writeNestedDoubleArrToAssoc(file, "expectedEnergyEvos", expectedEnergyEvo, 2, (int []) {simRepetitions, maxIterations}, maxIterations, 10);
	writeNestedDoubleArrToAssoc(file, "initParams", initParams, 2, (int []) {simRepetitions, numParams}, numParams, 10);
	//writeNestedDoubleArrToAssoc(file, "paramEvos", paramEvo, 3, (int []) {simRepetitions, numParams, maxIterations}, maxIterations, 10);
	writeNestedDoubleListToAssoc(file, "paramEvos", paramEvo, 3,  (int []) {simRepetitions, numParams, maxIterations}, 10);
	
	if (hamil.type == PAULI_TERMS) {
		writeStringToAssoc(file, "hamilType", "PAULI_TERMS");
		
		// write spectrum
		writeIntToAssoc(file, "spectrumSize", hamil.numAmps);
		writeDoubleArrToAssoc(file, "spectrum", eigvals, hamil.numAmps, 10);
		
		writeIntToAssoc(file, "uniqueSpectrumSize", numUniqEigvals);
		writeDoubleArrToAssoc(file, "uniqueSpectrum", uniqEigvals, numUniqEigvals, 10);
		writeIntArrToAssoc(file, "uniqueSpectrumDegeneracy", eigvalDegen, numUniqEigvals);
	
	}
	if (hamil.type == DIAGONAL) {
		
		// 3SAT equ
		writeStringToAssoc(file, "hamilType", "DIAGONAL");
		writeIntArrToAssoc(file, "3SATEqu", equ, numClauses*3);
		writeIntArrToAssoc(file, "3SATSol", sol, numQubits);
		
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
	for (int i=0; i < simRepetitions; i++) {
		for (int j=0; j < numParams; j++)
			free(paramEvo[i][j]);
		free(paramEvo[i]);
	}
	free(paramEvo);
	
	if (hamil.type == PAULI_TERMS) {
		freePauliHamilEigvals(eigvals, eigvecs, hamil.numAmps);
	}
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