/** @file 
 * A collect of ansatz circuits for variational chemistry or optimisation solving
 */
 
#include "ansatz_circuits.h"

#include <QuEST.h>
#include <complex.h>
#include <string.h>

#include "param_evolver.h"


#define LOW_DEPTH_ANSATZ_STR "lowDepth"
#define HARDWARE_EFFICIENT_ANSATZ_STR "hardEff"

#define DEBUG_ANSATZ_PRINT 0



void debugPrint(char* gate, double param, int qb) {
	
	//param =  fmod(param,2*M_PI);
	//if (param < 0)
	//	param += 2*M_PI;
	
	if (DEBUG_ANSATZ_PRINT) 
		printf("%s(%lf) | qb%d\n", gate, param, qb);
		
}


AnsatzCircuit getAnsatzFromString(char* string) {
	
	AnsatzCircuit ansatz = NULL;
	if (strcmp(string, LOW_DEPTH_ANSATZ_STR) == 0)
		ansatz = lowDepthChemistryAnsatzCircuit;
		
	else if (strcmp(string, HARDWARE_EFFICIENT_ANSATZ_STR) == 0)
		ansatz = hardwareEfficientChemistryAnsatzCircuit;
		
	return ansatz;
}

int getIdealNumParamsInAnsatz(AnsatzCircuit ansatz, int numQubits) {
	
	// 2 blocks (first block has 3 columns, subsequent blocks have 4)
	if (ansatz == hardwareEfficientChemistryAnsatzCircuit)
		return (3 + 1*4)*numQubits;
		
	// 2 blocks (once-off initial column) 
	if (ansatz == lowDepthChemistryAnsatzCircuit)
		return numQubits + 2*5*(numQubits/2 + (numQubits-1)/2);
		
	printf("ERROR: unrecognised ansatz circuit!\n");
	return -1;
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
 * Applies a hardware efficient parameterised ansatz circuit to the zero state, modifying
 * the wavefunction in qubits according to the values in params. This can be passed
 * to evolveParams in lieu of a custom ansatz circuit
 */
void hardwareEfficientChemistryAnsatzCircuit(EvolverMemory *mem, MultiQubit qubits, double* params, int numParams) {
	
	initState(mem, qubits);
	
	int block = 0;
	int paramInd = 0;
	
	// loop the following circuit until all params are featured
	while (paramInd < numParams) {
		
		// skip Rz column in the first block
		if (paramInd > 0) 
			for (int qb=0; qb < qubits.numQubits && paramInd < numParams; qb++) {
				debugPrint("RZ", params[paramInd], qb);
				rotateZ(qubits, qb, params[paramInd++]);
			}
		
		// Rx
		for (int qb=0; qb < qubits.numQubits && paramInd < numParams; qb++) {
			debugPrint("RX", params[paramInd], qb);
			rotateX(qubits, qb, params[paramInd++]);
		}
		
		// Rz
		for (int qb=0; qb < qubits.numQubits && paramInd < numParams; qb++) {
			debugPrint("RZ", params[paramInd], qb);
			rotateZ(qubits, qb, params[paramInd++]);
		}
		
		// C(Ry)
		for (int qb=0; qb < qubits.numQubits && paramInd < numParams; qb++) {
			debugPrint("CY", params[paramInd], qb);
			controlledRotateY(qubits, qb, (qb+block+1)%qubits.numQubits, params[paramInd++]);
		}
		
		block++;
	}
}

/**
 * Applies a low-depth parameterised ansatz circuit to the zero state, modifying
 * the wavefunction in qubits according to the values in params. This can be passed
 * to evolveParams in lieu of a custom ansatz circuit
 */
void lowDepthChemistryAnsatzSubcircuit(MultiQubit qubits, double param, int qb0, int qb1, int type0, int type1) {
	
	if (type0 == 0)
		hadamard(qubits, qb0);
	if (type0 == 1)
		rotateX(qubits, qb0, M_PI_2);
		
	if (type1 == 0)
		hadamard(qubits, qb1);
	if (type1 == 1)
		rotateX(qubits, qb1, M_PI_2);
	
	controlledNot(qubits, qb0, qb1);
	rotateZ(qubits, qb1, param);
	controlledNot(qubits, qb0, qb1);
	
	if (type0 == 0)
		hadamard(qubits, qb0);
	if (type0 == 1)
		rotateX(qubits, qb0, - M_PI_2);
		
	if (type1 == 0)
		hadamard(qubits, qb1);
	if (type1 == 1)
		rotateX(qubits, qb1, - M_PI_2);
}
void lowDepthChemistryAnsatzCircuit(EvolverMemory *mem, MultiQubit qubits, double* params, int numParams) {
	
	initState(mem, qubits);
	
	const int numTypes = 5;
	int types[5][2] = {{1,0}, {0,1}, {2,2}, {1,1}, {0,0}};
	
	int paramInd = 0;
	// initial Rx
	for (int qb=0; qb < qubits.numQubits && paramInd < numParams; qb++)
		rotateX(qubits, qb, params[paramInd++]);
	
	// loop the following circuit until all params are featured
	while (paramInd < numParams) {
		
		// outer sub-circuits
		for (int qb=0; qb < qubits.numQubits -1 && paramInd < numParams; qb += 2)
			for (int t=0; t < numTypes && paramInd < numParams; t++)
				lowDepthChemistryAnsatzSubcircuit(qubits, params[paramInd++], qb, qb+1, types[t][0], types[t][1]);
		
		// inner sub-circuits
		for (int qb=1; qb < qubits.numQubits -1 && paramInd < numParams; qb += 2) 
			for (int t=0; t < numTypes && paramInd < numParams; t++)
				lowDepthChemistryAnsatzSubcircuit(qubits, params[paramInd++], qb, qb+1, types[t][0], types[t][1]);
	}
}