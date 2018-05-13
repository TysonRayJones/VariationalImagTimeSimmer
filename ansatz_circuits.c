/** @file 
 * A collect of ansatz circuits for variational chemistry or optimisation solving
 */
 
#include "ansatz_circuits.h"

#include <QuEST.h>
#include <complex.h>

#include "param_evolver.h"


int getIdealNumParamsInAnsatz(void (*ansatzCircuit)(EvolverMemory *mem, MultiQubit, double*, int), int numQubits) {
	
	if (ansatzCircuit == hardwareEfficientChemistryAnsatzCircuit)
		return 3*numQubits;
		
	// but you can't double this: the first numQubits are overhead!
	if (ansatzCircuit == lowDepthChemistryAnsatzCircuit)
		return numQubits + 5*(numQubits/2 + (numQubits-1)/2);
		
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