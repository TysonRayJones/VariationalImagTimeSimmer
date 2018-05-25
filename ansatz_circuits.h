#ifndef ANSATZ_CIRCUITS_H_
#define ANSATZ_CIRCUITS_H_

#include <QuEST.h>
#include <complex.h>

#include "param_evolver.h"


typedef void (*AnsatzCircuit)(EvolverMemory *mem, MultiQubit, double*, int);

AnsatzCircuit getAnsatzFromString(char* string);

/**
 * Call after preparing qubits in the desired state to feed into
 * the ansatz at every evolution.
 * Can be recalled at any time during evolution to change the ansatz
 */
void setAnsatzInitState(EvolverMemory *mem, MultiQubit qubits);

void hardwareEfficientChemistryAnsatzCircuit(EvolverMemory *mem, MultiQubit qubits, double* params, int numParams);
void lowDepthChemistryAnsatzCircuit(EvolverMemory *mem, MultiQubit qubits, double* params, int numParams);

int getIdealNumParamsInAnsatz(AnsatzCircuit ansatz, int numQubits);


#endif // ANSATZ_CIRCUITS_H_