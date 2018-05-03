
#include "true_evolver.h"

#include <math.h>
#include <complex.h>

#include "QuEST/qubits.h"

void evolveWavefunction(MultiQubit qubits, double* hamil, double timeStep) {
	
	// DE: dy/ds = - H y
	
	// get fourth-order runge-kutta imag-time-Schrodinger-equ terms
	double complex k1[qubits.numAmps];
	for (long long int i=0LL; i < qubits.numAmps; i++) 
		k1[i] = - timeStep * hamil[i] * (getRealAmpEl(qubits, i) + I*getImagAmpEl(qubits, i));
	
	double complex k2[qubits.numAmps];
	for (long long int i=0LL; i < qubits.numAmps; i++) 
		k2[i] = - timeStep * hamil[i] * (getRealAmpEl(qubits, i) + I*getImagAmpEl(qubits, i) + k1[i]/2.0);
		
	double complex k3[qubits.numAmps];
	for (long long int i=0LL; i < qubits.numAmps; i++) 
		k3[i] = - timeStep * hamil[i] * (getRealAmpEl(qubits, i) + I*getImagAmpEl(qubits, i) + k2[i]/2.0);
	
	double complex k4[qubits.numAmps];
	for (long long int i=0LL; i < qubits.numAmps; i++) 
		k4[i] = - timeStep * hamil[i] * (getRealAmpEl(qubits, i) + I*getImagAmpEl(qubits, i) + k3[i]);
	
	// update wavefunction
	double norm = 0;
	for (long long int i=0LL; i < qubits.numAmps; i++) {
		double complex dPsi = (k1[i] + k4[i])/6.0 + (k2[i] + k3[i])/3.0;
		qubits.stateVec.real[i] += creal(dPsi);
		qubits.stateVec.imag[i] += cimag(dPsi);
		
		norm += pow(qubits.stateVec.real[i],2) + pow(qubits.stateVec.imag[i], 2);
	}
	
	// normalise wavefunction
	norm = sqrt(norm);
	for (long long int i=0LL; i < qubits.numAmps; i++) {
		qubits.stateVec.real[i] /= norm;
		qubits.stateVec.imag[i] /= norm;
	}
}