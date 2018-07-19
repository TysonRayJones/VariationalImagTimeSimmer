> The below instructions are out of date, since the interface has become more flexible.
> See [hamiltonian_builder.h](hamiltonian_builder.h) for creating a Hamiltonian, and see the 
> ```C
> /*
>  * PERFORM SIMULATION
>  */
> ```
> section of [variational_3SAT_solver.c](variational_3SAT_solver.c) to see an example of how to use [param_evolver.h](param_evolver.h).
> Feel free to contact me (via the email address advertised at [qtechtheory.org/tyson.jones](https://qtechtheory.org/person/tyson.jones/)) with questions and guidance.

Variational Imag-Time 3SAT Solver
========

# Setup

You'll need to install the GNU Scientific library; instructions [here](https://gist.github.com/TysonRayJones/af7bedcdb8dc59868c7966232b4da903).

The executable `Variational3SATSolver` has main in `variational_3SAT_solver.c` and includes
- QuEST
- `sat_generator` which generates the 3SAT problems **:^)**
- `hamiltonian_builder` which builds the Hamiltonians **:^)**
- `param_evolver` which evolves the parameters **:^)** with the variational imagianry time method
- `mmaformatter` which formats the output for importing into Mathematica

The QuEST `makefile` must be adapted for GSL. The changes are outlined [here](https://gist.github.com/TysonRayJones/af7bedcdb8dc59868c7966232b4da903) and my makefile is included in this repo.

# Running

THESE INSTRUCTIONS ARE OUT OF DATE: instead, run
```bash
./Variational3SATSolver help
```
OLD INSTRUCTIONS:

`Variational3SATSolver` accepts 9 command-line arguments:
- `num_bools`: 3SAT problem size (a random instance will be generated)
 - `num_params`: the number of parameterised gates in the ansatz circuit
 -  `rseed`: random seed for 3SAT generation 
 - `threshold`: the probability of the solution state which when achieved, parameter evolution ends
 - `timestep`: multiplies the change the parameters experience at each imag-time evolution step. Passing `0` will use the smallest stable time-step, ~**max(energyEigVal)**
 - `max_iters`: maximum iterations to perform before halting (if not already converged to `threshold`)
 -  `wrap_params`: whether (0 or 1) parameters will be wrapped-around to remain in **[0, 2 PI)** after each evolution step.
 - `deriv_accuracy`: order of finite difference deriv approx to use (1 to 4) in constructing A and C
 - `matrix_noise`: percent (0 to 1) noise to vary A and C elements by before solving for parameter evolution


# Customising

The parameters are updated by the method
```C
evolveOutcome evolveParams(
	evolverMemory *mem, 
	void (*ansatzCircuit)(MultiQubit, double*, int), 
	int (*inversionMethod)(evolverMemory*),
	MultiQubit qubits, 
	double* params, 
	Hamiltonian hamiltonian, 
	double timeStepSize, 
	int wrapParams
);
```
where `ansatzCircuit` and `inversionMethod` are functions you can define and pass, and `mem` is a data structure you must first create with
```C
prepareEvolverMemory(MultiQubit qubits, int numParams);
```
-------
 `ansatzCircuit` must be a function which accepts a `MultiQubit` instance, the parameter list and its length, and applies an ansatz circuit for those parameter values to `MultiQubit`. It's important this ansatz circuit sets the initial state of the qubits. For example
```C
void myAnsatz(MultiQubit qubits, double* params, int numParams) {
	
	initStateZero(&qubits);

	for (int p=0; p < numParams; p++)
		rotateX(qubits, p % qubits.numQubits, params[p]);
}
```
could be passed to ```evolveOutcome(..., myAnsatz, ...) ```
In lieu of your own code, `defaultAnsatzCircuit` can be passed.

-------

`inversionMethod` is called to numerically solve for the change in the parameters given the matrix equation.
The function must modify `mem->paramChange` (which is a GSL vector) using the `mem->matrA` and `mem->vecC` matrix/vector, to generate a solution to `matrA paramChange = vecC`, and should return `1` if it also numerically fails.

In lieu of your own code, you can pass `approxByLUDecomp`, `approxParamsByLeastSquares`, `approxParamsByRemovingVar`, `approxParamsByTSVD` or `approxParamsByTikhonov`; see the code-doc in `param_evolver.h` for details. I recommend the latter, which ensures `||paramChange||` is small.

Underdetermined equations of `paramChange` can be realised by changing `initStateZero` to `initStatePlus` in `defaultAnsatzCircuit`.

---------
`hamiltonian` must be an instance of `Hamiltonian`, declared in `hamiltonian_builder.c`, which you can construct via...

```C
Hamiltonian getRandom3SATHamil(int numBools, int **equ, int **sol, int *numClauses);
```

to generate a random 3SAT equation and accompanying diagonal Hamiltonian, or...

```C
Hamiltonian load3SATHamilFromFile(char *filename, int **equ, int **sol, int *numBools, int *numClauses);
```

to load a diagonal Hamiltonian from a 3SAT saved in a file, or...

```C
Hamiltonian loadPauliHamilFromFile(char *filename);
```

to load a Hamiltonian specified as a sum of pauli matrix products, with format
```CSV
coeff gateOnQubit1 gateOnQubit2 ... gateOnQubitN
coeff gateOnQubit1 gateOnQubit2 ... gateOnQubitN
```
where `coeff` is a double, `gateOnQubitn` is one of `0,1,2,3`, signifying an `I,X,Y,Z` gate on that qubit.

# Output

`Variational3SATSolver` (through `variational_3SAT_solver.c`) outputs a Mathematica association (using [this](https://gist.github.com/TysonRayJones/f25cb847aadd70aef3f5e0f8fae04947)) to file, which contains
- expected energy over iterations
- probability of the solution state over iterations
- evolution of every parameter
- which iterations underwent recovery numerical methods

and is simply read via 
```Mathematica 
SetDirectory @ NotebookDirectory[];
data = Get["wickSATdata.txt"];
```

An example `analysis.nb` notebook is included which plots the above data:
![mm](https://qtechtheory.org/wp-content/uploads/2018/04/example.png)
