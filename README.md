# Cerebral hemodynamic model.

A lumped parameter model and solver developed as part of my masters thesis. 

The model includes a description of detailed systemic and cerebral hemodynamics, as well as physiological control mechanisms. It is capable of simulating realistic hemodynamics in healthy conditions as well as atrial fibrillation.
June 3, 2021

## Usage.

First create a clone of the repository.
Then run the following commands.
```sh
cd cerebral-od-model 
make cbf
```

The program requires a directory named "inputs" at in the same directory as the binary. Inputs should contain files pnkNoise{run_index}.dat and expNoise{run_index}.dat

run the program as follows: 
```sh
./cbf run_index is_af cow_var hr_0
```

the arguments are as follows:

- run_index: the index of the simulation instance. Determines the name of input files used and output files generated. Accepted values are positive integers.
- is_af: determines whether simulation is run with healthy or af conditions. Accepted values are 0 or 1.
- cow_var: determines which cow variant is represented in the model. Accepted values are 0 - 5.
- hr_0: Intrinsic heart rate of the simulation. Accepted values are positive integers. Note: model may behave unusually if values are out of physiological range.


The code is originally developed in matlab by Tim Hunter and SRK, based on the Ursino and Heldt models, the code for which is openly available.
