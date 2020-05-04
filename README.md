Nutrient Solution Simulator and plant simulator
========================

# AUTHOR
>Byunghyun Ban
* CTO of Imagination Garden Inc. (2018~)
* CFO & CEO of Cheesecake Studio Inc. (2016 ~ 2017)
* CEO & CTO of Studio Mic Inc. (2015)
* CTO of New Page Inc. (2011~2013)
* Master's Degree @ Bio and Brain Engineering Department, KAIST (Korea Advanced Institute of Science and Technology)
* bhban@kaist.ac.kr
* halfbottle@sangsang.farm


# Citation
* For nutrient solution simulation only

> B. Ban, M. Lee and D. Ryu, "ODE Network Model for Nonlinear and Complex Agricultural Nutrient Solution System," 2019 International Conference on Information and Communication Technology Convergence (ICTC), Jeju Island, Korea (South), 2019, pp. 996-1001.

* For plant-included simulation (please cite both of them)

> B. Ban, M. Lee and D. Ryu, "ODE Network Model for Nonlinear and Complex Agricultural Nutrient Solution System," 2019 International Conference on Information and Communication Technology Convergence (ICTC), Jeju Island, Korea (South), 2019, pp. 996-1001.

> Ban, B. Mathematical Model for Secondary Transport of Cations at the Root of Plants. Preprints 2019, 2019090219 (doi: 10.20944/preprints201909.0219.v1).


## 1. Environments
* Python3
* Don't support any problems from python2 environment

## 2. Dependencies
* numpy library required.
> pip install numpy on bash.


## 3 . Preperation of input file
Please read the sample system file
> nutrient_system_ksp.txt

Each terms in any line has an indicator.

k is for reaction rate component

& is for the name of a component

\* is for each term for differential equation expression

$ is for initial molar density

@ is for molecular weight

ion ends with #

cation absorption equation starts with =

cation absorption equation has Q, B, M, N.

Q is Da/SF^1.5

N is the number of electron an cation lost.

The ions inside the plant has indicator 'p' in front of names.


Differential equations which contains reaction rate constant k should follow after k indicators.

To hanlde a component which appears on differnt differential equations; you don't have to combine those equations into one line. This model automatically does it.

## 4. Load the system
You may just import the module on you source code.

>import system_simulator as ss

## 5. Read the input file

> network = ss.Network(filename)

## 6. Synchronous update
>network.synchronous_update()

## 7. Convergence
>network.converge()

## 8. Show Results
>network.show_result()

## 9. Write Results
>network.write_result(filename)

## 10. Export network file as Cytoscae form
>network.export_cytoscape

It does not work when cation absorption model is applied.

## 11. Show list of physiological terms
>network.byunghyun_coefficients()


## 12. Build Training Data and Test Data
>python data_formation.py <result_dir> <network_dir> <normalization> <ISE_observable>

result_dir : Directory where the result .csv files are stored.

network_dir : Directory where the network files used for the simulation are stored.

normalization : A boolean value. If True, it normalizes the data.

  X -= np.min(X)
  X /= np.max(X)

  Y = np.log10(Y)

  Y -= np.min(y)

  Y /= np.max(Y)

 ISE_observable : A boolean value. If false, all components are processed. If True, only the components below are processed as X values.

 > (NH4+, NO3-, K+, Ca++, pH, TDS_of_nutrient_solution)
