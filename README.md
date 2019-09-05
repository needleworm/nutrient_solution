Nutrient Solution Simulator
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
Ban, B., Lee, M., & Ryu, D. (2019). ODE network model for nonlinear and complex agricultural nutrient solution system. arXiv preprint arXiv:1907.10800.


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
