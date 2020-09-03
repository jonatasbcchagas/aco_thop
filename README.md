# Ants can orienteer a thief in their robbery

This repository contains the source code and data associated to the paper ["Ants can orienteer a thief in their robbery"](https://www.sciencedirect.com/science/article/abs/pii/S0167637720301255) by Jonatas B. C. Chagas and Markus Wagner. The paper presents a Swarm Intelligence Based on Ant Colony Optimization (ACO) for solving the Thief Orienteering Problem (ThOP).

### Compiling the code

Before running our ACO algorithm, it is needed to compile its code. To this end, just run the following command:

```console
$ make
```

### Usage:

```console
$ ./acothop [parameters]

Parameters:

  -i, --inputfile       inputfile (ThOP format necessary)
  -o, --outputfile      outputfile
  -m, --ants            number of ants
  -a, --alpha           influence of pheromone trails
  -b, --beta            influence of heuristic information
  -e, --rho             pheromone trail evaporation
  -p, --ptries          number of tries to construct a packing plan from a give tour
  -t, --time            maximum time for each trial  
      --seed            seed for the random number generator
      --log             save an extra file (<outputfile>.log) with log messages

```

We provide a python script (see "run_all_experiments.py") for running all the computational experiments reported in our paper.
