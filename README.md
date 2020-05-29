# Ants can orienteer a thief in their robbery

This project contains the code of the Swarm Intelligence Based on Ant Colony Optimization (ACO), which is described in the paper ["Ants can orienteer a thief in their robbery - Jonatas B. C. Chagas and Markus Wagner"](https://), for solving the Thief Orienteering Problem (ThOP).

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
