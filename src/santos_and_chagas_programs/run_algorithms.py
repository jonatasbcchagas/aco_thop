#!/usr/bin/python
# -*- coding: utf-8 -*-

import itertools
import os


def launcher(algorithm, _tsp_base, _number_of_items_per_city, _knapsack_type, _knapsack_size, _maximum_travel_time):
    inputfile = "../instances/%s-thop/%s_%02d_%s_%02d_%02d.thop" % (_tsp_base, _tsp_base, _number_of_items_per_city, _knapsack_type, _knapsack_size, _maximum_travel_time)
    os.system("./%s_thop %s 10" % (algorithm, inputfile))

if __name__ == "__main__":

    algorithms = ["ils", "brkga", ]

    tsp_base = ["a280", "dsj1000", ]
    number_of_items_per_city = [1, 3, 5, 10, ]
    knapsack_type = ["bsc", "unc", "usw", ]
    knapsack_size = [1, 5, 10, ]
    maximum_travel_time = [1, 2, 3, ]
    
    for algorithm in algorithms:
        for _product in itertools.product(tsp_base, number_of_items_per_city, knapsack_type, knapsack_size, maximum_travel_time):
            _tsp_base, _number_of_items_per_city, _knapsack_type, _knapsack_size, _maximum_travel_time = _product
            launcher(algorithm, _tsp_base, _number_of_items_per_city, _knapsack_type, _knapsack_size, _maximum_travel_time)
