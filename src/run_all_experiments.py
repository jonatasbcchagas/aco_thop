#!/usr/bin/python
# -*- coding: utf-8 -*-

import itertools
import os
import multiprocessing
import argparse
import math

parameter_configurations = {
"eil51_01_bsc": {"--ants":  "100", "--alpha": "1.03", "--beta": "2.64", "--rho": "0.82", "--ptries": "2"}, 
"eil51_01_unc": {"--ants":  "200", "--alpha": "0.91", "--beta": "4.00", "--rho": "0.48", "--ptries": "1"}, 
"eil51_01_usw": {"--ants": "1000", "--alpha": "0.74", "--beta": "4.22", "--rho": "0.61", "--ptries": "3"}, 
"eil51_03_bsc": {"--ants":   "50", "--alpha": "1.00", "--beta": "3.81", "--rho": "0.37", "--ptries": "1"}, 
"eil51_03_unc": {"--ants":  "200", "--alpha": "0.93", "--beta": "4.12", "--rho": "0.69", "--ptries": "2"}, 
"eil51_03_usw": {"--ants":   "20", "--alpha": "1.00", "--beta": "3.87", "--rho": "0.49", "--ptries": "1"}, 
"eil51_05_bsc": {"--ants":  "200", "--alpha": "0.74", "--beta": "3.30", "--rho": "0.65", "--ptries": "4"}, 
"eil51_05_unc": {"--ants":   "50", "--alpha": "0.98", "--beta": "3.65", "--rho": "0.44", "--ptries": "2"}, 
"eil51_05_usw": {"--ants": "1000", "--alpha": "0.79", "--beta": "4.36", "--rho": "0.83", "--ptries": "1"}, 
"eil51_10_bsc": {"--ants":   "20", "--alpha": "1.15", "--beta": "4.94", "--rho": "0.88", "--ptries": "2"}, 
"eil51_10_unc": {"--ants":  "100", "--alpha": "0.86", "--beta": "4.34", "--rho": "0.74", "--ptries": "1"}, 
"eil51_10_usw": {"--ants":   "50", "--alpha": "0.93", "--beta": "6.15", "--rho": "0.69", "--ptries": "1"}, 
"pr107_01_bsc": {"--ants":  "200", "--alpha": "0.84", "--beta": "3.91", "--rho": "0.94", "--ptries": "5"}, 
"pr107_01_unc": {"--ants":  "100", "--alpha": "1.11", "--beta": "3.16", "--rho": "0.67", "--ptries": "2"}, 
"pr107_01_usw": {"--ants":   "50", "--alpha": "0.94", "--beta": "3.72", "--rho": "0.61", "--ptries": "3"}, 
"pr107_03_bsc": {"--ants":  "200", "--alpha": "0.79", "--beta": "3.50", "--rho": "0.28", "--ptries": "3"}, 
"pr107_03_unc": {"--ants":  "200", "--alpha": "0.78", "--beta": "3.31", "--rho": "0.74", "--ptries": "5"}, 
"pr107_03_usw": {"--ants":   "20", "--alpha": "0.84", "--beta": "3.73", "--rho": "0.62", "--ptries": "2"}, 
"pr107_05_bsc": {"--ants": "1000", "--alpha": "0.88", "--beta": "4.39", "--rho": "0.72", "--ptries": "3"}, 
"pr107_05_unc": {"--ants":  "500", "--alpha": "0.84", "--beta": "4.00", "--rho": "0.64", "--ptries": "1"}, 
"pr107_05_usw": {"--ants":  "200", "--alpha": "0.88", "--beta": "2.75", "--rho": "0.88", "--ptries": "4"}, 
"pr107_10_bsc": {"--ants":  "200", "--alpha": "0.75", "--beta": "3.36", "--rho": "0.70", "--ptries": "5"}, 
"pr107_10_unc": {"--ants":  "200", "--alpha": "0.97", "--beta": "3.98", "--rho": "0.69", "--ptries": "1"}, 
"pr107_10_usw": {"--ants":  "200", "--alpha": "0.83", "--beta": "3.76", "--rho": "0.49", "--ptries": "1"}, 
"a280_01_bsc": {"--ants":  "200", "--alpha": "0.75", "--beta": "6.28", "--rho": "0.61", "--ptries": "1"}, 
"a280_01_unc": {"--ants":  "500", "--alpha": "0.86", "--beta": "7.13", "--rho": "0.38", "--ptries": "1"}, 
"a280_01_usw": {"--ants":   "50", "--alpha": "0.96", "--beta": "4.15", "--rho": "0.48", "--ptries": "2"}, 
"a280_03_bsc": {"--ants":   "50", "--alpha": "0.75", "--beta": "6.16", "--rho": "0.40", "--ptries": "2"}, 
"a280_03_unc": {"--ants":  "200", "--alpha": "0.85", "--beta": "7.62", "--rho": "0.35", "--ptries": "1"}, 
"a280_03_usw": {"--ants":  "500", "--alpha": "0.88", "--beta": "5.28", "--rho": "0.48", "--ptries": "1"}, 
"a280_05_bsc": {"--ants":  "200", "--alpha": "0.82", "--beta": "7.91", "--rho": "0.37", "--ptries": "1"}, 
"a280_05_unc": {"--ants":  "200", "--alpha": "0.85", "--beta": "9.32", "--rho": "0.49", "--ptries": "1"}, 
"a280_05_usw": {"--ants":  "200", "--alpha": "0.93", "--beta": "4.86", "--rho": "0.26", "--ptries": "1"}, 
"a280_10_bsc": {"--ants":  "100", "--alpha": "0.74", "--beta": "6.34", "--rho": "0.60", "--ptries": "1"}, 
"a280_10_unc": {"--ants":   "50", "--alpha": "0.82", "--beta": "8.63", "--rho": "0.62", "--ptries": "1"}, 
"a280_10_usw": {"--ants":  "200", "--alpha": "0.78", "--beta": "6.79", "--rho": "0.58", "--ptries": "1"}, 
"dsj1000_01_bsc": {"--ants":  "500", "--alpha": "2.10", "--beta": "8.95", "--rho": "0.14", "--ptries": "1"}, 
"dsj1000_01_unc": {"--ants":  "200", "--alpha": "0.91", "--beta": "6.69", "--rho": "0.39", "--ptries": "1"}, 
"dsj1000_01_usw": {"--ants":  "500", "--alpha": "0.79", "--beta": "8.22", "--rho": "0.51", "--ptries": "1"}, 
"dsj1000_03_bsc": {"--ants":  "500", "--alpha": "2.21", "--beta": "6.67", "--rho": "0.27", "--ptries": "1"}, 
"dsj1000_03_unc": {"--ants":  "100", "--alpha": "2.74", "--beta": "5.89", "--rho": "0.16", "--ptries": "4"}, 
"dsj1000_03_usw": {"--ants":  "100", "--alpha": "0.82", "--beta": "8.71", "--rho": "0.35", "--ptries": "1"}, 
"dsj1000_05_bsc": {"--ants":  "100", "--alpha": "7.55", "--beta": "6.20", "--rho": "0.06", "--ptries": "3"}, 
"dsj1000_05_unc": {"--ants":  "200", "--alpha": "5.17", "--beta": "7.26", "--rho": "0.07", "--ptries": "2"}, 
"dsj1000_05_usw": {"--ants":  "100", "--alpha": "0.96", "--beta": "6.87", "--rho": "0.27", "--ptries": "1"}, 
"dsj1000_10_bsc": {"--ants":   "50", "--alpha": "1.12", "--beta": "7.27", "--rho": "0.25", "--ptries": "3"}, 
"dsj1000_10_unc": {"--ants":  "200", "--alpha": "0.88", "--beta": "7.34", "--rho": "0.43", "--ptries": "1"}, 
"dsj1000_10_usw": {"--ants":  "100", "--alpha": "0.99", "--beta": "7.99", "--rho": "0.34", "--ptries": "2"}, 
"general_parameter_configuration": {"--ants": "196", "--alpha": "1.24", "--beta": "5.46", "--rho": "0.51", "--ptries": "1"}
}

random_seeds = [ 269070,  99470, 126489, 644764, 547617, 642580,  73456, 462018, 858990, 756112, 
                 701531, 342080, 613485, 131654, 886148, 909040, 146518, 782904,   3075, 974703, 
                 170425, 531298, 253045, 488197, 394197, 519912, 606939, 480271, 117561, 900952, 
                 968235, 345118, 750253, 420440, 761205, 130467, 928803, 768798, 640300, 871462, 
                 639622,  90614, 187822, 594363, 193911, 846042, 680779, 344008, 759862, 661168, 
                 223420, 959508,  62985, 349296, 910428, 964420, 422964, 384194, 985214,  57575, 
                 639619,  90505, 435236, 465842, 102567, 189997, 741017, 611828, 699223, 335142, 
                  52119,  49256, 324523, 348215, 651525, 517999, 830566, 958538, 880422, 390645, 
                 148265, 807740, 934464, 524847, 408760, 668587, 257030, 751580,  90477, 594476, 
                 571216, 306614, 308010, 661191, 890429, 425031,  69108, 435783,  17725, 335928, ]

def launcher(tsp_base, number_of_items_per_city, knapsack_type, knapsack_size, maximum_travel_time, repetition, best_parameter_configurations=True, runtime_factor="1x"):
    inputfile = "../instances/%s-thop/%s_%02d_%s_%02d_%02d.thop" % (tsp_base, tsp_base, number_of_items_per_city, knapsack_type, knapsack_size, maximum_travel_time)
    outputfile = "../solutions/%s/%s-thop/%s_%02d_%s_%02d_%02d_%02d.thop.sol" % \
                ("acothop*" if best_parameter_configurations else "acothop", \
                tsp_base, tsp_base, number_of_items_per_city, knapsack_type, knapsack_size, maximum_travel_time, repetition)
    parameter_configuration_key = "%s_%02d_%s" % (tsp_base, number_of_items_per_city, knapsack_type) if best_parameter_configurations else "general_parameter_configuration"

    os.system("./acothop --mmas --tries 1 --seed %d --time %.1f --inputfile %s --outputfile %s %s --log" % (random_seeds[repetition], \
                                                                                                        float(runtime_factor.replace('x','')) * math.ceil((int(''.join(filter(lambda x: x.isdigit(), tsp_base))) - 2) * number_of_items_per_city / 10.0), \
                                                                                                        inputfile, outputfile, ' '.join("%s %s" % (k, v) for k, v in parameter_configurations[parameter_configuration_key].items())))

if __name__ == "__main__":

    tsp_base = ["eil51", "pr107", "a280", "dsj1000", ]
    number_of_items_per_city = [1, 3, 5, 10, ]
    knapsack_type = ["bsc", "unc", "usw", ]
    knapsack_size = [1, 5, 10, ]
    maximum_travel_time = [1, 2, 3, ]
    
    os.system("make clean")
    os.system("make")

    pool = multiprocessing.Pool(processes=max(1, multiprocessing.cpu_count() - 2))
    
    for _product in itertools.product(tsp_base, number_of_items_per_city, knapsack_type, knapsack_size, maximum_travel_time):
        _tsp_base, _number_of_items_per_city, _knapsack_type, _knapsack_size, _maximum_travel_time = _product
        for repetition in range(10):
            pool.apply_async(launcher, args=(_tsp_base, _number_of_items_per_city, _knapsack_type, _knapsack_size, _maximum_travel_time, repetition))
            
    for _product in itertools.product(tsp_base, number_of_items_per_city, knapsack_type, knapsack_size, maximum_travel_time):
        _tsp_base, _number_of_items_per_city, _knapsack_type, _knapsack_size, _maximum_travel_time = _product
        for repetition in range(10):
            pool.apply_async(launcher, args=(_tsp_base, _number_of_items_per_city, _knapsack_type, _knapsack_size, _maximum_travel_time, repetition, False))        

    pool.close()
    pool.join()
