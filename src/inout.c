/*

       AAAA    CCCC   OOOO   TTTTTT   SSSSS  PPPPP
      AA  AA  CC     OO  OO    TT    SS      PP  PP
      AAAAAA  CC     OO  OO    TT     SSSS   PPPPP
      AA  AA  CC     OO  OO    TT        SS  PP
      AA  AA   CCCC   OOOO     TT    SSSSS   PP

######################################################
##########    ACO algorithms for the TSP    ##########
######################################################

      Version: 1.0
      File:    InOut.c
      Author:  Thomas Stuetzle
      Purpose: mainly input / output / statistic routines
      Check:   README and gpl.txt
      Copyright (C) 2002  Thomas Stuetzle
 */

/***************************************************************************

    Program's name: acotsp

    Ant Colony Optimization algorithms (AS, ACS, EAS, RAS, MMAS, BWAS) for the 
    symmetric TSP 

    Copyright (C) 2004  Thomas Stuetzle

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    email: stuetzle no@spam ulb.ac.be
    mail address: Universite libre de Bruxelles
                  IRIDIA, CP 194/6
                  Av. F. Roosevelt 50
                  B-1050 Brussels
          Belgium

 ***************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include "inout.h"
#include "thop.h"
#include "timer.h"
#include "utilities.h"
#include "ants.h"
#include "ls.h"
#include "parse.h"

long int *best_in_try;
long int *best_found_at;
double *time_best_found;
double *time_total_run;

long int n_try; /* try counter */
long int n_tours; /* counter of number constructed tours */
long int iteration; /* iteration counter */
long int restart_iteration; /* remember iteration when restart was done if any */
double restart_time; /* remember time when restart was done if any */
long int max_tries; /* maximum number of independent tries */
long int max_tours; /* maximum number of tour constructions in one try */
long int max_packing_tries; /* number of tries to construct a good packing plan from a give tour */
long int seed;

double lambda; /* Parameter to determine branching factor */
double branch_fac; /* If branching factor < branch_fac => update trails */

double max_time; /* maximal allowed run time of a try  */
double time_used; /* time used until some given event */
double time_passed; /* time passed until some moment*/
long int optimal; /* optimal solution or bound to find */

double mean_ants; /* average tour length */
double stddev_ants; /* stddev of tour lengths */
double branching_factor; /* average node branching factor when searching */
double found_branching; /* branching factor when best solution is found */

long int found_best; /* iteration in which best solution is found */
long int restart_found_best; /* iteration in which restart-best solution is found */

/* ------------------------------------------------------------------------ */

FILE *log_file;

char input_name_buf[LINE_BUF_LEN];
char output_name_buf[LINE_BUF_LEN];
int opt;
long int log_flag; /* --log was given in the command-line.  */
long int output_flag; 
long int calibration_mode;

void init_program(long int argc, char * argv[])
/*    
      FUNCTION:       initialize the program, 
      INPUT:          program arguments, needed for parsing commandline
      OUTPUT:         none
      COMMENTS:       
 */
{

    char temp_buffer[LINE_BUF_LEN];

    /*printf(PROG_ID_STR);*/
    
    set_default_parameters();
    setbuf(stdout, NULL);
    parse_commandline(argc, argv);

    assert(max_tries <= MAXIMUM_NO_TRIES);

    best_in_try = calloc(max_tries, sizeof(long int));
    best_found_at = calloc(max_tries, sizeof(long int));
    time_best_found = calloc(max_tries, sizeof(double));
    time_total_run = calloc(max_tries, sizeof(double));

    TRACE(printf("read problem data  ..\n\n");)
    read_thop_instance(input_name_buf, &instance.nodeptr, &instance.itemptr);
    TRACE(printf("\n .. done\n\n");)
    
    if ( max_time < 0 ) {
        /* Change default parameter max_time for ceil(number of items * 0.1) */
        max_time = ceil(instance.m / 10.0);
    }

    if (n_ants < 0) n_ants = instance.n;
    /* default setting for elitist_ants is 0; if EAS is applied and
       option elitist_ants is not used, we set the default to
       elitist_ants = instance.n */
    if (eas_flag && elitist_ants <= 0) elitist_ants = instance.n;

    nn_ls = MIN(instance.n - 1, nn_ls);

    assert(n_ants < MAX_ANTS - 1);
    assert(nn_ants < MAX_NEIGHBOURS);
    assert(nn_ants > 0);
    assert(nn_ls > 0);

    if (!log_flag) {
        sprintf(temp_buffer, "%s.log", output_name_buf);
        log_file = fopen(temp_buffer, "w");
    } else {
        log_file = NULL;
    }

    instance.distance = compute_distances();
    
    write_params();
    
    allocate_ants();
}

void exit_program(void)
/*    
      FUNCTION:       save some final statistical information on a trial once it finishes
      INPUT:          none
      OUTPUT:         none
      COMMENTS:       
 */
{
    if (log_file) {        
        fprintf(log_file,"\n\n");
        long int ntry = 0;
        for(; ntry < max_tries; ntry++) {
            fprintf(log_file, "try %10ld,        best %10ld,        found at iteration %10ld,        found at time %10.2f\n", ntry, instance.UB + 1 - best_in_try[ntry], best_found_at[ntry], time_best_found[ntry]);
            fflush(log_file);
        }
    }
    
    long int profit = instance.UB + 1 - global_best_ant->fitness;

    if ( calibration_mode ) printf("%ld\n", -profit);    
    else printf("Best solution: %ld\n", profit);
    
    if (output_flag) save_best_thop_solution();
}

void init_try( long int ntry ) 
/*    
      FUNCTION: initilialize variables appropriately when starting a trial
      INPUT:    trial number
      OUTPUT:   none
      COMMENTS: none
 */
{

    TRACE ( printf("INITIALIZE TRIAL\n"); );

    start_timers();
    time_used = elapsed_time( VIRTUAL );
    time_passed = time_used;

    /* Initialize variables concerning statistics etc. */

    n_tours      = 1;
    iteration    = 1;
    restart_iteration = 1;
    lambda       = 0.05;
    best_so_far_ant->fitness = INFTY;
    found_best   = 0;

    /* Initialize the Pheromone trails, only if ACS is used, pheromones
       have to be initialized differently */
    if ( !(acs_flag || mmas_flag || bwas_flag) ) {
        trail_0 = 1. / ( (rho) * nn_tour() );
        /* in the original papers on Ant System, Elitist Ant System, and
           Rank-based Ant System it is not exactly defined what the
           initial value of the pheromones is. Here we set it to some
           small constant, analogously as done in MAX-MIN Ant System.  
         */
        init_pheromone_trails( trail_0 );
    }
    if ( bwas_flag ) {
        trail_0 = 1. / ( (double) instance.n * (double) nn_tour()) ;
        init_pheromone_trails( trail_0 );
    }
    if ( mmas_flag ) {
        trail_max = 1. / ( (rho) * nn_tour() );
        trail_min = trail_max / ( 2. * instance.n );
        init_pheromone_trails( trail_max );
    }
    if ( acs_flag ) {
        trail_0 = 1. / ( (double) instance.n * (double) nn_tour( ) ) ;
        init_pheromone_trails( trail_0 );
    }

    /* Calculate combined information pheromone times heuristic information */
    compute_total_information();

    if (log_file) fprintf(log_file,"\nbegin try %li \n",ntry);
}

void exit_try(long int ntry)
/*    
      FUNCTION:       save some statistical information on a trial once it finishes
      INPUT:          trial number
      OUTPUT:         none
      COMMENTS:       
 */
{    
    best_in_try[ntry] = best_so_far_ant->fitness;
    best_found_at[ntry] = found_best;
    time_best_found[ntry] = time_used;
    time_total_run[ntry] = elapsed_time(VIRTUAL);
    
    if (best_so_far_ant->fitness < global_best_ant->fitness) {
        copy_from_to( best_so_far_ant, global_best_ant );
    }
        
    if (log_file) fprintf(log_file,"end try %li \n",ntry);
}

void read_thop_instance(const char *input_file_name, struct point **nodeptr, struct item **itemptr)
/*    
      FUNCTION: parse and read instance file
      INPUT:    instance name
      OUTPUT:   list of coordinates for all nodes
      COMMENTS: Instance files have to be in TSPLIB format, otherwise procedure fails
 */
{
    FILE *input_file;
    char buf[LINE_BUF_LEN];
    long int i, j, k;

    input_file = fopen(input_file_name, "r");
    if (input_file == NULL) {
        fprintf(stderr, "No instance file specified, abort\n");
        exit(1);
    }
    assert(input_file != NULL);
    /*printf("\nreading thop-file %s ... \n\n", input_file_name);*/

    fscanf(input_file, "PROBLEM NAME: %s\n", buf);
    fscanf(input_file, "KNAPSACK DATA TYPE: %[^\n]\n", instance.knapsack_data_type);
    fscanf(input_file, "DIMENSION: %ld\n", &instance.n); ++instance.n;
    assert(instance.n > 3 && instance.n < 6000);
    fscanf(input_file, "NUMBER OF ITEMS: %ld\n", &instance.m);
    fscanf(input_file, "CAPACITY OF KNAPSACK: %ld\n", &instance.capacity_of_knapsack);
    fscanf(input_file, "MAX TIME: %lf\n", &instance.max_time);
    fscanf(input_file, "MIN SPEED: %lf\n", &instance.min_speed);
    fscanf(input_file, "MAX SPEED: %lf\n", &instance.max_speed);
    fscanf(input_file, "EDGE_WEIGHT_TYPE: %s\n", buf);
    if (strcmp("EUC_2D", buf) == 0) distance = round_distance;
    else if (strcmp("CEIL_2D", buf) == 0) distance = ceil_distance;
    else if (strcmp("GEO", buf) == 0) distance = geo_distance;
    else if (strcmp("ATT", buf) == 0) distance = att_distance;
    fgets(buf, LINE_BUF_LEN, input_file); /* NODE_COORD_SECTION  (INDEX, X, Y): */

    if (( * nodeptr = malloc(sizeof(struct point) * (instance.n))) == NULL)
        exit(EXIT_FAILURE);
    else {
        for (i = 0; i < instance.n - 1; i++) {
            fscanf(input_file, "%ld %lf %lf\n", & j, &(*nodeptr)[i].x, &(*nodeptr)[i].y);
        }
    }
    TRACE(printf("number of cities is %ld\n", n);)

    fgets(buf, LINE_BUF_LEN, input_file); /* ITEMS SECTION    (INDEX, PROFIT, WEIGHT, ASSIGNED NODE NUMBER): */

    if (( * itemptr = malloc(sizeof(struct item) * instance.m)) == NULL)
        exit(EXIT_FAILURE);
    else {
        for (i = 0; i < instance.m; i++) {
            fscanf(input_file, "%ld %ld %ld %ld\n", &j, &(*itemptr)[i].profit, &(*itemptr)[i].weight, &(*itemptr)[i].id_city);
            (*itemptr)[i].id_city -= 1;
        }
    }
    
    double *item_vector = malloc(instance.m * sizeof(double));
    long int *help_vector = malloc(instance.m * sizeof(long int));
    
    for ( j = 0 ; j < instance.m ; j++ ) {
        item_vector[j] = ( -1.0 * (*itemptr)[j].profit ) / (*itemptr)[j].weight;
        help_vector[j] = j;
    }

    sort2_double(item_vector, help_vector, 0, instance.m - 1);
    
    instance.UB = 0;
    long int _w = 0;
    for ( k = 0 ; k < instance.m ; k++ ) {
        j = help_vector[k];
        if ( _w + (*itemptr)[j].weight <= instance.capacity_of_knapsack ) {
            _w += (*itemptr)[j].weight;
            instance.UB += (*itemptr)[j].profit;
        }
        else {
            instance.UB += ceil((instance.capacity_of_knapsack - _w) / (double) (*itemptr)[j].weight * (*itemptr)[j].profit);
            break;
        }
    }
    
    free(item_vector);
    free(help_vector);

    TRACE(printf("number of items is %ld\n", instance.m);)
    TRACE(printf("\n... done\n");)
    
    fclose(input_file);
}

void set_default_parameters(void)
/*    
      FUNCTION: set default parameter settings
      INPUT:    none
      OUTPUT:   none
      COMMENTS: none
 */
{
    ls_flag = 0; /* per default run with no local search*/
    dlb_flag = TRUE; /* apply don't look bits in local search */
    nn_ls = 20; /* use fixed radius search in the 20 nearest neighbours */
    n_ants = 25; /* number of ants */
    nn_ants = 20; /* number of nearest neighbours in tour construction */
    alpha = 1.0;
    beta = 2.0;
    rho = 0.5;
    q_0 = 0.0;
    max_tries = 1;
    max_tours = 0;
    max_packing_tries = 1;
    seed = (long int) time(NULL);
    max_time = -1;
    optimal = 1;
    branch_fac = 1.00001;
    u_gb = INFTY;
    as_flag = FALSE;
    eas_flag = FALSE;
    ras_flag = FALSE;
    mmas_flag = TRUE;
    bwas_flag = FALSE;
    acs_flag = FALSE;
    ras_ranks = 0;
    elitist_ants = 0;
}

void set_default_as_parameters(void) {
    assert(as_flag);
    n_ants = -1; /* number of ants (-1 means instance size) */
    nn_ants = 20; /* number of nearest neighbours in tour construction */
    alpha = 1.0;
    beta = 2.0;
    rho = 0.5;
    q_0 = 0.0;
    ras_ranks = 0;
    elitist_ants = 0;
}

void set_default_eas_parameters(void) {
    assert(eas_flag);
    n_ants = -1; /* number of ants (-1 means instance size) */
    nn_ants = 20; /* number of nearest neighbours in tour construction */
    alpha = 1.0;
    beta = 2.0;
    rho = 0.5;
    q_0 = 0.0;
    ras_ranks = 0;
    elitist_ants = n_ants;
}

void set_default_ras_parameters(void) {
    assert(ras_flag);
    n_ants = -1; /* number of ants (-1 means instance size) */
    nn_ants = 20; /* number of nearest neighbours in tour construction */
    alpha = 1.0;
    beta = 2.0;
    rho = 0.1;
    q_0 = 0.0;
    ras_ranks = 6;
    elitist_ants = 0;
}

void set_default_bwas_parameters(void) {
    assert(bwas_flag);
    n_ants = -1; /* number of ants (-1 means instance size) */
    nn_ants = 20; /* number of nearest neighbours in tour construction */
    alpha = 1.0;
    beta = 2.0;
    rho = 0.1;
    q_0 = 0.0;
    ras_ranks = 0;
    elitist_ants = 0;
}

void set_default_mmas_parameters(void) {
    assert(mmas_flag);
    n_ants = -1; /* number of ants (-1 means instance size) */
    nn_ants = 20; /* number of nearest neighbours in tour construction */
    alpha = 1.0;
    beta = 2.0;
    rho = 0.02;
    q_0 = 0.0;
    ras_ranks = 0;
    elitist_ants = 0;
}

void set_default_acs_parameters(void) {
    assert(acs_flag);
    n_ants = 10; /* number of ants (-1 means instance size) */
    nn_ants = 20; /* number of nearest neighbours in tour construction */
    alpha = 1.0;
    beta = 2.0;
    rho = 0.1;
    q_0 = 0.9;
    ras_ranks = 0;
    elitist_ants = 0;
}

void set_default_ls_parameters(void) {
    assert(ls_flag);
    dlb_flag = TRUE; /* apply don't look bits in local search */
    nn_ls = 20; /* use fixed radius search in the 20 nearest neighbours */
    
    n_ants = 25; /* number of ants */
    nn_ants = 20; /* number of nearest neighbours in tour construction */
    alpha = 1.0;
    beta = 2.0;
    rho = 0.5;
    q_0 = 0.0;

    if (mmas_flag) {
        n_ants = 25;
        rho = 0.2;
        q_0 = 0.00;
    } else if (acs_flag) {
        n_ants = 10;
        rho = 0.1;
        q_0 = 0.98;
    } else if (eas_flag) {
        elitist_ants = n_ants;
    }
}

void save_best_thop_solution(void) 
{
    int i, first_print;
    char *visited = calloc(instance.n, sizeof(char));
    visited[0] = visited[instance.n - 2] = 1;

    long int profit = 0.0;

    for (i = 0; i < instance.m; i++) {
        if ( global_best_ant->packing_plan[i] ) {
            visited[instance.itemptr[i].id_city] = 1;
            profit += instance.itemptr[i].profit;
        }
    }
    
    FILE *sol_file = fopen(output_name_buf, "w");
    
    first_print = TRUE;
    fprintf(sol_file, "[");
    for (i = 1; i < instance.n - 2; i++) {
        if ( visited[global_best_ant->tour[i]] ) {
            if ( first_print == TRUE ) {
                first_print = FALSE;
                fprintf(sol_file, "%d", global_best_ant->tour[i] + 1);
            }
            else fprintf(sol_file, ",%d", global_best_ant->tour[i] + 1);
        }
    }
    fprintf(sol_file, "]\n[");

    first_print = TRUE;
    for (i = 0; i < instance.m; i++) {
        if ( global_best_ant->packing_plan[i] ) {
            if ( first_print == TRUE ) {
                first_print = FALSE;            
                fprintf(sol_file, "%d", i+1);
            }
            else fprintf(sol_file, ",%d", i+1);
        }
    }
    fprintf(sol_file, "]\n");
    free(visited);
    
    fclose(sol_file);
}

void write_report(void)
/*    
      FUNCTION: output some info about trial (best-so-far solution quality, time)
      INPUT:    none
      OUTPUT:   none
      COMMENTS: none
 */
{
    if (log_file) {
        fprintf(log_file, "best %10ld,        iteration: %10ld,        time %10.2f\n", instance.UB + 1 - best_so_far_ant->fitness, iteration, elapsed_time(VIRTUAL));
        fflush(log_file);
    }
}


void write_params(void)
/*    
      FUNCTION:       writes chosen parameter settings in log file 
      INPUT:          none
      OUTPUT:         none
 */
{
    if (log_file) {           
        fprintf(log_file, "Parameter-settings: \n\n");
        fprintf(log_file, "--inputfile          %s\n", input_name_buf);
        fprintf(log_file, "--outputfile         %s\n", output_name_buf);
        fprintf(log_file, "--tries              %ld\n", max_tries);
        fprintf(log_file, "--tours              %ld\n", max_tours);
        fprintf(log_file, "--ptries             %ld\n", max_packing_tries);    
        fprintf(log_file, "--time               %.2f\n", max_time);
        fprintf(log_file, "--seed               %ld\n", seed);
        fprintf(log_file, "--optimum            %ld\n", optimal);
        fprintf(log_file, "--ants               %ld\n", n_ants);
        fprintf(log_file, "--nnants             %ld\n", nn_ants);
        fprintf(log_file, "--alpha              %.2f\n", alpha);
        fprintf(log_file, "--beta               %.2f\n", beta);
        fprintf(log_file, "--rho                %.2f\n", rho);
        fprintf(log_file, "--q0                 %.2f\n", q_0);
        fprintf(log_file, "--elitistants        %ld\n", elitist_ants);
        fprintf(log_file, "--rasranks           %ld\n", ras_ranks);
        fprintf(log_file, "--localsearch        %ld\n", ls_flag);
        fprintf(log_file, "--nnls               %ld\n", nn_ls);
        fprintf(log_file, "--dlb                %ld\n", dlb_flag);
        fprintf(log_file, "--as                 %ld\n", as_flag);
        fprintf(log_file, "--eas                %ld\n", eas_flag);
        fprintf(log_file, "--ras                %ld\n", ras_flag);
        fprintf(log_file, "--mmas               %ld\n", mmas_flag);
        fprintf(log_file, "--bwas               %ld\n", bwas_flag);
        fprintf(log_file, "--acs                %ld\n\n", acs_flag);
    }
}
