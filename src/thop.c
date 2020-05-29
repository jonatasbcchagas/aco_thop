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
      File:    TSP.c
      Author:  Thomas Stuetzle
      Purpose: TSP related procedures, distance computation, neighbour lists
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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <assert.h>

#include "inout.h"
#include "thop.h"
#include "ants.h"
#include "ls.h"
#include "utilities.h"

#define M_PI 3.14159265358979323846264

struct problem instance;

static double dtrunc (double x)
{
    int k;

    k = (int) x;
    x = (double) k;
    return x;
}

long int  (*distance)(long int, long int);  /* function pointer */

/*    
      FUNCTION: the following four functions implement different ways of 
                computing distances for TSPLIB instances
      INPUT:    two node indices
      OUTPUT:   distance between the two nodes
 */

long int round_distance (long int i, long int j) 
/*    
      FUNCTION: compute Euclidean distances between two nodes rounded to next 
                integer for TSPLIB instances
      INPUT:    two node indices
      OUTPUT:   distance between the two nodes
      COMMENTS: for the definition of how to compute this distance see TSPLIB
 */
{
    double xd = instance.nodeptr[i].x - instance.nodeptr[j].x;
    double yd = instance.nodeptr[i].y - instance.nodeptr[j].y;
    double r  = sqrt(xd*xd + yd*yd) + 0.5;

    return (long int) r;
}

long int ceil_distance (long int i, long int j) 
/*    
      FUNCTION: compute ceiling distance between two nodes rounded to next 
                integer for TSPLIB instances
      INPUT:    two node indices
      OUTPUT:   distance between the two nodes
      COMMENTS: for the definition of how to compute this distance see TSPLIB
 */
{
    double xd = instance.nodeptr[i].x - instance.nodeptr[j].x;
    double yd = instance.nodeptr[i].y - instance.nodeptr[j].y;
    double r  = sqrt(xd*xd + yd*yd);

    return (long int)(ceil (r));
}

long int geo_distance (long int i, long int j) 
/*    
      FUNCTION: compute geometric distance between two nodes rounded to next 
                integer for TSPLIB instances
      INPUT:    two node indices
      OUTPUT:   distance between the two nodes
      COMMENTS: adapted from concorde code
                for the definition of how to compute this distance see TSPLIB
 */
{
    double deg, min;
    double lati, latj, longi, longj;
    double q1, q2, q3;
    long int dd;
    double x1 = instance.nodeptr[i].x, x2 = instance.nodeptr[j].x,
            y1 = instance.nodeptr[i].y, y2 = instance.nodeptr[j].y;

    deg = dtrunc (x1);
    min = x1 - deg;
    lati = M_PI * (deg + 5.0 * min / 3.0) / 180.0;
    deg = dtrunc (x2);
    min = x2 - deg;
    latj = M_PI * (deg + 5.0 * min / 3.0) / 180.0;

    deg = dtrunc (y1);
    min = y1 - deg;
    longi = M_PI * (deg + 5.0 * min / 3.0) / 180.0;
    deg = dtrunc (y2);
    min = y2 - deg;
    longj = M_PI * (deg + 5.0 * min / 3.0) / 180.0;

    q1 = cos (longi - longj);
    q2 = cos (lati - latj);
    q3 = cos (lati + latj);
    dd = (int) (6378.388 * acos (0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0);
    return dd;

}

long int att_distance (long int i, long int j) 
/*    
      FUNCTION: compute ATT distance between two nodes rounded to next 
                integer for TSPLIB instances
      INPUT:    two node indices
      OUTPUT:   distance between the two nodes
      COMMENTS: for the definition of how to compute this distance see TSPLIB
 */
{
    double xd = instance.nodeptr[i].x - instance.nodeptr[j].x;
    double yd = instance.nodeptr[i].y - instance.nodeptr[j].y;
    double rij = sqrt ((xd * xd + yd * yd) / 10.0);
    double tij = dtrunc (rij);
    long int dij;

    if (tij < rij)
        dij = (int) tij + 1;
    else
        dij = (int) tij;
    return dij;
}



long int** compute_distances(void)
/*    
      FUNCTION: computes the matrix of all intercity distances
      INPUT:    none
      OUTPUT:   pointer to distance matrix, has to be freed when program stops
 */
{
    long int     i, j;
    long int     **matrix;

    if((matrix = malloc(sizeof(long int) * instance.n * instance.n + sizeof(long int *) * instance.n)) == NULL){
        fprintf(stderr,"Out of memory, exit.");
        exit(1);
    }
    for ( i = 0 ; i < instance.n ; i++ ) {
        matrix[i] = (long int *)(matrix + instance.n) + i * instance.n;
    }

    int max_distance = 0;
    for ( i = 0 ; i < instance.n - 1 ; i++ ) {
        for ( j = 0  ; j < instance.n - 1 ; j++ ) {
            matrix[i][j] = distance(i, j);
            if ( matrix[i][j] > max_distance ) max_distance = matrix[i][j];
        }
    }

    for ( i = 0; i < instance.n ; i++ ) {
        matrix[i][instance.n - 1] = matrix[instance.n - 1][i] = max_distance * ( instance.n - 1 );
    }
    matrix[0][instance.n - 1] = matrix[instance.n - 1][0] = 0;
    matrix[instance.n - 2][instance.n - 1] = matrix[instance.n - 1][instance.n - 2] = 0;

    return matrix;
}


long int** compute_nn_lists( void )
/*    
      FUNCTION: computes nearest neighbor lists of depth nn for each city
      INPUT:    none
      OUTPUT:   pointer to the nearest neighbor lists
 */
{
    long int i, node, nn;
    long int *distance_vector;
    long int *help_vector;
    long int **m_nnear;

    TRACE ( printf("\n computing nearest neighbor lists, "); )

    nn = MAX(nn_ls,nn_ants);
    if ( nn >= instance.n )
        nn = instance.n - 1;
    DEBUG ( assert( instance.n > nn ); )

    TRACE ( printf("nn = %ld ... \n",nn); )

    if((m_nnear = malloc(sizeof(long int) * instance.n * nn + instance.n * sizeof(long int *))) == NULL){
        exit(EXIT_FAILURE);
    }
    distance_vector = calloc(instance.n, sizeof(long int));
    help_vector = calloc(instance.n, sizeof(long int));

    for ( node = 0 ; node < instance.n ; node++ ) {  /* compute cnd-sets for all node */
        m_nnear[node] = (long int *)(m_nnear + instance.n) + node * nn;

        for ( i = 0 ; i < instance.n ; i++ ) {  /* Copy distances from nodes to the others */
            distance_vector[i] = instance.distance[node][i];
            help_vector[i] = i;
        }
        distance_vector[node] = LONG_MAX;  /* city is not nearest neighbour */
        sort2(distance_vector, help_vector, 0, instance.n - 1);
        for ( i = 0 ; i < nn ; i++ ) {
            m_nnear[node][i] = help_vector[i];
        }
    }
    free(distance_vector);
    free(help_vector);
    TRACE ( printf("\n    .. done\n"); )

    return m_nnear;
}

long int compute_fitness( long int *t,  char *p ) 
/*    
      FUNCTION: compute the fitness of the ThOP solution generated from tour t
      INPUT:    pointer to tour t and pointer to packing plan p
      OUTPUT:   fitness of the ThOP solution generated from tour t
 */
{
    int i, j, k;
    double par_a, par_b, par_c, par_sum;
    long int prev_city, curr_city;
    double _total_time;
    long int _total_weight, total_weight, total_profit;    
    int violate_max_time;
    
    const double v = ( instance.max_speed - instance.min_speed ) / instance.capacity_of_knapsack;
    
    long int *distance_accumulated = calloc (instance.n, sizeof(long int));

    long int total_distance = 0;

    for ( i = 0 ; i < instance.n ; i++ ) {
        distance_accumulated[t[i]] = total_distance;
        total_distance += instance.distance[t[i]][t[i+1]];        
    }
    
    double *item_vector = malloc(instance.m * sizeof(double));
    long int *help_vector = malloc(instance.m * sizeof(long int));
    
    long int *profit_accumulated = malloc(instance.n * sizeof(long int));
    long int *weight_accumulated = malloc(instance.n * sizeof(long int));
    
    long int best_packing_plan_profit = 0;
    char *tmp_packing_plan = malloc(instance.m * sizeof(char));
    
    long int _try;
        
    for( _try = 0; _try < max_packing_tries; _try++) {
        
        for ( i = 0 ; i < instance.n ; i++ ) {
            profit_accumulated[i] = weight_accumulated[i] = 0;
        }
        
        par_a = ran01( &seed );  /* uniform random number between [0.0, 1.0] */
        par_b = ran01( &seed );  /* uniform random number between [0.0, 1.0] */
        par_c = ran01( &seed );  /* uniform random number between [0.0, 1.0] */

        par_sum = (par_a + par_b + par_c); 
        par_a /= par_sum; par_b /= par_sum; par_c /= par_sum;
                
        for ( j = 0 ; j < instance.m ; j++ ) {
            item_vector[j] = ( 
                                -1.0 * pow(instance.itemptr[j].profit, par_a)
                            ) 
                            / 
                            ( 
                                pow(instance.itemptr[j].weight, par_b) * pow((distance_accumulated[instance.n - 2] - distance_accumulated[instance.itemptr[j].id_city]), par_c)
                            );
            help_vector[j] = j;
        }

        sort2_double(item_vector, help_vector, 0, instance.m - 1);
    
        total_weight = 0, total_profit = 0;            
        
        for ( k = 0 ; k < instance.m ; k++ ) {

            j = help_vector[k];
            tmp_packing_plan[j] = 0;
            
            if ( total_weight + instance.itemptr[j].weight > instance.capacity_of_knapsack ) continue;

            profit_accumulated[instance.itemptr[j].id_city] += instance.itemptr[j].profit;
            weight_accumulated[instance.itemptr[j].id_city] += instance.itemptr[j].weight; 
            
            violate_max_time = FALSE;
            _total_time = _total_weight = 0;
            prev_city = 0;
            for ( i = 1 ; i < instance.n - 1; i++ ) {
                curr_city = t[i];
                if ( weight_accumulated[curr_city] == 0 && curr_city != instance.n - 2) continue;
                _total_time += instance.distance[prev_city][curr_city] / ( instance.max_speed - v * _total_weight );    
                if ( _total_time - EPSILON > instance.max_time ) {
                    violate_max_time = TRUE; break;
                }                
                _total_weight += weight_accumulated[curr_city];
                prev_city = curr_city;
            }

            if ( violate_max_time == FALSE) {
                total_profit += instance.itemptr[j].profit;     
                total_weight += instance.itemptr[j].weight;
                tmp_packing_plan[j] = 1;
            }
            else {
                profit_accumulated[instance.itemptr[j].id_city] -= instance.itemptr[j].profit;
                weight_accumulated[instance.itemptr[j].id_city] -= instance.itemptr[j].weight; 
            }
        }
        
        if ( total_profit > best_packing_plan_profit) {
            best_packing_plan_profit = total_profit;
            for ( j = 0 ; j < instance.m ; j++ ) {
                p[j] = tmp_packing_plan[j];
            }            
        }
    }
    
    free(distance_accumulated);
    free(item_vector);
    free(help_vector);
    free(profit_accumulated);
    free(weight_accumulated);
    free(tmp_packing_plan);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
    return instance.UB + 1 - best_packing_plan_profit;
}
