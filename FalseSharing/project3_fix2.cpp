/******************************************************************************
*  Author: Deirdre Moran
*  Program: project3_fix2.cpp
*  Date: 5/7/2018
*  Descripton: False sharing implementation, fix2
*  NOTE:  Code was partially created by Michael Bailey http://web.engr.oregonstate.edu/~mjb
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

struct s
{
		float value;
		//int pad[NUM];
} Array[4];

int main( int argc, char *argv[ ] )
{
#ifndef _OPENMP
         fprintf( stderr, "OpenMP is not available\n" );
         return 1;
#endif
         omp_set_num_threads( NUMTHREADS );
         double mcps;
         int someBigNumber = 1000000000;

 		 double time0 = omp_get_wtime( );

         #pragma omp parallel for
         for( int i = 0; i < 4; i++ )
         {
			 	// create localized temp variable in each core's stack,
			 	//so little to no cache line conflict
			 	float tmp = Array[i].value;
                 for( int j = 0; j < someBigNumber; j++ )
                 {
                          tmp = tmp + 2.;
                 }
        		 Array[i].value = tmp;
        }
     	 double time1 = omp_get_wtime( );
         printf("Execution time: %f\n", time1 - time0);
         mcps = ((float)(someBigNumber)* 4/(time1-time0)/1000000.);
	 	 printf("MegaCounts Per Second %10.21lf\n", mcps);
//printf("Address = %p\n", (void*) Array);
}


