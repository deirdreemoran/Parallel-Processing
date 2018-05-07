/******************************************************************************
*  Author: Deirdre Moran
*  Program: project3_fix1.cpp
*  Date: 5/7/2018
*  Descripton: False sharing implementation
*  NOTE:  Code was partially created by Michael Bailey http://web.engr.oregonstate.edu/~mjb
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

struct s
{
		float value;
		int pad[NUM];
} Array[4];

int main( int argc, char *argv[ ] )
{
#ifndef _OPENMP
         fprintf( stderr, "OpenMP is not available\n" );
         return 1;
#endif
         omp_set_num_threads( NUMTHREADS );
         int someBigNumber = 1000000000;
		 double mcps;
 		 double time0 = omp_get_wtime( );

         #pragma omp parallel for
         for( int i = 0; i < 4; i++ )
         {
                 for( int j = 0; j < someBigNumber; j++ )
                 {
                          Array[ i ].value = Array[ i ].value + 2.;
                 }
        }
 		 double time1 = omp_get_wtime( );
         printf("Execution time: %f\n", time1 - time0);
		 mcps = ((float)(someBigNumber)*4/(time1-time0)/1000000.);
	 	 printf("MegaCounts Per Second %10.21lf\n", mcps);
//printf("Address = %p\n", (void*) Array );
}


