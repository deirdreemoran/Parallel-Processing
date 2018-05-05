/******************************************************************************
*  Author: Deirdre Moran
*  Program: project0.cpp
*  Date: 4/5/2018
*  Description:   Simple OpenMP Experiment - Large number arrays are
*				multiplied together to test parallel processing execution times.
*******************************************************************************/
#include <omp.h>
#include <stdio.h>
#include <math.h>

#define NUMT	        4  	// 1 or 4 threads
#define ARRAYSIZE       1000000
#define NUMTRIES        40

float A[ARRAYSIZE];
float B[ARRAYSIZE];
float C[ARRAYSIZE];

int main()
{
	#ifndef _OPENMP
		fprintf( stderr, "OpenMP is not supported here -- sorry.\n" );
		return 1;
	#endif

		omp_set_num_threads( NUMT );
		fprintf( stderr, "Using %d threads\n", NUMT );
		double totalExecutionTime = 0;
		int i, t, k = 0;
		double maxMegaMults = 0.;
		double sumMegaMults = 0.;
		double time0;
		double time1;
		printf( "Execution times:\n", time0, time1 );
		
    for(  t = 0; t < NUMTRIES; t++ )
		{
			time0 = omp_get_wtime( );

			#pragma omp parallel for
			for( i = 0; i < ARRAYSIZE; i++ )
			{
				C[i] = A[i] * B[i];
			}

			time1 = omp_get_wtime( );
			printf( "	time0 %lf time1 %lf  \n", time0, time1 );
			totalExecutionTime += (time1 - time0);
			k++;
			double megaMults = (double)ARRAYSIZE/(time1-time0)/1000000.;
			sumMegaMults += megaMults;
			if( megaMults > maxMegaMults )
				maxMegaMults = megaMults;
		}

		double avgMegaMults = sumMegaMults/(double)NUMTRIES;
		printf( "Peak Performance = %8.2lf MegaMults/Sec\n", maxMegaMults );
		printf( "Average Performance = %8.2lf MegaMults/Sec\n", avgMegaMults );
		printf( "Average Execution Time = %lf Secs\n", totalExecutionTime/k);



		return 0;
}
