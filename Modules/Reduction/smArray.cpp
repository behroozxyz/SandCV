/*=========================================================================

  Module    : Solid Mechanics
  File      :
  Copyright : (C)opyright 2007++
              See COPYRIGHT statement in top level directory.
  Authors   : D. Millan, A. Rosolen
  Modified  : 
  Purpose   :
  Date      :
  Version   :
  Changes   :

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.
=========================================================================*/
#include "smArray.h"

// Constructor
smArray::smArray()
{
}

//Destructor
smArray::~smArray()
{
}

// Checks if the arrays w1 in wa2 are the same
template <class T> 
	inline bool smEqual_OMP ( int N, const T *wa1 , T *wa2 )
{
	register int i, j;
	int size, P;
	bool *equal[NUM_THREADS];
	
	P = omp_get_num_threads();
	if ( P == 1)
		return smArray::Equal(N, wa1, wa2);
	
	if ( sizeof(T) == 4 )
		size = LINE_SIZE_I;
	else if ( sizeof(T) == 4 )
		size = LINE_SIZE_F;
	else if ( sizeof(T) == 8 )
		size = LINE_SIZE_D;
	else
	{
		printf("ERROR::smArray::Equal_OMP() Only int, float and double types are allowed\n");
		exit(1);
	}
	
	for (i=0; i<NUM_THREADS; i++)
		equal[i] = new bool[size];
	for (i=0; i<NUM_THREADS; i++)
		equal[i][0] = true;
	
	// In absence of "chunk", each thread executes approx. N/P chunks 
	// 	for a loop of length N and P threads
#ifdef _OPENMP
#pragma omp parallel default(none) private(i,j) shared(N, wa1, wa2, equal)
#endif
	{
		j = omp_get_thread_num();
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
		for( i=0; i<N; i++ )
		{
			if ( wa1[i] != wa2[i] )
			{
				equal[j][0] = false;
				i = N;
			}
		}
	} //-- End of parallel for
	
	// checking the partition results
	for ( i=0; i<NUM_THREADS; i++ )
	{
		if ( equal[i][0] == false )
			return false;
	}
	
	// checking the boundary of each partition 
	size = N/P;
	for ( i=size-1; i<(N-1); i+=size )
	{
		if ( wa1[i] != wa2[i] )
			return false;
	}
	
	return true;
}
bool smArray::Equal_OMP ( int N, const    int *wa1,    int *wa2 ) { return smEqual_OMP(N, wa1, wa2); }
bool smArray::Equal_OMP ( int N, const  float *wa1,  float *wa2 ) { return smEqual_OMP(N, wa1, wa2); }
bool smArray::Equal_OMP ( int N, const double *wa1, double *wa2 ) { return smEqual_OMP(N, wa1, wa2); }


//  Cumulative sum of the work array wa
template <class T>
		inline bool smIsSorted_OMP ( int N, const T *wa )
{
	register int i, j;
	int size, P;
	bool *isSort[NUM_THREADS];
	
	P    = omp_get_num_threads();
	if ( P == 1)
		return smArray::IsSorted(N, wa);
	
	if ( sizeof(T) == 4 )
		size = LINE_SIZE_I;
	else if ( sizeof(T) == 4 )
		size = LINE_SIZE_F;
	else if ( sizeof(T) == 8 )
		size = LINE_SIZE_D;
	else
	{
		printf("ERROR::smArray::IsSorted_OMP() Only int, float and double types are allowed\n");
		exit(1);
	}
	
	for (i=0; i<NUM_THREADS; i++)
		isSort[i] = new bool[size];
	for (i=0; i<NUM_THREADS; i++)
		isSort[i][0] = true;
	
	// In absence of "chunk", each thread executes approx. N/P chunks 
	// 	for a loop of length N and P threads
#ifdef _OPENMP
#pragma omp parallel default(none) private(i,j) shared(N, wa, isSort)
#endif
	{
		j = omp_get_thread_num();
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
		for( i=0; i<(N-1); i++ )
		{
			if ( wa[i+1] < wa[i] )
			{
				isSort[j][0] = false;
				i = N;
			}
		}
	} //-- End of parallel for
	
	// checking the partition results
	for ( i=0; i<NUM_THREADS; i++ )
	{
		if ( isSort[i][0] == false )
			return false;
	}
	
	// checking the boundary of each partition 
	size = N/P;
	for ( i=size-1; i<(N-1); i+=size )
	{
		if ( wa[i+1] < wa[i] )
			return false;
	}
	
	return true;
}
bool  smArray::IsSorted_OMP ( int N,    cint *wa ) { return smIsSorted_OMP(N, wa); }
bool  smArray::IsSorted_OMP ( int N,  cfloat *wa ) { return smIsSorted_OMP(N, wa); }
bool  smArray::IsSorted_OMP ( int N, cdouble *wa ) { return smIsSorted_OMP(N, wa); }

