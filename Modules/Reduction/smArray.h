

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
#ifndef __smArray_h
#define __smArray_h

// Standard C++ header files
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>

#include "smUtilities.h"



/** \class smArray 
 * \brief smArray is provides methods to perform common operations on arrays. These 
 * include providing operations such as copying, initializing, checking, and reverse algorithm
 * 
 * \ingroup Utilities
 */
class SOLMEC_EXPORT smArray : public smUtilities
{
	public :
	
		// ----------------------------------------------------
		// Typedefs
		//-----------------------------------------------------
		typedef smArray Self;
		typedef smSmartPointer<Self> Pointer;	
		
		// ----------------------------------------------------
		// Methods
		smNewMacro(Self); //Build-> :: New()
		
		/** Return the name of this class as a string. */
		const char *GetNameOfClass() const
		{return "smArray";}
		
		/** The Copy() function copies the elements from wa1 to wa2. */
		//@{
		template <class T>
				inline static void Copy3D ( const T *wa1, T *wa2 )
		{
			wa2[0] = wa1[0];
			wa2[1] = wa1[1];
			wa2[2] = wa1[2];
			return;
		}
		template <class T>
				inline static void Copy ( int N, const T *wa1, T *wa2 )
		{
			register int i;
			for ( i=0; i < N; i++ )
				wa2[i] = wa1[i];
			return;
		}
        template <class T>
                inline static void Copy ( int N, int M, const T *wa1, T *wa2 )
        {
            register int i, j, row;
            for ( i=0; i < N; i++ )
            {
                row = i*M;
                for ( j=0; j<M; j++ )
                    wa2[row+j] = wa1[row+j];
            }
            return;
        }
		template <class T>
				inline static void Copy ( int N, std::vector<T> wa1, T *wa2 )
		{
			register int i;
			for ( i=0; i < N; i++ )
				wa2[i] = wa1[i];
			return;
		}
		template <class T>
				inline static void Copy_OMP ( int N, const T *wa1, T *wa2 )
		{
			register int i;
#ifdef _OPENMP
#pragma omp parallel for default(none) private(i) shared(N, wa1, wa2) schedule(static)
#endif
			for ( i=0; i < N; i++ )
				wa2[i] = wa1[i];
			return;
		}
		//@}
		/** This function returns true if the elements in two arrays are the same. */
		template <class T>
				inline static bool Equal ( int N, const T *wa1, const T *wa2 )
		{
			register int i;
			for ( i=0; i < N; i++ )
				if ( wa2[i] != wa1[i] )
					return false;
			return true;
		}
		static bool Equal_OMP ( int N,    cint *wa1,    int *wa2 );
		static bool Equal_OMP ( int N,  cfloat *wa1,  float *wa2 );
		static bool Equal_OMP ( int N, cdouble *wa1, double *wa2 );
		

		/** The function Fill() assigns val to all of the elements of wa. */
		//@{
		template <class T>
				inline static void Fill ( int N, T val, T *wa )
		{
			register int i;
			for ( i=0; i < N; i++ )
				wa[i] = val;
			return;
		}
		template <class T>
				inline static void Fill_OMP ( int N, T val, T *wa )
		{
			register int i;
#ifdef _OPENMP
#pragma omp parallel for default(none) private(i) shared(N, val, wa) schedule(static)
#endif
			for ( i=0; i < N; i++ )
				wa[i] = val;
			return;
		}
		template <class T>
				inline static void Fill3D ( T val, T *wa )
		{
			wa[0] = val;
			wa[1] = val;
			wa[2] = val;
			return;
		}
		//@}
		/** This function returns true if the elements in wa are sorted in ascending order
		 * and false in other ways. */
		template <class T>
				inline static bool IsSorted ( int N, const T *wa )
		{
			register int i;
			for ( i=0; i < N-1; i++ )
				if ( wa[i+1] < wa[i] )
					return false;
			return true;
		}
		static bool IsSorted_OMP ( int N,    cint *wa );
		static bool IsSorted_OMP ( int N,  cfloat *wa );
		static bool IsSorted_OMP ( int N, cdouble *wa );

		
		/** The Reverse() algorithm reverses the order of elements in the range. */
		template <class T>
				inline static void Reverse ( int N, T *wa )
		{
			register int i;
			T temp;

			for ( i = 0; i < N/2; i++ )
			{
				temp      = wa[i];
				wa[i]     = wa[N-1-i];
				wa[N-1-i] = temp;
			}
			return;
		}

		/** Array swaping, the positions a[i], and a[j] are swaped */
		template <class T>
				inline static void Swap( T* a, const long i, const long j )
		{ 
			T tmp = a[i];
			a[i]  = a[j];
			a[j]  = tmp;
		}

		/** Returns the largest element in a range [implemented from C++ Algorithms]*/
		template <class T>
				inline static T Max ( int N, T *wa )
		{
			return (*std::max_element(wa, wa+N));
		}
		/** Returns the smallest element in a range [implemented from C++ Algorithms]*/
		template <class T>
				inline static T Min ( int N, T *wa )
		{
			return (*std::min_element(wa, wa+N));
		}
		
		/**Program to find the position of a given element in a working array, \e wa, using Binary Search.
		 * If an element matching \e skey is found, the return value is its position in \e wa.
		 * Otherwise, the return value is the length of the array, \e n. */
		template <class T>
				inline static T BinarySearch(int n, const T *wa, T skey)
		{
			register int low, high, middle;
			low = 0;
			high= n-1;
			while(low<=high)
			{
				middle=(low+high)/2;
				if(skey==wa[middle])
				{
					//printf("Element found at %d",middle);
					return middle;
				}
				else
				{
					if(skey>wa[middle])
						low=middle+1;
					else
						high=middle-1;
				}
			}
			return n;
		}
	protected:
		
		/** Constructor */
		smArray(void);
 
		/** Destructor */
		virtual ~smArray();
		
	private:
		smArray(const Self&);			//purposely not implemented
		void operator=(const Self&);	//purposely not implemented
};

#endif // __smArray_h
