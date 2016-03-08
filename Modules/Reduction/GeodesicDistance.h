/*=========================================================================

  Module    : Solid Mechanics
  File      :
  Copyright : (C)opyright 2011++
              See COPYRIGHT statement in top level directory.
  Authors   : D. Millan
  Modified  :
  Purpose   : Dijkstra algorithm to compute the graph geodesic.
  Date      :
  Version   :
  Changes   :

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.
=========================================================================*/
#pragma once

#include "smUtilities.h"
#include "smArray.h"
#include "smMath.h"
#include "smNNAdjacencyStructure.h"
#include "smFibonacciHeap.h"
#include "smHeapNode.h"
#include "smSparseMatrix.h"
#include "smSort.h"

//Intel Math Kernel Library
#include <mkl.h>

//C++ Standard Libraries
#include <vector>
#include <memory>
#include <exception>
#include <stdexcept>
#include <algorithm>
#include <functional>
#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>
#include <numeric>
using namespace std;

/** \class smGeodesicDistance
 * \brief This class uses the Dijkstra algorithm to compute the geodesic distance through a graph.
 *
 * In the mathematical field of graph theory, the distance between two vertices in a graph
 * is the number of edges or lenght in a shortest path connecting them. This is also known as the
 * geodesic distance because it is the length of the graph geodesic between those two vertices.
 * If there is no path connecting the two vertices, i.e., if they belong to different
 * connected components, then conventionally the distance is defined as infinite.
 *
 * See <a href="http://en.wikipedia.org/wiki/Distance_(graph_theory)">Geodesic graph distance</a>
 * 
 * Based on ISOMAP code which can be found at <a href="http://isomap.stanford.edu">Isomap</a>
 * 
 * \ingroup Utilities
*/
class GeodesicDistance
{
	public :
		
		// ----------------------------------------------------
		// Typedefs
		//-----------------------------------------------------
		typedef smGeodesicDistance Self;
		typedef smSmartPointer<Self> Pointer;
        typedef smEuclidianDistance Euclidian;
        typedef smNNAdjacencyStructure< Euclidian > NNAdjacencyStructure;
        typedef NNAdjacencyStructure::Pointer       NNAdjacencyPointer;
        
		// ----------------------------------------------------
		// Methods
		smNewMacro(Self); //Build-> :: New()

        /** This is the input for this class. This function sets the reference nearest
         * neigbord adjacency structure for each node point */
        void SetAdjacencyStructure( const NNAdjacencyPointer nears )
        {   this->nears = nears;    }
        
		/** Get the error for the Dijkstra algorithm \n
		 * 0 : means all OK \n
		 * 1 : input settings are wrong,
         * 2 : The number of non-zeros of the geodesic distance matrix is bad defined (nnz<1)
		 */
		int GetError( ) const
		{	return this->ierr;      }
		
		/** Bring this algorithm's outputs up-to-date. */
		virtual void Update();

        //@{
        /** Maximum number of nearest neighbors stored in the Geodesic Distance Matrix.*/
        void SetNumberOfNearestNeighbours ( int knn )
        {   this->kNN   = knn;      }
        int  GetNumberOfNearestNeighbors ( )
        {   return this->kNN;       }
        //@}
        
        //@{
        /** By DEFAULT the output distance matrix have the identity entry which always have
         * a zero value dist(i,i)=0.\n
         * It means that each row contains the self index identifier i.e. row I have non zero
         * values in many columns and a zero entry in the column I.\n
         * In order to avoid the self index entry this function is provided.*/
        void SetSelfIdentifierON ( )
        {   this->selfId = true;      }
        void SetSelfIdentifierOFF( )
        {   this->selfId = false;     }
        //@}

        /** Only the UPPER part of the distance matrix in CSR format is stored (default ON)
         * This option only works when the full Geodesic Distance matrix is computed, kNN=nPts-1.
         * Also in this case the data stored are averaged such that this is a symmetric matrix.
         *
         * Note that for a distance matrix from a graph the condition of symmetry is not
         * always fulfilled.
         *
         * Also with this option less memory is used, then larger set of data can be processed.  */
        //@{
        void SetUpperStorageON()
        {   this->upper = true;     }
        void SetUpperStorageOFF()
        {   this->upper = false;    }
        //@}

        /**Main output of this class.\n
         * The Geodesic Distance Matrix in compressed sparse by row (CSR) format.
         * For each point/vertex/row the geodesic graph distance for the first \b kNN
         * nearest neighbors is computed by using the Dijkstra's algorithm.
         */
        SparseMatrix *GetGeodesicDistanceMatrix()
        {   return &Gmat;           }

	protected:
        bool    selfId;     ///< each row contains the proper index identifier (ON), dist(i,i)=0
        bool    upper;      ///< only the upper part is stored (Full matrix case)
        int     ierr;       ///< info error identifier
        int     nPts;       ///< number of points/vertes in the graph, rows of the matrix distance
        int     kNN;        ///< maximum number of nearest neighbors to compute the Geodesic distance
        cint    *in_nodes;  ///< cumulative number of points by row is stored (internal array of nears)
        cint    *jn_nodes;  ///< adjacency structure with the point indices of the sparse matrix of Euclidian distance (internal array of nears)
        cdouble *dn_nodes;  ///< sparse matrix with the Euclidian distance (internal array of nears)

        SparseMatrix Gmat; ///< main output, geodesic distances from the sources to its k-nn

        NNAdjacencyPointer nears;   ///< nearest neighbors adjacency object.

		/** Method to delete the preallocated member variables */
		virtual void Clear();
		
		/** Check if the parameters and the data of the nodes are set.*/
		int CheckSettings();
		
	    /** Constructor */
		smGeodesicDistance(void);
 
		/** Destructor */
		virtual ~smGeodesicDistance();

	private:
        int   kNN_;     ///<internal variable, it contains the maximum number of nearest neighbors to compute the Geodesic distance

        void FullDijkstra( int s, int *Jn, double *Dn );
        void ShortestFullPath();

        void PartialDijkstra( int s, int *rangeID, int *mark_neigh, int *Jn, double *Dn );
        void ShortestPartialPath();


		smGeodesicDistance(const Self&);    //purposely not implemented
		void operator=(const Self&);        //purposely not implemented
};
