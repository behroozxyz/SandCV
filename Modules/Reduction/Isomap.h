/*=========================================================================

Module    : Solid Mechanics
File      : Isomap.h
Copyright : (C)opyright 201++
            See COPYRIGHT statement in top level directory.
Authors   : D. Millan
Modified  :
Purpose   :Isometric Mapping ISOMAP
Date      :
Version   :
Changes   :

    This software is distributed WITHOUT ANY WARRANTY; without even 
    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
    PURPOSE.  See the above copyright notices for more information.
=========================================================================*/

#pragma once

#include <time.h>
#include <sys/timeb.h>
#include "mkl_lapacke.h"

#include "DimensionReduction.h"
#include "SparseMatrix.h"
#include "GeodesicDistance.h"
#include "ANNSearch.h"

/** \class Isomap
* \brief The Isomap class performs Isometric embedding of a point set \b X_n in D-dimension.
* 
* 
* 
* \b References: \n
* [1] <a href="http://isomap.stanford.edu">Isomap</a>
* 
* \ingroup Numerics
*/
class Isomap
{
    public :
    
        /** Default Constructor */
        Isomap() = default;
        
        /** Constructors with initialization */
        Isomap(int knn, int lowDim, int highDim, const vector<double> &xNodes);
        
        /** Destructor */
        ~Isomap() = default;
        
        /** Copy constructor */
        Isomap(const Isomap& inISO);
        
        // Nearest neighbor search
        vector<int>  _idxNN;    ///< Nearest neighbor indices
        ANNSearch    _searcher; ///< Nearest neighbor search object
        
        
        /** This function returns the error identifier
        * \return ierr  0: all computations finished successfully \n
        *              11: Eigenvalues are NaN or Inf \n
        *              12: Eigenvalues are negatives \n
        *              13: Values from the eigenvectors are NaN or Inf
        */
        int GetError() const
        {	return this->ierr;			}
        
        /** Bring the algorithm's outputs up-to-date. */
        void Update();
        
    protected:
        SparseMatrix *M;	  ///< sparse matrix, actually we store in M the (I-W)'(I-W) values

        int    *jn_nodes;  ///< nearest nodes adjacency list structure with the nearest node neighbors
        int    *in_nodes;  ///< nearest adjacency structure pointer indexes for the nodes [nPts+1]
        double *dn_nodes;  ///< nodes distances in high dimension R^D for the nearest neighbors (sparse matrix)
        
        /** Check if parameters and node/sample point coordinates are properly set. */
        int CheckSettings();
        
        /** Method to delete the preallocated member variables */
        void Clear();

        /** Step 1: Geodesic Graph Distance is computed by using the Dijkstra's algorithm */
        void GeodesicDistance();

        /** Step 2: Double Centering of a distance matrix: Classical Multidimensional Scaling (MDS)
        * Center the Geodesic Graph matrix by subtracting the mean of each column, row,
        * adding the overall mean, and multiplying by -0.5
        *
        * DD = -0.5*[D.^2 - sum(D.^2)'*ones(1,N)/N - ones(N,1)*sum(D.^2)/N + sum(sum(D.^2))/(N^2)]
        */
        void DoubleCentering();
        
        /** Step 3: Solve Eigenproblem: dspevx (mkl_lapacke routine)
         *  Computes selected eigenvalues and eigenvectors of a real symmetric matrix in packed storage.
         */
        void EigenvalueProblem();
        
        /** Step 4: Compute node coordinates (y_nodes) in the lower-dimensional space */
        void Embedding();
        
    private:
        double *eigVal;   ///<array with the eigenvalues
        double *eigVec;   ///<array with the egenvectors, these are stored consecutively [v1 v2 ... vd]
        
        GeodesicDistance *geoDist;
};
