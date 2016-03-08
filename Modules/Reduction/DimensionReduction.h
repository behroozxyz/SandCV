/*=========================================================================

  Module    : Solid Mechanics
  File      : smDimensionReduction.h
  Copyright : (C)opyright 2009++
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

#ifndef __smDimensionReduction_h
#define __smDimensionReduction_h

// Standard C header files
#include <cmath>
#include <stdio.h>
#include <stdlib.h>

// Standard C++ header files
#include <iterator>

// SolMec lib
#include "smArray.h"
#include "smAtriaSearch.h"
#include "smMath.h"
#include "smLinearAlgebra.h"
#include "smMetric.h"
#include "smNNAdjacencyStructure.h"
#include "smNumerics.h"
#include "smPointSet.h"

/** \class smDimensionReduction
 * \brief This class should be provide the necessary methods for performing Linear and
 * Nonlinear Dimention Reduction computations.
 * The intend methods are: \n
 * \b Linear \b Methods:
 *  - PCA: Principal Component Analysis [1]
 *  - MDS: Multidimensional Scaling [2]
 *  - ICA: Independent Component Analysis [3]
 * 
 * \b Nonlinear \b Methods:
 *  - ISOMAP: isometric feature mapping [4] <a href="http://waldron.stanford.edu/~isomap/">link</a>
 *  - LLE: local linear embedding [5,6] with [6] <a href="http://www.cs.toronto.edu/~roweis/lle/">link</a>
 * 		-# HLLE: Hessian Eigenmaps [8] <a href="http://www-stat.stanford.edu/~donoho/index.html">link</a>
 * 		-# RLLE: Robust Locally Linear Embedding [10] <a href="http://www.cse.ust.hk/~hongch/research.htm">link</a>
 * 		-# MLLE: Modified Local Linear Embedding [11] <a href="">link</a>
 * 		-# RSLLE: Robust and Stable Locally Linear Embedding [9,12] <a href="">link</a>
 * 
 * \note This is an abstract class. The abstract classes are never instanciated.
 * For this reason, the New() operation will be never called and smNewMacro is meanless.
 * 
 * \b References: \n
 * [1] I. T. Jolliffe. Principal Component Analysis. Springer-Berlag, New York, 1986.
 * 
 * [2] T. Cox, M. Cox, Multidimensional Scaling (Chapman & Hall, London, 1994).
 * 
 * [3] A. Hyvarinen. "Survey on independent component analysis". Neural Computing Surveys,
 * Vol 2, 94-128, 1999.
 * 
 * [4] J Tenenbaum, V. De Silva and J. Langford. "A global geometric framework for
 * nonlinear dimension reduction". Science, Vol 290, 2319–2323, 2000.
 * <a href="http://www-clmc.usc.edu/publications/T/tenenbaum-Science2000.pdf">pdf</a>
 * 
 * [5] S. Roweis and L. Saul. "Nonlinear dimensionality reduction by locally linear
 * embedding. Science, Vol 290, 2323–2326, 2000.
 * <a href="http://www.sciencemag.org/cgi/reprint/290/5500/2323.pdf">pdf</a>
 * 
 * [6] D. De Ridder, R. Duin. "Locally Linear Embedding for Classification", Technical
 * Report PH-2002-1, Delf University of Technology, pattern Recognition Group.
 * <a href="http://ict.ewi.tudelft.nl/pub/dick/ph-2002-01.pdf">pdf</a>
 * 
 * [7] L. Saul and S. Roweis. "Think globally, fit locally: unsupervised learning of
 * nonlinear manifolds". Journal of Machine Learning Research, Vol 4, 119-155, 2003.
 * 
 * [8] D. Donoho and C. Grimes. "Hessian eigenmaps: Locally linear embedding techniques
 * for high-dimensional data". Proceedings of National Academy of Science, Vol 100, Nro 10,
 * 5591-5596, 2003.
 * 
 * [9] J. Park, Z. Zhang, H. Zha and R. Kasturi. "Local Smoothing for Manifold Learning".
 * IEEE Computer Society Conference on Computer Vision and Pattern Recognition (CVPR 04),
 * Vol 2, 452-459, 2004.
 * 
 * [10] H. Chang, D. Yeung. "Robust Locally Linear Embedding, Pattern Recognition", Vol 39,
 * Nro 6, 1053-1065, 2006.
 * 
 * [11] Z. Zhang, J. Wang. "Modified Local Linear Embedding Using Multiple Weights", Advanced in Neural
 * Information Porcessing Systems, 19, MIT Press, 2007.
 * <a href="http://books.nips.cc/papers/files/nips19/NIPS2006_0065.pdf">pdf</a>
 * 
 * [12] J. Wang. "Robust and Stable Locally Linear Embedding", IEEE 2008.
 * 
 * \todo ICA, MDS, RSLLE
 * \todo Add a class to compute intrinsic dimensionality of the data at this neighborhood
 * 
 * \ingroup Numerics
*/
class SOLMEC_EXPORT smDimensionReduction : public smNumerics
{
	public:
		// ----------------------------------------------------
		// Typedefs
		//-----------------------------------------------------
		typedef smDimensionReduction Self;
		typedef smSmartPointer<Self> Pointer;
		typedef smEuclidianDistance Euclidian;
		typedef smNNAdjacencyStructure< Euclidian > NNAdjacencyStructure;
		typedef NNAdjacencyStructure::Pointer		NNAdjacencyPointer;
		typedef smPointSet::Pointer 				PointSetPointer;

		/** Return the name of this class as a string. */
		const char *GetNameOfClass() const
		{	return "smDimensionReduction";		}
		
		//@{
		/** Spatial embedded dimension for the output projection, tDim \< nDim */
		void SetEmbeddedDimension( int tDim )
		{	this->tDim = tDim;				}
		int  GetEmbeddedDimension( ) const
		{	return this->tDim;				}
		//@}
		
		/** Number of nearest neighbors to use.
		 * The \e k neighbors span a space of dimensionality at most \e k-1.		*/
		void SetNumberOfNearestNeighbours ( int knn )
		{	this->kNN   = knn;		}


		//@{
		/** Set input: node point set data in the HIGH D-dimension. */
		virtual void SetPointSetNodes( PointSetPointer nodes )
		{	this->nodes = nodes;	}
		virtual PointSetPointer GetPointSetNodes()
		{	return this->nodes;		}
		//@}
		//@{
		/** Set input: samples point set data in the HIGH D-dimension */
		virtual void SetPointSetSamples( PointSetPointer samples )
		{	this->samples = samples;	}
		virtual PointSetPointer GetPointSetSamples()
		{	return this->samples;		}
		//@}

		
		//@{
		/** Spatial dimension for the input point set (training set)*/
		void SetDimension( int nDim)
		{	this->nDim = nDim;				}
		int  GetDimension( ) const
		{	return this->nDim;				}
		//@}
		//@{
		/**  Set/Get number of node points (training set)*/
		void SetNumberOfNodePoints( cint nPts )
		{	this->nPts = nPts;				}
		int  GetNumberOfNodePoints( ) const
		{	return this->nPts;				}
		//@}
		
		//@{
		/** The node points coordinates are set/get (training set)*/
		void SetNodePoints( cdouble *x_nodes )
		{	this->x_nodes  =  x_nodes;		}
		cdouble *GetNodePoints( ) const
		{	return this->x_nodes;			}
		//@}


        //@{
        /** The low-dimensional embedded node points coordinates are set/get (training set)*/
        void SetEmbeddedNodePoints( double *xi_nodes )
        {	this->xi_nodes  =  xi_nodes;	}
        double *GetEmbeddedNodePoints( ) const
        {	return this->xi_nodes;			}
        //@}

		//@{
		/**  Set/Get the number of sample points */
		void SetNumberOfSamplePoints( cint sPts )
		{	this->sPts = sPts;				}
		int  GetNumberOfSamplePoints( ) const
		{	return this->sPts;				}
        //@}

		//@{
		/** The sample points coordinates */
		void SetSamplePoints( cdouble *x_samples )
		{	this->x_samples =  x_samples;	}
		cdouble *GetSamplePoints( ) const
		{	return this->x_samples;			}
		//@}

		//@{
		/** This function sets the flag 'nodeProj' which indicates if the nodes will be
		 * projected onto the tangent space(default ON=true)*/
		void SetNodesProjectionON ()
		{	this->nodeProj = true;		}
		void SetNodesProjectionOFF()
		{	this->nodeProj = false;		}
		//@}

		//@{
		/** This function sets the flag 'sampProj' which indicates if the samples will be
		 * projected onto the tangent space(default OFF=false).
		 * 
		 * A new, previously unseen sample x can be mapped by finding its \e k nearest
		 * neighbours in \b X_n , calculating new reconstruction weights \e w_j and
		 * mapping it using \f$ y = \sum_j w_j y_j \f$.
		 * 
		 * \b WARNING this is not a parametric mapping; still, the axes of the mapped
		 * space can be seen as corresponding to latent variables in the original data.
		 */
		void SetSamplesProjectionON ()
		{	this->sampProj = true;		}
		void SetSamplesProjectionOFF()
		{	this->sampProj = false;		}
		//@}


		/** \b MAIN \b OUTPUT. The low-dimensional embedded coordinates for the samples from a dimension
		 * reduction are returned*/
		cdouble *GetSamplesEmbedding()
		{	return this->xi_samples;	}
		
		/** This function returns the error identifier
		 * \return ierr 0: all computations finished succesfully \n
		 * 				1: an input parameter either has not been set or has been set bad
		*/
		virtual int GetError() const
		{	return this->ierr;			}
		
		/** Bring the algorithm's outputs up-to-date. */
		virtual void Update() = 0;

	protected:
		/** @name  Data Inputs: */
		bool	nodeProj;   ///< flag to indicate that nodes will be projected (default ON=true)
		bool    sampProj;   ///< flag to indicate that samples will be defined and then projected to

		int     nDim;		///< spatial dimension ( = 0 )
		int     tDim;		///< spatial embedded dimension ( = 0 )
		int     nPts;		///< number of node points ( = 0 )
		int     sPts;		///< number of sample points ( = 0 )

		int     kNN;		///< number of nearest neighbors

		cdouble *x_nodes;   ///< spatial coordinates for the node points ( default=NULL )
		cdouble *x_samples; ///< spatial coordinates for the sample points ( default=NULL )
		
		PointSetPointer  nodes;   ///< poinset for the node training set D-dim
		PointSetPointer  samples; ///< poinset for the sample test set D-dim

		/** @name Data Outputs: */
		int     ierr;       ///< error identifier
		double  *xi_nodes;   ///< projected nodes (training set) in the low-dimensional space
		double  *xi_samples; ///< projected samples in the low-dimensional space

		/** @name Member internal variables and objects: */
		/** Here the graph structures for the nodes and samples are defined. Its are defined
		 * in this way in order to account different number of nearest neighbors at each
		 * node/sample (for a future no too far... we hope) */
		cint    *jn_knn;    ///< k-nn nodes adjacency list structure with the nearest neighbors
		cint    *in_knn;    ///< k-nn adjacency strucutre pointer indexes for the nodes [nPts+1]
		cint    *js_knn;    ///< k-nn samples adjacency list structure with the nearest neighbors
		cint    *is_knn;    ///< k-nn adjacency strucutre pointer indexes for the samples [sPts+1]
		
		NNAdjacencyPointer nears;	///< nearest neighbors adjacency object.

		/** @name Member methods: */
		/** Method to delete the preallocated member variables */
		virtual void Clear() = 0;

		/** Check if parameters and node/sample point coordinates are propely set. */
		virtual int CheckSettings() = 0;

		smDimensionReduction(void);
		virtual ~smDimensionReduction(void);
	private:
		smDimensionReduction(const Self&);  //purposely not implemented
		void operator=(const Self&);	    //purposely not implemented
};

#endif // _smDimensionReduction_h



