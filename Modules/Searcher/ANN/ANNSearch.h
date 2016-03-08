/*=========================================================================
 * File         : ANNSearch.h
 * Module       : Searcher
 * Copyright    : (C)opyright 2013
 *                See COPYRIGHT statement in top level directory.
 * Authors      : Behrooz Hashemian
 * Created      : Wed Jun 24, 2013  07:40PM
 * Last modified: Mon Jul 22, 2013  09:28AM
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS ``AS IS''
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *=========================================================================*/
#pragma once

//Aproximate Nearest Neighbor Search
#include "ANN.h"
//C++ Standard Libraries
#include <memory>
#include <vector>
using namespace std;

/** \class ANNSearch
 * \brief This class computes and stores the Approximate Nearest Neighbor (ANN)
 * searching three for a given input data set. Such that find the indices and
 * distances of the k-nearest neighbors of a query point.
 *
 * \param [in] dim          Spatial dimension
 * \param [in] knn          Number of requested nearest neighbor
 * \param [in] nDataPoints  Number of node points
 * \param [in] dataVector   A vector of data points
 * \param [in] dataArray    An array of data points
 *
 * \param [out] nnIndex     Indices of k-nearest neighbors
 * \param [out] nnDists     Euclidean distances of k-nearest neighbors 
*/
class ANNSearch
{
public :
    /** Default Constructor */
    ANNSearch() = default;

    /** Constructor with data point initialization */
    ANNSearch( double radius, int dim, const vector<double>& dataVector );
    ANNSearch( int knn, int dim, const vector<double>& dataVector );
    ANNSearch( int knn, int dim, int nDataPints, const double* dataArray );

    /** Destructor */
    ~ANNSearch();

    /** Copy constructor (explicitly deleted) */
    ANNSearch( const ANNSearch& ) = delete;

    /** Assignmet operator (explicitly deleted) */
    void operator=( const ANNSearch& ) = delete;

    /**  Get number of node points */
    int GetNumberOfNodePoints(  ) const { return _nDataPoints; }

    /** The node points coordinates are set */
    const vector<double> GetNodePoints( ) const { return _dataVector; }

    /** K-nearest neighbors searcher */
    const vector<int>& FindNearestNeighbors( const vector<double> &queryPoint );
    const int*         FindNearestNeighbors( double *queryInput );

    /** Fixed-radious nearest neighbors searcher */
    const vector<int>& FindNearestNeighborsFR( const vector<double> &queryPoint );

    /** Returns the indices for the k-nearest neighbors to x_sample (Euclidian metric)*/
    const vector<int>& GetNearestIndices() { return _nnIndex; }

    /** Returns the Euclidian distances between x_sample and its k-nearest neighbors */
    const vector<double>& GetNearestDistances() { return _nnDists; }

protected:
    /** @name Data Inputs: */
    int _nDataPoints = 0;   ///< number of node points
    int _dim         = 0;   ///< spatial dimension
    int _knn         = 0;   ///< number of requested Nearest Neighbours
    double _sqRadius = 0.0; ///< squared radius, for fixed-radius search

    vector<double> _dataVector; ///< spatial coordinates for the node points ( default=NULL )

    /** @name ANN objects */
    // ANNpointArray
    vector<double*>        _dataPoints;  ///< array of data points
    vector<double>         _queryPoint;  ///< query point (test sample point)
    vector<int>            _nnIndex;     ///< nearest neighbors indices
    vector<double>         _nnDists;     ///< nearest neighbors distances
    shared_ptr<ANNkd_tree> _tree;        ///< search structure

    /** Update the variables and make the class ready to use. */
    void Update() throw( exception );

    /** Check if parameters and point coordinates are properly set. */
    bool CheckSettings();
};
