/*=========================================================================
 * File         : ANNSearch.cpp
 * Module       : Searcher
 * Copyright    : (C)opyright 2013
 *                See COPYRIGHT statement in top level directory.
 * Authors      : Behrooz Hashemian
 * Created      : Wed Jun 24, 2013  07:40PM
 * Last modified: Thu Jul 25, 2013  05:40PM
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
#include "ANNSearch.h"
#include <algorithm>
#include <utility>
using namespace std;

ANNSearch::ANNSearch( double radius, int dim, const vector<double>& dataVector )
    : _dim(dim),
      _dataVector(dataVector)
{
    _sqRadius = radius*radius;
    _nDataPoints = (double)_dataVector.size() / (double)_dim;
    _nnIndex.reserve( _nDataPoints ); // Reserve the highest amount needed memory for near neighbor indices
    _nnDists.reserve( _nDataPoints ); // Reserve the highest amount needed memory for near neighbor distances
    this->Update();
    cout<<"Fixed radius Search..."<<endl;
}

ANNSearch::ANNSearch( int knn, int dim, const vector<double>& dataVector )
    : _dim(dim),
      _knn(knn),
      _dataVector(dataVector)
{
    _nDataPoints = _dataVector.size() / _dim;
    _nnIndex.resize( _knn );  // Allocate memory for near neighbor indices
    _nnDists.resize( _knn );  // Allocate memory for near neighbor distances
    this->Update();
    cout<<"KNN Search..."<<endl;
}


ANNSearch::ANNSearch( int knn, int dim, int nDataPoints, const double* dataArray )
{
    _knn = knn;
    _dim = dim;
    _nDataPoints = nDataPoints;
    _dataVector.assign(dataArray,dataArray+_nDataPoints*_dim);
    _nnIndex.resize( _knn );  // allocate near neighbor indices
    _nnDists.resize( _knn );  // allocate near neighbor distances
    this->Update();
}


ANNSearch::~ANNSearch()
{
    for ( int i = 0; i < _nDataPoints; ++i )
        delete[] _dataPoints.at(i);
    annClose();
}


void ANNSearch::Update()
throw( exception )
{
    if ( this->CheckSettings() )
    {
        cerr<<__func__<<" >> CheckSettings went wrong."<<endl;
        throw exception();
    }
    int i, j;
    // allocate data points using smart pointers
    _dataPoints.resize( _nDataPoints );
    for ( i = 0; i < _nDataPoints; ++i )
        _dataPoints.at(i) = new double [_dim];


    // Convert the input data points vector to ANN data array (_dataPts)
    for( i = 0; i < _nDataPoints; ++i )
        for( j = 0; j < _dim; ++j )
            _dataPoints.at(i)[j] = _dataVector.at( i*_dim + j );

    // build search structure (ANN tree)
    _tree = make_shared<ANNkd_tree>( _dataPoints.data(), _nDataPoints, _dim );

    cout<<"The tree has been updated with "<<endl<<
          "\tdim   = "<<_tree->theDim() <<endl<<
          "\tnData = "<<_tree->nPoints()<<endl<<
          "\tknn   = "<<_knn<<endl<<
          "\tradius= "<<sqrt(_sqRadius)<<endl;
}


bool ANNSearch::CheckSettings()
{
    bool bError = false;
    if (_dim < 1)
    {
        cerr<<__func__<<" >>  Spatial dimension is "<<_dim<<" < 2 "<<endl;
        bError = true;
    }
    if (_nDataPoints < 1)
    {
        cerr<<__func__<<" >>  Number of data points is "<<_nDataPoints<<" < 2 "<<endl;
        bError = true;
    }
    if (_dataVector.empty())
    {
        cerr<<__func__<<" >>  Data point vector is empty"<<endl;
        bError = true;
    }
    if (_sqRadius > 0.0)
    {
        cout<<"ANN range search is set."<<endl;
    }
    else if (_knn < 1)
    {
        cerr<<__func__<<" >>  Number of requested near neighbors is "<<_knn<<" < 1 "<<endl;
        bError = true;
    }
    else if ( _knn > _nDataPoints)
    {
        cerr<<__func__<<" >>  Number of requested nearest neighbors ("<<_knn
              <<") is bigger than number of data points("<<_nDataPoints<<") "<<endl;
        bError = true;
    }
    return bError;
}


const vector<int> &ANNSearch::FindNearestNeighbors(const vector<double>& queryPoint)
{
    // Find its k-nearest neighbors with ANN
    _tree->annkSearch( const_cast<double *>(queryPoint.data()), // query point
                       _knn,               // number of nearest neighbors
                       _nnIndex.data(),    // nearest neighbors (returned)
                       _nnDists.data());   // squared distance  (returned)
    // ANN use the squared of the distance as distance measure
    //  so we should take the square root of the output distance
    //     transform(_nnDists.begin(), _nnDists.end(),
    //              _nnDists.begin(), (double(*)(double)) sqrt);
    return _nnIndex;
}


const int* ANNSearch::FindNearestNeighbors(double* queryPoint)
{
    _tree->annkSearch( queryPoint,       // query point
                       _knn,             // number of nearest neighbors
                       _nnIndex.data(),  // nearest neighbors (returned)
                       _nnDists.data()); // squared distance (returned)
    return _nnIndex.data();
}


const vector<int> &ANNSearch::FindNearestNeighborsFR(const vector<double> &queryPoint)
{
    _knn = _tree->annkFRSearch( const_cast<double *>(queryPoint.data()), // query point
                                _sqRadius, // squared radious
                                0, NULL, NULL);
//    cout<<"-knn- = "<<_knn<<endl;
    _nnIndex.resize(_knn);
    _nnDists.resize(_knn);
    // Find all the neighbor points in the given radius
    _tree->annkSearch( const_cast<double *>(queryPoint.data()), // query point
                       _knn,             // number of nearest neighbors
                       _nnIndex.data(),  // nearest neighbors (returned)
                       _nnDists.data()); // squared distance  (returned)
    return _nnIndex;
}
