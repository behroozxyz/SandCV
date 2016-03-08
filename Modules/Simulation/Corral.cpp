/*=========================================================================
 * File         : Corral.cpp
 * Module       : Simulation
 * Copyright    : (C)opyright 2013
 *                See COPYRIGHT statement in top level directory.
 * Authors      : Behrooz Hashemian
 * Created      : Sun Aug 11, 2013  02:15PM
 * Last modified: Sun Aug 11, 2013  02:15PM
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
#include "Corral.h"

Corral::Corral(const vector<double> &xiCorral, int dim, double spacing , double corralHeight)
    : _xiMarkers(xiCorral),
      _dim(dim),
      _spacing(spacing),
      _corralHeight(corralHeight)
{
    _nMarkers = _xiMarkers.size() / _dim;
    _xiMarkers   = xiCorral;
    this->GenerateCorralMarkers();
    this->Update();
}


Corral::Corral(string filename , double corralHeight)
    :  _corralHeight(corralHeight)
{
    _reader.ReadPointSet(filename, _xiMarkers, _nMarkers, _dim);
    _spacing        = _xiMarkers.at(0);
    _nInsideMarkers = _xiMarkers.at(_dim);
    _xiMarkers.erase(_xiMarkers.begin(), _xiMarkers.begin()+2*_dim);
    _nMarkers      -= 2;
    this->Update();
}


void Corral::Update( )
throw( invalid_argument )
{
    if ( this->CheckSettings() )
    {
        cerr<<__func__<<" >> CheckSettings went wrong."<<endl;
        throw invalid_argument("One or more inputs are wrong.");
    }
    _betaPU  = _gammaPU/(_spacing*_spacing);
    _rangePU = sqrt(-log(_tolPU)/_betaPU);
    cout<<"Corral is initialized."<<endl<<
          "\tDimension                = "<<_dim<<endl<<
          "\tspacing                  = "<<_spacing<<endl<<
          "\tbeta                     = "<<_betaPU<<endl<<
          "\trange                    = "<<_rangePU<<endl<<
          "\tNumber of    all markers = "<<_nMarkers<<endl<<
          "\tNumber of inside markers = "<<_nInsideMarkers<<endl;
    _markerSearch = make_shared<ANNSearch>(_rangePU, _dim, _xiMarkers);
    _dZ  .resize(_dim);
    _diff.resize(_dim);

}


bool Corral::CheckSettings( )
{
    bool bError = false;

    if (_dim < 1)
    {
        cerr<<__func__<<" >>  Corral dimension is less than 1("<<_dim<<")."<<endl;
        bError = true;
    }

    if (_spacing <= 0.0 )
    {
        cerr<<__func__<<" >>  Corral spacing is wrongly set (spacing ="<<_spacing<<")."<<endl;
        bError = true;
    }

    return bError;
}


void Corral::CalculateCorral( const vector<double> &xiNew,
                              vector<double> &gradient)
{
    _wIn = 0.0;
    fill( gradient.begin(), gradient.end(), 0.0 );
    vector<int> markID = _markerSearch->FindNearestNeighborsFR(xiNew);

    if ( markID.size() < 1 )
    {
        cerr<<"There is no PU marker near this point : \n xi = ";
        for (auto iter : xiNew)
            cerr<<iter<<" ";
        cerr<<endl;
    }
    else
    {
        _Z     = 0.0;
        fill( _dZ.begin(), _dZ.end(), 0.0 );
        for ( const auto &mid : markID )
        {
            transform(xiNew.begin(), xiNew.end(), _xiMarkers.begin()+_dim*mid,
                      _diff.begin(), std::minus<double>() );
            _w  = exp( -_betaPU * cblas_ddot( _dim, _diff.data(), 1, _diff.data(), 1 ) );
            _Z += _w;
            cblas_daxpy( _dim, -2.0*_betaPU*_w, _diff.data(), 1, _dZ.data(), 1 );

            if (mid >= _nInsideMarkers)
            {
                _wIn += _w;
                cblas_daxpy( _dim, -2.0 *_betaPU*_w, _diff.data(), 1, gradient.data(), 1 );
            }
        }

        for ( int i = 0; i < _dim; ++i)
            gradient.at(i) = _corralHeight * ( _wIn *_dZ.at(i)/_Z - gradient.at(i) ) / _Z;
    }
}


void Corral::GenerateCorralMarkers()
{
    int iD;
    vector<double> nodeMin(_dim, +numeric_limits<double>::max()); // minumum of node points in each direction
    vector<double> nodeMax(_dim, -numeric_limits<double>::max()); // maximum of node points in each direction
    vector<double> nodeRange(_dim, 0.0);                          // range   of node points in each direction
    //______________________________________________________________________________
    // Find maximum, minimum and range of nodes in each dimension
    for ( auto iter = _xiMarkers.begin(); iter != _xiMarkers.end(); )
        for ( iD = 0; iD < _dim; ++iD, ++iter )
        {
            nodeMin.at(iD) = min( *iter, nodeMin.at(iD) );
            nodeMax.at(iD) = max( *iter, nodeMax.at(iD) );
        }
    transform( nodeMax.begin(), nodeMax.end(), nodeMin.begin(), nodeRange.begin(), minus<double>() );

    cout<<"min nodes = [";
    for( auto iter : nodeMin )
        cout<<iter<<"  ";
    cout<<"]"<<endl;
    cout<<"max nodes = [";
    for( auto iter : nodeMax )
        cout<<iter<<"  ";
    cout<<"]"<<endl;
    //______________________________________________________________________________
    // TODO: Implement corral marker generation
    cout<<"Warning: Corral marker genration has not been implemented."<<endl;
}


