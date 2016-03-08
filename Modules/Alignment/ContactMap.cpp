/*=========================================================================
 * File         : ContactMap.h
 * Module       : Alignment
 * Copyright    : (C)opyright 2011-2013
 *                See COPYRIGHT statement in top level directory.
 * Authors      : Behrooz Hashemian
 * Created      : Fri Nov 22, 2013  01:47PM
 * Last modified: Fri Nov 22, 2013  01:47PM
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
#include "ContactMap.h"

// initialize _rij, _d0
ContactMap::ContactMap(int sdim, int nAtoms, double d0)
    : _sdim(sdim),
      _nAtoms(nAtoms)
{
    _d0inv = 1.0/(d0*d0);
    _cmap .resize( _nAtoms*(_nAtoms-1)/2               , 0.0 );
    _dcmap.resize( _nAtoms*(_nAtoms-1)/2 *_nAtoms*_sdim, 0.0 );
    _rij  .resize( _sdim                               , 0.0 );
    cout<<"nAtoms = "<<_nAtoms<<endl;
    cout<<"sdim   = "<<_sdim<<endl;
    cout<<"d0inv  = "<<_d0inv<<endl;
}



ContactMap::ContactMap(const ContactMap& inCMap)
{
    _sdim   = inCMap._sdim;
    _nAtoms = inCMap._nAtoms;
    _d0inv  = inCMap._d0inv;
    _cmap   = inCMap._cmap;
    _dcmap  = inCMap._dcmap;
    _rij    = inCMap._rij;
}


vector<double> ContactMap::ComputeAlignment( const vector<double> &x_new_full )
{
    int iN,jN;
    int iC = 0;
    for (iN = 0; iN < _nAtoms; ++iN )
        for (jN = iN+1; jN < _nAtoms; ++jN )
        {
            transform( x_new_full.begin()+iN*_sdim, x_new_full.begin()+(iN+1)*_sdim,
                       x_new_full.begin()+jN*_sdim, _rij.begin(), std::minus<double>() );
            _d = 0.0;
            for ( auto &item : _rij )
                _d += item*item;
            _d *= _d0inv;
            _cmap.at(iC) = (1.0-pow(_d,3))/(1.0-pow(_d,5));
            ++iC;
        }
    return _cmap;
}

vector<double> ContactMap::ComputeDerivatives(const vector<double> &x_new_full)
{
    int iN,jN,iD;
    int iC = 0;
    double dtemp;
    for (iN = 0; iN < _nAtoms; ++iN )
        for (jN = iN+1; jN < _nAtoms; ++jN )
        {
            transform( x_new_full.begin()+iN*_sdim, x_new_full.begin()+(iN+1)*_sdim,
                       x_new_full.begin()+jN*_sdim, _rij.begin(), std::minus<double>() );
            _d = 0.0;
            for ( auto &item : _rij )
                _d += item*item;
            _d *= _d0inv;
            _cmap.at(iC) = (1.0-pow(_d,3))/(1.0-pow(_d,5));
            for(iD = 0; iD < _sdim; ++iD)
            {
                dtemp = 2.0*_rij.at(iD)*_d0inv * _d*_d * (5.0*_d*_d*_cmap.at(iC)-3.0) / (1.0-pow(_d,5)) ;
                _dcmap.at(iC*_nAtoms*_sdim+iN*_sdim+iD) = +dtemp;
                _dcmap.at(iC*_nAtoms*_sdim+jN*_sdim+iD) = -dtemp;
            }
            ++iC;
        }
    return _dcmap;
}
