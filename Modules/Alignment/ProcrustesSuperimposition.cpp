/*=========================================================================
 * File         : ProcrustesSuperimposition.h
 * Module       : Alignment
 * Copyright    : (C)opyright 2011-2013
 *                See COPYRIGHT statement in top level directory.
 * Authors      : Behrooz Hashemian
 * Created      : Wed Jun 12, 2011  05:52PM
 * Last modified: Wed Jun 12, 2013  08:12PM
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
#include "ProcrustesSuperimposition.h"

ProcrustesSuperimposition::ProcrustesSuperimposition( int dim,
                                                      const vector<double> &x_ref_full )
{
    _sdim = dim;
    this->SetReferenceConfiguration( x_ref_full );
}


ProcrustesSuperimposition::ProcrustesSuperimposition( int dim,
                                                      const vector<double> &x_ref_full,
                                                      const vector<int> &involvedAtoms )
{
    _sdim = dim;
    this->SetReferenceConfiguration( x_ref_full, involvedAtoms );
}

ProcrustesSuperimposition::ProcrustesSuperimposition(const ProcrustesSuperimposition & inPS)
{

    _sdim = inPS._sdim;
    this->SetReferenceConfiguration( inPS._x_ref_full, inPS._involvedAtoms );
}


int ProcrustesSuperimposition::SetReferenceConfiguration( const vector<double> &x_points_reference )
{
    _x_ref_full    = x_points_reference;
    _nFullAtoms    = _x_ref_full.size() / _sdim;
    _nFullAtomsInv = 1.0 / (double)_nFullAtoms;
    _nAtoms        = _nFullAtoms;
    _nAtomsInv     = 1.0 / (double)_nAtoms;

    /* Construct the involved atoms for the sake of copy constructor */
    _involvedAtoms.resize(_nAtoms);
    iota(_involvedAtoms.begin(),_involvedAtoms.end(),0); //1:nAtoms

    return Update();
}


int ProcrustesSuperimposition::SetReferenceConfiguration(const vector<double> &x_points_reference , const vector<int> &involvedAtomes )
{
    _x_ref_full     = x_points_reference;
    _nFullAtoms     = _x_ref_full.size() / _sdim;
    _nFullAtomsInv  = 1.0 / (double)_nFullAtoms;
    _involvedAtoms  = involvedAtomes;
    _nAtoms         = _involvedAtoms.size();
    _nAtomsInv      = 1.0 / (double)_nAtoms;
    return Update();
}


int ProcrustesSuperimposition::Update( )
throw( exception )
{
    if ( this->CheckSettings() )
    {
        cerr<<__func__<<" >> CheckSettings went wrong."<<endl;
        throw exception();
    }
    this->UpdateActivePoints();
    _mean_x_ref.resize(_sdim);
    _mean_x_new.resize(_sdim);
    _hdim = _nFullAtoms*_sdim;
    // Alignment
    _ipiv      .resize(_sdim);
    _M         .resize(_sdim*_sdim);
    _U         .resize(_sdim*_sdim);
    _S         .resize(_sdim);
    _Vt        .resize(_sdim*_sdim);
    _Rot       .resize(_sdim*_sdim);
    _Trans     .resize(_hdim);
    _x_aligned .resize(_hdim);
    // Derivatives
    _Omega     .resize(_sdim*_sdim);
    _dRdr      .resize(_sdim*_sdim);
    _dCdr      .resize(_sdim*_sdim);
    _dAdr      .resize(_hdim*_hdim);

    return 0;
}


bool ProcrustesSuperimposition::CheckSettings( )
{
    bool bError = false;
    if ( _sdim < 2 )
    {
        cerr<<__func__<<" >>  Spatial dimension is "<<_sdim<<" < 2 "<<endl;
        bError = true;
    }
    if ( _nFullAtoms < 2 )
    {
        cerr<<__func__<<" >>  Number of nodes either is not defined or is less than 2 ( nAtoms = "<<_nFullAtoms<<" )."<<endl;
        bError = true;
    }
    return bError;
}


int ProcrustesSuperimposition::UpdateActivePoints( )
{

    if (_nAtoms == _nFullAtoms)
        _x_ref = _x_ref_full;
    else
    {
        int i,j;
        _x_ref.resize( _nAtoms*_sdim );
        for ( i = 0; i < _nAtoms; ++i )
            for ( j = 0; j < _sdim; ++j )
                _x_ref.at(i*_sdim + j) = _x_ref_full.at(_involvedAtoms.at(i)*_sdim + j);
        //vector<double>::iterator activeIterPos;
        //for ( int i = 0; i < _nAtoms; ++i )
        //{
        //    activeIterPos =  _x_ref_full.begin() + _involvedAtoms.at(i)*_dim;
        //    copy( activeIterPos, activeIterPos+_dim, _x_ref.begin()+i*_dim );
        //}
    }
    _mean_x_ref = this->MakeZeroMean( _x_ref );
    return 0;
}


int ProcrustesSuperimposition::UpdateActivePoints( vector<double> x_new_full )
{
    if (_nAtoms == _nFullAtoms)
        _x_new = x_new_full;
    else
    {
        int i,j;
        _x_new.resize( _nAtoms*_sdim );
        for ( i = 0; i < _nAtoms; ++i )
            for ( j = 0; j < _sdim; ++j )
                _x_new.at(i*_sdim + j) = x_new_full.at(_involvedAtoms.at(i)*_sdim + j);
        //vector<double>::iterator activeIterPos;
        //for ( int i = 0; i < _nAtoms; ++i )
        //{
        //    activeIterPos =  x_new_full.begin() + _involvedAtoms.at(i)*_dim;
        //    copy( activeIterPos, activeIterPos+_dim, _x_new.begin()+i*_dim );
        //}
    }
    _mean_x_new = this->MakeZeroMean( _x_new );

    return 0;
}


vector<double> ProcrustesSuperimposition::MakeZeroMean( vector<double> &x_points )
{
    int i,j;
    int    nAtoms    = x_points.size()/_sdim;
    double nAtomsInv = 1.0 / double(nAtoms);
    vector<double> meanValue(_sdim, 0.0);

    for ( i = 0; i < nAtoms; ++i )
        for ( j = 0; j < _sdim; ++j )
            meanValue.at(j) += x_points.at(i*_sdim+j);

    for ( auto &item : meanValue )
        item *= nAtomsInv;

    for ( i = 0; i < nAtoms; ++i )
        for ( j = 0; j < _sdim; ++j )
            x_points.at(i*_sdim+j) -= meanValue.at(j);

    return meanValue;
}


vector<double> ProcrustesSuperimposition::ComputeAlignment( const vector<double>& x_new_full )
{
    double detS;

    // Select a subset of x_new and put into _x_new
    // and then make it zero mean
    this->UpdateActivePoints( x_new_full );

    // Initialize translation matrix, tran = repmat(meanXref, nAtoms, 1)
    // meanXnewMatrix = repmat(meanXnew, nAtoms_full, 1) , use _x_aligned variable as a temporary storage
    for ( int iN = 0; iN < _nFullAtoms; ++iN )
    {
        copy( _mean_x_ref.begin(), _mean_x_ref.end(), _Trans.begin()    +iN*_sdim );
        copy( _mean_x_new.begin(), _mean_x_new.end(), _x_aligned.begin()+iN*_sdim );
    }

    // Covariance matrix: Cov = Xref' * Xnew (+ 0.0*A)
    cblas_dgemm ( CblasRowMajor, CblasTrans, CblasNoTrans, _sdim, _sdim, _nAtoms,
                  1.0, _x_ref.data(), _sdim, _x_new.data(), _sdim, 0.0 , _M.data(), _sdim );

    // Singular Value decomposition: Cov = U*S*Vt
    LAPACKE_dgesdd ( LAPACK_ROW_MAJOR , _job, _sdim, _sdim, _M.data(), _sdim,
                     _S.data(), _U.data(), _sdim, _Vt.data(), _sdim );

    // Calculate preliminary rotation matrix: Rot =  V * U' = Vt' * U';
    cblas_dgemm ( CblasRowMajor, CblasTrans, CblasTrans, _sdim, _sdim, _sdim,
                  1.0, _Vt.data(), _sdim, _U.data(), _sdim, 0.0, _Rot.data(), _sdim );

    // Determinant of rotation matrix (VU')
    detS = _Rot.at(1)*_Rot.at(5)*_Rot.at(6) + _Rot.at(2)*_Rot.at(3)*_Rot.at(7) + _Rot.at(0)*_Rot.at(4)*_Rot.at(8)
         - _Rot.at(2)*_Rot.at(4)*_Rot.at(6) - _Rot.at(0)*_Rot.at(5)*_Rot.at(7) - _Rot.at(1)*_Rot.at(3)*_Rot.at(8);

    // If the determinant is negative, multiply last row of Vt with -1
    if (detS < 0.0)
    {
        for( int iD = _sdim*(_sdim-1); iD <_sdim*_sdim; ++iD)
            _Vt.at(iD) *= -1.0;
        _S.at(_sdim-1) *= -1.0;

        // Calculate optimum rotation matrix: Rot =  V * U' = Vt' * U';
        cblas_dgemm ( CblasRowMajor, CblasTrans, CblasTrans, _sdim, _sdim, _sdim,
                      1.0, _Vt.data(), _sdim, _U.data(), _sdim, 0.0, _Rot.data(), _sdim );
    }

    // Optimum translation matrix: Trans = rep( meanXref - meanXnew*T , nAtoms, dim );
    // Here _x_aligned is a temporary storage of mean values of x_new
    cblas_dgemm ( CblasRowMajor, CblasNoTrans, CblasNoTrans, _nFullAtoms, _sdim, _sdim,
                  -1.0, _x_aligned.data(), _sdim, _Rot.data(), _sdim, 1.0, _Trans.data(), _sdim );

    // Rotate and translate the configuration, x_aligned = x_new * R + T
    _x_aligned = _Trans;
    cblas_dgemm ( CblasRowMajor, CblasNoTrans, CblasNoTrans, _nFullAtoms, _sdim, _sdim,
                  1.0 , x_new_full.data(), _sdim, _Rot.data(), _sdim, 1.0, _x_aligned.data(), _sdim );

    return _x_aligned;
}

vector<double> ProcrustesSuperimposition::ComputeDerivatives( vector<double> x_new_full )
{
    size_t i, j, k, m, n, p, q, s, t;

    this->ComputeAlignment( x_new_full );
    fill( _dAdr.begin(), _dAdr.end(), 0.0 );
    if ( _S.end() != unique( _S.begin(), _S.end() ) ) // FIXIT: This criterion is not correct.
    {
        cout<<"Derivatives are not Computed: There are at least two identical singular values."<<endl;
        _ierr = 1;
    }
    else
    {
        if ( _nAtoms == _nFullAtoms )
        {
            for ( i = 0; i < _nFullAtoms; ++i )
            {
                for ( j = 0; j < _sdim; ++j )
                {
                    fill(_dCdr.begin(), _dCdr.end(), 0.0);
                    for ( m = 0; m < _sdim; ++m )
                        for ( n = 0; n < _sdim; ++n )
                            _dCdr.at( m*_sdim + n ) = _x_ref.at( i*_sdim + m ) * delta(n,j);

                    fill(_Omega.begin(),_Omega.end(), 0.0);
                    for ( p = 0; p < _sdim; ++p )
                        for ( q = 0; q < _sdim; ++q )
                        {
                            for ( m = 0; m < _sdim; ++m )
                                for ( n = 0; n < _sdim; ++n )
                                    _Omega.at( p*_sdim + q ) +=
                                            ( _Vt.at(p*_sdim+m) * _dCdr.at(n*_sdim+m) * _U .at(n*_sdim+q) -
                                              _U .at(m*_sdim+p) * _dCdr.at(m*_sdim+n) * _Vt.at(q*_sdim+n) );
                            _Omega.at( p*_sdim + q ) /= _S.at(p) + _S.at(q);
                        }
                    fill(_dRdr.begin(),_dRdr.end(), 0.0);
                    for ( t = 0; t < _sdim; ++t )
                    {
                        for ( k = 0; k < _sdim; ++k )
                            for ( p = 0; p < _sdim; ++p )
                                for ( q = 0; q < _sdim; ++q )
                                    _dRdr.at( k*_sdim + t ) += _Vt.at(p*_sdim+k) * _Omega.at(p*_sdim+q) * _U.at(t*_sdim+q);

                        for ( s = 0; s < _nFullAtoms; ++s )
                        {
                            for ( k = 0; k < _sdim; ++k )
                                _dAdr.at( s*_sdim*_hdim + t*_hdim  + i*_sdim + j )
                                        += _x_new.at( s*_sdim + k ) *_dRdr.at( k*_sdim + t );
                            _dAdr.at( s*_sdim*_hdim + t*_hdim  + i*_sdim + j )
                                    += _Rot.at(j*_sdim+t) * ( this->delta(i,s) - _nFullAtomsInv );
                        }
                    }
                }
            }
        }
        else
        {
            for ( i = 0; i < _nFullAtoms; ++i )
                for ( j = 0; j < _sdim; ++j )
                    for ( s = 0; s < _nFullAtoms; ++s )
                        for ( t = 0; t < _sdim; ++t )
                            _dAdr.at( s*_sdim*_hdim + t*_hdim  + i*_sdim + j )
                                    = _Rot.at(j*_sdim+t) * this->delta(i,s);

            for ( i = 0; i < _nAtoms; ++i )
            {
                for ( j = 0; j < _sdim; ++j )
                {
                    for ( m = 0; m < _sdim; ++m )
                        for ( n = 0; n < _sdim; ++n )
                            _dCdr.at( m*_sdim + n ) = _x_ref.at( i*_sdim + m ) * delta(n,j);
                    fill(_Omega.begin(),_Omega.end(), 0.0);
                    for ( p = 0; p < _sdim; ++p )
                        for ( q = 0; q < _sdim; ++q )
                        {
                            for ( m = 0; m < _sdim; ++m )
                                for ( n = 0; n < _sdim; ++n )
                                    _Omega.at( p*_sdim + q ) +=
                                            ( _Vt.at(p*_sdim+m) * _dCdr.at(n*_sdim+m) * _U .at(n*_sdim+q) -
                                              _U .at(m*_sdim+p) * _dCdr.at(m*_sdim+n) * _Vt.at(q*_sdim+n) );
                            _Omega.at( p*_sdim + q ) /= _S.at(p) + _S.at(q);
                        }
                    fill(_dRdr.begin(),_dRdr.end(), 0.0);
                    for ( t = 0; t < _sdim; ++t )
                    {
                        for ( k = 0; k < _sdim; ++k )
                            for ( p = 0; p < _sdim; ++p )
                                for ( q = 0; q < _sdim; ++q )
                                    _dRdr.at( k*_sdim + t ) += _Vt.at(p*_sdim+k) * _Omega.at(p*_sdim+q) * _U.at(t*_sdim+q);

                        for ( s = 0; s < _nFullAtoms; ++s )
                        {
                            for ( k = 0; k < _sdim; ++k )
                                _dAdr.at( s*_sdim*_hdim + t*_hdim  + _involvedAtoms.at(i)*_sdim + j )
                                        +=  ( x_new_full.at( s*_sdim + k ) - _mean_x_new.at(k) ) *_dRdr.at( k*_sdim + t );
                            _dAdr.at( s*_sdim*_hdim + t*_hdim  + _involvedAtoms.at(i)*_sdim + j )
                                    -= _Rot.at(j*_sdim+t) * _nAtomsInv ;
                        }
                    }
                }
            }
        }
    }
    return _dAdr;
}

