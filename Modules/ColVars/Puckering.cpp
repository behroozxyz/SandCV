/*=========================================================================
 * File         : Puckering.cpp
 * Module       : Simulation
 * Copyright    : (C)opyright 2013
 *                See COPYRIGHT statement in top level directory.
 * Authors      : Behrooz Hashemian
 * Created      : Mon Dec 16, 2013  03:40PM
 * Last modified: Mon Dec 16, 2013  03:41PM
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
#include "Puckering.h"

Puckering::Puckering() :
    _sdim(3),
    _ldim(2),
    _nAtoms(6),
    _sinx ({0.0, 0.5*sqrt(3.0), 0.5*sqrt(3.0), 0.0,-0.5*sqrt(3.0),-0.5*sqrt(3.0)}), // sin( pi/3*(j-1))
    _cosx ({1.0, 0.5          ,-0.5          ,-1.0,-0.5          , 0.5          }), // cos( pi/3*(j-1))
    _sin2x({0.0, 0.5*sqrt(3.0),-0.5*sqrt(3.0), 0.0, 0.5*sqrt(3.0),-0.5*sqrt(3.0)}), // sin(2pi/3*(j-1))
    _cos2x({1.0,-0.5          ,-0.5          , 1.0,-0.5          ,-0.5          }), // cos(2pi/3*(j-1))
    _signx({1.0,-1.0          , 1.0          ,-1.0, 1.0          ,-1.0          })  // (-1)^(j-1)
{
    _hdim      = _nAtoms*_sdim;
    _nAtomsInv = 1.0 / double(_nAtoms);
    /* Build centering matrix  */
    _Delta.resize(_nAtoms*_nAtoms, 0.0);
    for (auto &iter : _Delta)
        iter = -_nAtomsInv;
    for (int iN = 0; iN < _nAtoms; ++iN)
        _Delta.at(iN*_nAtoms+iN) += 1.0;

    _z     .resize(_nAtoms      , 0.0);
    _R     .resize(_hdim        , 0.0);
    _R1    .resize(_sdim        , 0.0);
    _R2    .resize(_sdim        , 0.0);
    _R1xR2 .resize(_sdim        , 0.0);
    _dR1   .resize(_nAtoms      , 0.0);
    _dR2   .resize(_nAtoms      , 0.0);
    _R2xRj .resize(_hdim        , 0.0);
    _RjxR1 .resize(_hdim        , 0.0);

    _dz    .resize(_nAtoms*_hdim, 0.0);
    _dQ    .resize(_hdim        , 0.0);
    _dTheta.resize(_hdim        , 0.0);
    _dPhi  .resize(_hdim        , 0.0);

    _value   .resize(_ldim      , 0.0);
    _Jacobian.resize(_ldim*_hdim, 0.0);
}



vector<double> Puckering::CalculateCV(const vector<double> &X)
{
    this->CalculateGradZ(X);

    int jN;
    _qx = 0.0;
    _qy = 0.0;
    _qz = 0.0;
    for( jN = 0; jN < _nAtoms; ++jN )
    {
        _qx -= _z.at(jN) * _sin2x.at(jN);
        _qy += _z.at(jN) * _cos2x.at(jN);
        _qz += _z.at(jN) * _signx.at(jN);
    }
//    _qx *= sqrt(2.0*_nAtomsInv);
//    _qy *= sqrt(2.0*_nAtomsInv);
//    _qz *= sqrt(    _nAtomsInv);

    this->CalculateGradQ();
    this->CalculateGradTheta();
    this->CalculateGradPhi();

    _value.at(0) = _Theta;
    _value.at(1) = _Phi;
//    _value.at(2) = _Q;

    copy(_dTheta.begin(), _dTheta.end(), _Jacobian.begin());
    copy(_dPhi  .begin(), _dPhi  .end(), _Jacobian.begin()+_hdim);
//    copy(_dQ    .begin(), _dQ    .end(), _Jacobian.begin()+_dTheta.size()+_dPhi.size());

    return _value;
}


void Puckering::CalculateZ(const vector<double> &X)
{
    int jN, iN, iD;
    /* Coordinates of ring elements with respect to the geometrical center */
    fill(_R.begin(), _R.end(), 0.0);
    for ( jN = 0; jN < _nAtoms; ++jN)
        for (iD = 0; iD < _sdim; ++iD)
            for ( iN = 0; iN < _nAtoms; ++iN)
                _R.at(jN*_sdim+iD) += X.at(iN*_sdim+iD)*_Delta.at(iN*_nAtoms+jN);


    fill(_R1.begin(), _R1.end(), 0.0);
    fill(_R2.begin(), _R2.end(), 0.0);

    /* R1 = sinx*R  and  R2 = sinx*R */
    for (iD = 0; iD < _sdim; ++iD)
        for ( jN = 0; jN < _nAtoms; ++jN)
        {
            _R1.at(iD) += _sinx.at(jN)*_R.at(jN*_sdim+iD);
            _R2.at(iD) += _cosx.at(jN)*_R.at(jN*_sdim+iD);
        }

    /* R1xR2 : Cross product of R1 and R2  */
    outer_product(_R1.begin(), _R2.begin(), _R1xR2.begin());

    /* |R1xR2| : L2 norm of cross product of R1 and R2  */
    _norm = sqrt( _R1xR2.at(0)*_R1xR2.at(0) + _R1xR2.at(1)*_R1xR2.at(1) + _R1xR2.at(2)*_R1xR2.at(2));

    /* Normal vector to the mean plane is n = R1xR2/|R1xR2|,
     * so z = R.n */
    for ( jN = 0; jN < _nAtoms; ++jN)
        _z.at(jN) = inner_product( _R1xR2.begin(), _R1xR2.end(), _R.begin()+jN*_sdim, 0.0) / _norm;
}


void Puckering::CalculateGradZ(const vector<double> &X)
{
    /* Calculate puckering variable z_j */
    this->CalculateZ(X);

    int iN, jN, iD;

    fill(_dR1.begin(), _dR1.end(), 0.0);
    fill(_dR2.begin(), _dR2.end(), 0.0);

    /* Calculate derivatives of R' and R'' */
    for ( iN = 0; iN < _nAtoms; ++iN )
        for ( jN = 0; jN < _nAtoms; ++jN )
        {
            _dR1.at(iN) += _Delta.at(iN*_nAtoms+jN)*_sinx.at(jN);
            _dR2.at(iN) += _Delta.at(iN*_nAtoms+jN)*_cosx.at(jN);
        }

    /* Calculate some dot products */
    _R1R1 = inner_product(_R1.begin(), _R1.end(), _R1.begin(), 0.0);
    _R2R2 = inner_product(_R2.begin(), _R2.end(), _R2.begin(), 0.0);
    _R1R2 = inner_product(_R1.begin(), _R1.end(), _R2.begin(), 0.0);

    /* Calculate some cross products */
    _R1R1 = inner_product(_R1.begin(), _R1.end(), _R1.begin(), 0.0);
    int jj = 0;
    for ( jN = 0; jN < _nAtoms; ++jN )
    {
        jj = jN*_sdim;
        outer_product(_R2.begin()   , _R .begin()+jj, _R2xRj.begin()+jj);
        outer_product(_R .begin()+jj, _R1.begin()   , _RjxR1.begin()+jj);
    }

    /* Calculate derivativs of z_j */
    double norm2Inv = 0.5/(_norm*_norm);
    double normInv  = 1.0/_norm;
    for ( iD = 0; iD < _sdim; ++iD )
        for ( iN = 0; iN < _nAtoms; ++iN )
        {
            _temp1 = 2.0 * (  _dR1.at(iN)*( _R1.at(iD)*_R2R2 - _R2.at(iD)*_R1R2 )
                                     +_dR2.at(iN)*( _R2.at(iD)*_R1R1 - _R1.at(iD)*_R1R2 ) );

            for ( jN = 0; jN < _nAtoms; ++jN )
            {
                _temp2 = _Delta.at(iN*_nAtoms+jN) *_R1xR2.at(iD)
                        + _dR1.at(iN)*_R2xRj.at(jN*_sdim+iD)
                        + _dR2.at(iN)*_RjxR1.at(jN*_sdim+iD);
                // effective computation of grad_i z_j
                _dz.at(iN*_hdim+jN*_sdim+iD) = _temp2*normInv -_z.at(jN)*_temp1*norm2Inv ;
            }
        }
}


double Puckering::CalculateQ()
{
    _Q = sqrt(inner_product(_z.begin(), _z.end(), _z.begin(), 0.0));
    return _Q;
}


vector<double> Puckering::CalculateGradQ()
{
    /* Calculate puckering coordinate Q*/
    this->CalculateQ();

    int iN, jN, iD;
    fill(_dQ.begin(), _dQ.end(), 0.0);

    for( iN = 0; iN < _nAtoms; ++iN )
        for( iD = 0; iD < _sdim; ++iD)
        {
            for( jN = 0; jN < _nAtoms; ++jN )
                _dQ.at(iN*_sdim+iD) += _z.at(jN) * _dz.at(iN*_hdim+jN*_sdim+iD);
            _dQ.at(iN*_sdim+iD) /= _Q;
        }

    return _dQ;
}


double Puckering::CalculateTheta()
{
    _Theta = acos( sqrt(_nAtomsInv)*_qz / _Q );
    return _Theta;
}


vector<double> Puckering::CalculateGradTheta()
{
    /* Calculate puckering coordinate Theta */
    this->CalculateTheta();

    int iN, jN, iD;
    double qzQi    = _qz/_Q;
    double tempInv =  1.0/sqrt(double(_nAtoms)*_Q*_Q-_qz*_qz);
    fill(_dTheta.begin(), _dTheta.end(), 0.0);

    for( iN = 0; iN < _nAtoms; ++iN )
        for( iD = 0; iD < _sdim; ++iD)
        {
            _dqz = 0.0;
            for ( jN = 0; jN < _nAtoms; ++jN)
                _dqz += _dz.at(iN*_hdim+jN*_sdim+iD) * _signx.at(jN);
//            dqz *= sqrt(_nAtomsInv);

            _dTheta.at(iN*_sdim+iD) = (_dQ.at(iN*_sdim+iD)*qzQi - _dqz ) * tempInv;
        }

    return _dTheta;
}


double Puckering::CalculatePhi()
{
    _Phi = atan2(_qx,_qy);

    /* Make Phi in the range of [0,2pi) */
    if (_Phi < 0.0)
        _Phi += 2.0*M_PI;

    return _Phi;
}


vector<double> Puckering::CalculateGradPhi()
{
    /* Calculate puckering coordinate Phi */
    this->CalculatePhi();

    int iN, jN, iD;
    double tempInv = 1.0 / (_qy*_qy+_qx*_qx);
    fill(_dPhi.begin(), _dPhi.end(), 0.0);

    for( iN = 0; iN < _nAtoms; ++iN )
        for( iD = 0; iD < _sdim; ++iD)
        {
            _dqx = 0.0;
            _dqy = 0.0;
            for( jN = 0; jN < _nAtoms; ++jN )
            {
                _dqx -= _dz.at(iN*_hdim+jN*_sdim+iD) * _sin2x.at(jN);
                _dqy += _dz.at(iN*_hdim+jN*_sdim+iD) * _cos2x.at(jN);
            }
//            dqx *= sqrt(2.0*_nAtomsInv);
//            dqy *= sqrt(2.0*_nAtomsInv);

            _dPhi.at(iN*_sdim+iD) = ( _dqx*_qy -_dqy*_qx ) * tempInv;
        }

    return _dPhi;
}



inline void Puckering::outer_product(const vector<double>::const_iterator &iA,
                                     const vector<double>::const_iterator &iB,
                                     const vector<double>::iterator       &iC)
{
    *(iC+0) = *(iA+1)*(*(iB+2)) - *(iA+2)*(*(iB+1));
    *(iC+1) = *(iA+2)*(*(iB+0)) - *(iA+0)*(*(iB+2));
    *(iC+2) = *(iA+0)*(*(iB+1)) - *(iA+1)*(*(iB+0));
}

