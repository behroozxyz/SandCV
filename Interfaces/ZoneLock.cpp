/*=========================================================================
 * File         : ZoneLock.cpp
 * Module       : Swig
 * Copyright    : (C)opyright 2013
 *                See COPYRIGHT statement in top level directory.
 * Authors      : Behrooz Hashemian
 * Created      : Thu Dec 19, 2013  07:44PM
 * Last modified: Thu Dec 19, 2013  07:44PM
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
#include "ZoneLock.h"

ZoneLock::ZoneLock(int nAtoms, int sdim, int ldim,
                    char *xRefFilename, int nInvolvedAtoms, int *involvedAtoms,
                    char *zoneCenterFilename, double kBias, double bBias)
    :
      _ldim(ldim),
      _nAtoms(nAtoms),
      _hdim(sdim*nAtoms),
      _kBias(kBias),
      _bBias(bBias)
{
    int nLines1 = 0;
    int nCols1  = 0;

    _bApplyBias = _kBias > 0.0 ? true : false;

    vector<int>    vInvolvedAtoms(involvedAtoms,involvedAtoms+nInvolvedAtoms);
    vector<double> xRef;
    PointSetIO     reader;

    //______________________________________________________________________________
    // Read xRef
    reader.ReadPointSet(xRefFilename, xRef, nLines1, nCols1);
    if (nLines1 < 1)
        cerr<<"xRefFilename("<<xRefFilename<<") is empty."<<endl;
    else if (nLines1 > 1)
        cerr<<"xRefFilename("<<xRefFilename<<") includes more than one configuration."<<endl;

    _alignment = make_shared<ProcrustesSuperimposition>(sdim,xRef,vInvolvedAtoms);

    //______________________________________________________________________________
    // Read Zone Centers
    reader.ReadPointSet(zoneCenterFilename, _zPosition, _nZones, nCols1);
   if (_nZones < 1)
        cerr<<"zoneCenterFilename("<<zoneCenterFilename<<") is empty."<<endl;
   if (nCols1 != _hdim)
       cerr<<"zoneCenterFilename("<<zoneCenterFilename<<"):"<<
             "Dimension of zone centers("<<nCols1<<") does not match with hdim ("<<_hdim<<")."<<endl;

   _diff .resize(_hdim);
   _force.resize(_hdim);
}


void ZoneLock::CalculateForces( int stepNum, double *masses, double *positions, double *forces )
{
    /* Change the masses and positions array to the STL vectors */
    vector<double> vMasses  (masses   , masses   +_nAtoms);
    vector<double> vPosition(positions, positions+_hdim  );
    _alignment->ComputeDerivatives(vPosition);
    vector<double> aPosition( _alignment->GetAlignedX() );

    /* Compute the CV value of the current configuration and its related Jacobian */
    if(_bApplyBias)
    {
        fill(_force.begin(),_force.end(), 0.0);
        for (int iZ = 0; iZ < _nZones; ++iZ)
        {
            transform(aPosition.begin(), aPosition.end(), _zPosition.begin()+iZ*_hdim, _diff.begin(), std::minus<double>() );
            _dist2 = inner_product(_diff.begin(),_diff.end(),_diff.begin(),0.0);
            cblas_daxpy(_hdim,  _kBias*exp(-_bBias*_dist2)/sqrt(_dist2), _diff.data(), 1, _force.data(), 1);
        }
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 1, _hdim, _hdim, 1.0,
                    _force.data(), _hdim, _alignment->GetdAdX().data(), _hdim, 0.0, forces, _hdim);

    }
    else
        fill( forces, forces+_hdim, 0.0 );

}
