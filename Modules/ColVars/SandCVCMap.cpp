/*=========================================================================
 * File         : SandCVCMap.cpp
 * Module       : Simulation
 * Copyright    : (C)opyright 2013
 *                See COPYRIGHT statement in top level directory.
 * Authors      : Behrooz Hashemian
 * Created      : Sun Aug 11, 2013  02:54PM
 * Last modified: Sun Aug 11, 2013  02:54PM
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
#include "SandCVCMap.h"
SandCVCMap::SandCVCMap(int sdim, int ldim, int nAtoms, double spacing, double d0,
               const vector<double> &xiNodes, const vector<double> &xNodes,
               const vector<double> &xiSeeds)
    : _ldim(ldim),
      _hdim(nAtoms*(nAtoms-1)/2),
      _rdim(sdim*nAtoms),
      _Jacobian((ldim+1)*_rdim, 0.0),    //+1 is for the distance as an optional extra CV
      _alignment(sdim,nAtoms,d0),
      _projection(ldim,_hdim,spacing,xiNodes,xNodes,xiSeeds)
{
}


SandCVCMap::SandCVCMap(int sdim, int ldim, int nAtoms, double spacing, double d0,
               const vector<double> &xiNodes, const vector<double> &xNodes,
               const vector<double> &xiSeeds,
               const vector<double> &xiMarkers)
    : _ldim(ldim),
      _hdim(nAtoms*(nAtoms-1)/2),
      _rdim(sdim*nAtoms),
      _Jacobian((ldim+1)*_rdim, 0.0),    //+1 is for the distance as an optional extra CV
      _alignment(sdim,nAtoms,d0),
      _projection(ldim,_hdim,spacing,xiNodes,xNodes,xiSeeds,xiMarkers)
{
}


SandCVCMap::SandCVCMap(const SandCVCMap &inSandCVCMap)
    : _ldim(inSandCVCMap._ldim),
      _hdim(inSandCVCMap._hdim),
      _rdim(inSandCVCMap._rdim),
      _Jacobian(inSandCVCMap._Jacobian),
      _alignment(inSandCVCMap._alignment),
      _projection(inSandCVCMap._projection)
{
}


void SandCVCMap::CalculateCV(const vector<double> &xNew, int stepNum)
throw(runtime_error)
{
    _bUseSeed = _stepNum == stepNum-1 ? false : true;

    try
    {
        _alignment.ComputeDerivatives(xNew);
        _value = _projection.FindClosestPoint( _alignment.GetAlignedX(), _bUseSeed );
        _cdist = _projection.CalculateDistance();
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, _ldim, _rdim, _hdim, 1.0,
                    _projection.GetJacobian().data(), _hdim,
                    _alignment.GetdAdX().data(), _rdim,
                    0.0, _Jacobian.data(), _rdim);
        _stepNum  = stepNum;
    }
    catch(const runtime_error &e)
    {
        cerr<<__func__<<" >> "<<e.what()<<endl;
        throw;
    }
}


void SandCVCMap::CalculateCV(const vector<double> &xNew)
throw(runtime_error)
{
    try
    {
        _alignment.ComputeDerivatives(xNew);
        _value = _projection.FindClosestPoint( _alignment.GetAlignedX() );
        _cdist = _projection.CalculateDistance();
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, _ldim, _rdim, _hdim, 1.0,
                    _projection.GetJacobian().data(), _hdim,
                    _alignment.GetdAdX().data(), _rdim,
                    0.0, _Jacobian.data(), _rdim);
    }
    catch(const runtime_error &e)
    {
        cerr<<__func__<<" >> "<<e.what()<<endl;
        throw;
    }
}

void SandCVCMap::CalculateCVdist(const vector<double> &xNew)
throw(runtime_error)
{
    try
    {
        _alignment.ComputeDerivatives(xNew);
        _value = _projection.FindClosestPoint( _alignment.GetAlignedX() );
        _cdist = _projection.CalculateDistance();
        _value.push_back( _cdist );
        vector<double> tempJacobian = _projection.GetJacobian();
        tempJacobian.insert(tempJacobian.end(),
                            _projection.GetDistanceJacobian().begin() ,
                            _projection.GetDistanceJacobian().end  () );
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, _ldim+1, _rdim, _hdim, 1.0,
                    tempJacobian.data(), _hdim,
                    _alignment.GetdAdX().data(), _rdim,
                    0.0, _Jacobian.data(), _rdim);
    }
    catch(const runtime_error &e)
    {
        cerr<<__func__<<" >> "<<e.what()<<endl;
        throw;
    }
}
