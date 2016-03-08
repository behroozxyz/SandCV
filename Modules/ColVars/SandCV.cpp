/*=========================================================================
 * File         : SandCV.cpp
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
#include "SandCV.h"
SandCV::SandCV(int sdim, int ldim, int hdim, double spacing,
               const vector<double> &xRef   , const vector<int> &activeList,
               const vector<double> &xiNodes, const vector<double> &xNodes,
               const vector<double> &xiSeeds)
    : _ldim(ldim),
      _hdim(hdim),
      _Jacobian((ldim+1)*hdim, 0.0),    //+1 is for the distance as an optional extra CV
      _alignment(sdim,xRef,activeList),
      _projection(ldim,hdim,spacing,xiNodes,xNodes,xiSeeds)
{
}

SandCV::SandCV(int sdim, int ldim, int hdim, double spacing,
               const vector<double> &xRef   , const vector<int> &activeList,
               const vector<double> &xiNodes, const vector<double> &xNodes,
               const vector<double> &xiSeeds,
               const vector<double> &xiMarkers)
    : _ldim(ldim),
      _hdim(hdim),
      _Jacobian((ldim+1)*hdim, 0.0),    //+1 is for the distance as an optional extra CV
      _alignment(sdim,xRef,activeList),
      _projection(ldim,hdim,spacing,xiNodes,xNodes,xiSeeds,xiMarkers)
{
}

SandCV::SandCV(const SandCV &inSandCV)
    : _ldim(inSandCV._ldim),
      _hdim(inSandCV._hdim),
      _Jacobian(inSandCV._Jacobian),
      _alignment(inSandCV._alignment),
      _projection(inSandCV._projection)
{
}


void SandCV::CalculateCV(const vector<double> &xNew, int stepNum)
throw(runtime_error)
{
    _bUseSeed = _stepNum == stepNum-1 ? false : true;

    try
    {
        _alignment.ComputeDerivatives(xNew);
        _value = _projection.FindClosestPoint( _alignment.GetAlignedX(), _bUseSeed );
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, _ldim, _hdim, _hdim, 1.0,
                    _projection.GetJacobian().data(), _hdim,
                    _alignment.GetdAdX().data(), _hdim,
                    0.0, _Jacobian.data(), _hdim);
        _stepNum  = stepNum;
    }
    catch(const runtime_error &e)
    {
        cerr<<__func__<<" >> "<<e.what()<<endl;
        throw;
    }
}


void SandCV::CalculateCV(const vector<double> &xNew)
throw(runtime_error)
{
    try
    {
        _alignment.ComputeDerivatives(xNew);
        _value = _projection.FindClosestPoint( _alignment.GetAlignedX() );
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, _ldim, _hdim, _hdim, 1.0,
                    _projection.GetJacobian().data(), _hdim,
                    _alignment.GetdAdX().data(), _hdim,
                    0.0, _Jacobian.data(), _hdim);
    }
    catch(const runtime_error &e)
    {
        cerr<<__func__<<" >> "<<e.what()<<endl;
        throw;
    }
}

void SandCV::CalculateCVdist(const vector<double> &xNew)
throw(runtime_error)
{
    try
    {
        _alignment.ComputeDerivatives(xNew);
        _value = _projection.FindClosestPoint( _alignment.GetAlignedX() );
        _value.push_back( _projection.CalculateDistance() );
        vector<double> tempJacobian = _projection.GetJacobian();
        tempJacobian.insert(tempJacobian.end(),
                            _projection.GetDistanceJacobian().begin() ,
                            _projection.GetDistanceJacobian().end  () );
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, _ldim+1, _hdim, _hdim, 1.0,
                    tempJacobian.data(), _hdim, _alignment.GetdAdX().data(), _hdim,
                    0.0, _Jacobian.data(), _hdim);
    }
    catch(const runtime_error &e)
    {
        cerr<<__func__<<" >> "<<e.what()<<endl;
        throw;
    }
}
