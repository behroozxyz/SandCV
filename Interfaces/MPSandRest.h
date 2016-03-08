/*=========================================================================
 * File         : MPSandRest.h
 * Module       : Swig
 * Copyright    : (C)opyright 2013
 *                See COPYRIGHT statement in top level directory.
 * Authors      : Behrooz Hashemian
 * Created      : Thu Oct 17, 2013  11:50AM
 * Last modified: Thu Oct 17, 2013  11:50AM
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
#ifndef MPSandRest_h_
#define MPSandRest_h_

//Sand Libraries
#include "SandCV.h"
#include "PointSetIO.h"

//Intel Math Kernel Library
#include <mkl.h>

//C++ Standard Libraries
#include <map>
#include <set>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <exception>
#include <stdexcept>
#include <memory>
#include <cmath>
#include <limits>
#include <numeric>

class MPSandRest
{
public:
    MPSandRest(int nPartitions,
               // simulation information
               int nAtoms, int sdim, int logFrequency,
               // alignment
               char *xRefFilename, int nInvolvedAtoms, int *involvedAtoms,
               // parametrization
               int ldim, int bDistance,
               char *xiNodesFilename, char *xNodesFilename,
               char *xiSeedsFilename,
               // Markers
               char *xMarkersFilename,
               // Bias force
               char *xCenterFilename, double kBias, double rBias,
               double *xiMin, double *xiMax, double gammaPU, double spacingPU, double tol0PU);

    ~MPSandRest();

    void CalculateForces(int stepNum, double phiAngle, double psiAngle, double *masses, double *positions, double *forces);


    int TransferLambda();

    double                     _kBias;
    double                     _rBias;
    double                     _betaPU;
    vector<shared_ptr<SandCV>> _vsandcv;
    vector<double>             _xMarkers;
    vector<double>             _mark2part;
    vector<double>             _lambdaMaster;
    vector<double>             _lambdaSlave;
    vector<double>             _xMaster;
    vector<double>             _dxMaster;
    vector<double>             _xSlave;
    vector<double>             _dxSlave;
    vector<double>             _normDiff;
    vector<double>             _dPhi;
    vector<double>             _dPhiOverlap;
    vector<double>             _xiOverlap;
    vector<double>             _xCenter;
    vector<vector<int>>        _markerNeighbors;
    shared_ptr<ANNSearch>                 _closestMarkerSearch;
    shared_ptr<ProcrustesSuperimposition> _alignment;

    ofstream _fileXiOverlap;
    ofstream _filePsiOverlap;
    ofstream _fileDPhiOverlap;
    ofstream _fileDiffOverlap;
    ofstream _fileDihedrals;

    bool _bApplyABF    = true;
    bool _bDistance    = false;
    int  _logFrequency = 1;
    int  _ldim         = 0;
    int  _hdim         = 0;
    int  _nAtoms       = 0;
    int  _nPartitions  = 0;


    // Parition of Unity
    int            _pid = -1;
    double         _Z   = 0.0;
    double         _wA  = 0.0;
    vector<double> _diff;
    vector<double> _psi;

private:
};

#endif
