/*=========================================================================
 * File         : SandCMapABF.h
 * Module       : Swig
 * Copyright    : (C)opyright 2013
 *                See COPYRIGHT statement in top level directory.
 * Authors      : Behrooz Hashemian
 * Created      : Tue Aug 13, 2013  11:49AM
 * Last modified: Thu Aug 22, 2013  01:12PM
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
#ifndef SandCMapABF_h_
#define SandCMapABF_h_

//Sand Libraries
#include "ABF.h"
#include "SandCVCMap.h"
#include "PointSetIO.h"

//C++ Standard Libraries
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <exception>
#include <stdexcept>
#include <memory>

class SandCMapABF
{
public:
    SandCMapABF( int bApplyABF,
             // simulation information
             int nAtoms, int sdim, double timeStep,
             int logFrequency, int histFrequency,
             char *coreName,
             // alignment
             char *xRefFilename, int nInvolvedAtoms, int *involvedAtoms,
             // parametrization
             int ldim, int bDistance,
             char *xiNodesFilename, char *xNodesFilename,
             char *xiSeedsFilename,
             // corral potential
             double corralHeight, char *corralFilename,
             // histogram
             int bLoadHistogram , double *xiMin, double *xiMax, int *nBins
             );

    SandCMapABF( int bApplyABF,
             // simulation information
             int nAtoms, int sdim, double timeStep,
             int logFrequency, int histFrequency,
             char *coreName,
             // alignment
             char *xRefFilename, int nInvolvedAtoms, int *involvedAtoms,
             // parametrization
             int ldim, int bDistance,
             char *xiNodesFilename, char *xNodesFilename,
             char *xiSeedsFilename,
             // semi-harmonic potential
             double kWall, double *lowerWall, double *upperWall,
             // histogram
             int bLoadHistogram , double *xiMin, double *xiMax, int *nBins
             );

    ~SandCMapABF() = default;

    void CalculateForces(int stepNum, double *masses, double *positions, double *forces);

    void UpdateSand(int nAtoms, int sdim,
                         char *xRefFilename, int nInvolvedAtoms, int *involvedAtoms,
                         int ldim, char *xiNodesFilename, char *xNodesFilename,
                         char *xiSeedsFilename, const vector<double> &vXiMin, const vector<double> &vXiMax );


    shared_ptr<SandCV> _sandcv;
    shared_ptr<ABF   > _abf;
    bool   _bApplyABF = true;
    bool   _bDistance = false;
    int    _nAtoms    = 0;
    int    _hdim      = 0;
};

#endif
