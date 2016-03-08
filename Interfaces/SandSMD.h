/*=========================================================================
 * File         : SandSMD.h
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
#ifndef SandSMD_h_
#define SandSMD_h_

//Sand Libraries
#include "SMD.h"
#include "SandCV.h"
//#include "PointSetIO.h"
//C++ Standard Libraries
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <exception>
#include <stdexcept>

class SandSMD
{
public:
    SandSMD(int bApplySMD,
                 // simulation information
                 int nAtoms, int sdim, double timeStep,
                 int logFrequency,
                 char *coreName,
                 // alignment
                 char *xRefFilename, int nInvolvedAtoms, int *involvedAtoms,
                 // parametrization
                 int ldim, int bDistance,
                 char *xiNodesFilename, char *xNodesFilename,
                 char *xiSeedsFilename,
                 // Pulling SMD
                 int bWithReturn,
                 double pullingRate, double pullingConstant, double pullingMaxLength,
                 double kWall, double *lowerWall, double *upperWall,
                 double xiMin, double xiMax
                 );

    ~SandSMD();

    void CalculateForces(double *positions, double *forces);

    SandCV *sandcv;
    SMD    *smd;
    bool   _bApplySMD = true;
    bool   _bDistance = false;
    int    _hdim      = 0;
};

#endif
