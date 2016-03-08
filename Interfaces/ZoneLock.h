/*=========================================================================
 * File         : ZoneLock.h
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
#ifndef ZoneLock_h_
#define ZoneLock_h_

//Sand Libraries
#include "SandCV.h"
#include "PointSetIO.h"

//Intel Math Kernel Library
#include <mkl.h>

//C++ Standard Libraries
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <exception>
#include <stdexcept>
#include <memory>
#include <numeric>

class ZoneLock
{
public:
    ZoneLock(int nAtoms, int sdim, int ldim,
             char *xRefFilename, int nInvolvedAtoms, int *involvedAtoms,
             char *zoneCenterFilename, double kBias, double bBias);

    ~ZoneLock() = default;

    void CalculateForces(int stepNum, double *masses, double *positions, double *forces);

    shared_ptr<ProcrustesSuperimposition> _alignment;
    vector<double> _zPosition;
    vector<double> _force;
    vector<double> _diff;
    bool   _bApplyBias = true;
    int    _nAtoms     = 0;
    int    _hdim       = 0;
    int    _ldim       = 0;
    int    _nZones     = 0;
    double _dist2      = 0.0;
    double _bBias      = 0.0;
    double _kBias      = 0.0;
};

#endif
