/*=========================================================================
 * File         : SMD.h
 * Module       : Simulation
 * Copyright    : (C)opyright 2013
 *                See COPYRIGHT statement in top level directory.
 * Authors      : Behrooz Hashemian
 * Created      : Mon Oct 07, 2013  12:55PM
 * Last modified: Mon Oct 07, 2013  12:55PM
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
#pragma once

//Sand Libraries
#include "Corral.h"
//Intel Math Kernel Library
#include "mkl.h"
//C++ Standard Libraries
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <exception>
#include <stdexcept>
#include <algorithm>
#include <iterator>
//#include <functional>
//#include <cmath>
//#include <numeric>
//#include <limits>
using namespace std;

/** \brief Steered Molecular Dynamics(SMD)
 *
 */
class SMD
{
public :
    /** Default Constructor */
    SMD() = default;

    /** Constructor with initialization */
    SMD(bool bApplySMD, bool bDistance, bool bWithReturn, int ldim, int sdim, int nAtoms, double timeStep, \
         double pullingRate, double pullingConstant, double pullingMaxLength,
         int logFrequency,  string coreName,
         double kWall, const vector<double> &lowerWall, const vector<double> &upperWall,
         double xiMin, double xiMax);

    /** Default Destructor */
    ~SMD();

    /** Copy constructor (explicitly deleted) */
    SMD( const SMD& ) = delete;

    /** Assignmet operator (explicitly deleted) */
    void operator=( const SMD& ) = delete;

    /** Calculate the collective forces */
    vector<double> &CalculateForce(const vector<double> &position, const vector<double> &xiVal, const vector<double> &xiJac, double cdist = 0.0);

protected :
    bool   _bApplySMD;
    bool   _bDistance;
    int    _sdim;
    int    _ldim;
    int    _bdim;  /// bDistance ? ldim + 1 : ldim
    int    _hdim;
    int    _nAtoms;
    int    _nSteps = 0;
    int    _logFrequency;
    double _timeStep;
    double _cdist;

    /* Boundary potentials */
    bool   _bBoundaryPotential = false; ///< apply a potential at the boundary either corral or semi-harmonic
    bool   _bCorral            = false; ///< apply corral potential at the boundary
    bool   _bHarmonic          = false; ///< apply semi-hamonic potential at the boundaries
    double _corralHeight       = 1.0;   ///< corral height in the simulation energy unit [kcal/mol].
    double _corralWidth        = 1.0;   ///< corral width in the bin size unit.
    double _kWall              = 1.0;   ///< spring stiffness of the semi-hamonic potential
    double _diffLower;                  ///< differce between xi value and the lower wall
    double _diffUpper;                  ///< differce between xi value and the upper wall
    vector<double> _lowerWall;          ///< lower wall tos apply semi-hamonic potential
    vector<double> _upperWall;          ///< upper wall to apply semi-hamonic potential
    vector<double> _xiBoundaryForce;    ///< current  boundary force at the xi
    Corral         _corral;             ///< corral potential object

    /* Atomistic variables*/
    vector<double> _positionOld;    ///< atomic pisitions at the previous time step
    vector<double> _displacement;   ///< displacement of each atom at each time step
    vector<double> _force;          ///< forces to be exreted on the individual atoms
    double         _work;           ///< work done at each time step
    double         _workTotal = 0.0;///< total work done from start of simulation

    /* CV's variables*/
    double _xiMin;           ///< minimum value of \f$\xi\f$ */
    double _xiMax;           ///< maximum value of \f$\xi\f$ */
    double _xiRange;         ///< range of \f$\xi\f$ */
    double _xiPulling;       ///< \f$'xi\f$ value of other end of pulling spring
    vector<double> _xiVal;   ///< \f$\xi\f$ value */
    vector<double> _xiJac;   ///< \f$\xi\f$ Jacobian */
    vector<double> _xiForce; ///< \f$\xi\f$* force to be exerted on the system */

    /* Pulling variable */
    bool   _bWithReturn   = false;
    bool   _bForwardMode  = true;
    bool   _bBackwardMode = false;
    double _pullingRate;
    double _pullingConstant;
    double _pullingDisplacement;
    double _restTimeMiddleStart  = 0.0;
    double _restTimeMiddleFinish = 0.0;
    double _restTimeEnd          = 0.0;
    double _pullingLength;
    double _pullingMaxLength;

    //@{
    /** Variables for writing logs and histograms*/
    string   _coreName;
    string   _logName;
    ofstream _logFile;
    ofstream _histFile;
    //@}

    /** Calculate boundary potential force */
    void CalculateBoundaryForces();

    /** Write log */
    void WriteLogs() throw(runtime_error);

    /** Check if parameters and point coordinates are properly set. */
    bool CheckSettings();

    /** Update the variables and make the class ready to use. */
    void Update() throw( invalid_argument );
private :
};
