/*=========================================================================
 * File         : ABF.h
 * Module       : Simulation
 * Copyright    : (C)opyright 2013
 *                See COPYRIGHT statement in top level directory.
 * Authors      : Behrooz Hashemian
 * Created      : Tue Aug 13, 2013  11:49AM
 * Last modified: Tue Aug 13, 2013  11:50AM
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
#include <numeric>
//#include <functional>
//#include <cmath>
//#include <limits>
using namespace std;

/** \brief Adaptive biasing force (ABF) method to enhance sampling of molecular
 * dynamics simulation and calculate free energy, see \cite Darve2008 for details.
 *
 * In free energy calculations based on thermodynamic integration, it is necessary
 * to compute the derivatives of the free energy as a function of one (scalar case)
 * or several (vector case) order parameters.
 *
 */
class ABF
{
public :
    /** Default Constructor */
    ABF() = default;

    /** Constructor with initialization */
    ABF(bool bApplyABF, bool bLoadHistogram, int ldim, int sdim, int nAtoms, double timeStep,
        int logFrequency, int histFrequency,
        double corralHeight, string corralFilename,
        string coreName, const vector<double> &xiMin, const vector<double> &xiMax, const vector<int> &nBins, int fullSample = 200, int rampSample = 100);
    ABF(bool bApplyABF, bool bLoadHistogram, int ldim, int sdim, int nAtoms, double timeStep,
        int logFrequency, int histFrequency,
        double kWall, const vector<double> &lowerWall, const vector<double> &upperWall,
        string coreName, const vector<double> &xiMin, const vector<double> &xiMax, const vector<int> &nBins, int fullSample = 200, int rampSample = 100);

    /** Default Destructor */
    ~ABF();

    /** Copy constructor (explicitly deleted) */
    ABF( const ABF& ) = delete;

    /** Assignmet operator (explicitly deleted) */
    void operator=( const ABF& ) = delete;

    /** Calculate the collective forces */
    vector<double> &CalculateForce(const int &stepNum,
                                   const vector<double> &mass, const vector<double> &position,
                                   const vector<double> &xiVal, const vector<double> &xiJac,
                                   bool bMaster = true,
                                   const vector<double> &oldLambda = vector<double>(0))
    throw(runtime_error);

    /** Retrieve the value of biaisng force in the master ABF */
    const vector<double>& GetLambda( ) const { return _binLambda; }

    /** Write histograms with postfix appended to their filename */
    void WriteHistograms( string postfix ) throw(runtime_error);

protected :
    bool   _bMaster;
    bool   _bApplyABF;
    bool   _bLoadHistogram;
    int    _ldim;
    int    _sdim;
    int    _hdim;
    int    _nAtoms;
    int    _nSteps  = 0;
    int    _stepNum = 0;
    int    _logFrequency;
    int    _histFrequency;
    double _timeStepInv;


    /* Boundary potentials */
    bool   _bBoundaryPotential = false; ///< apply a potential at the boundary either corral or semi-harmonic
    bool   _bCorral            = false; ///< apply corral potential at the boundary
    bool   _bHarmonic          = false; ///< apply semi-hamonic potential at the boundaries
    double _kWall              = 1.0;   ///< spring stiffness of the semi-hamonic potential
    double _diffLower;                  ///< differce between xi value and the lower wall
    double _diffUpper;                  ///< differce between xi value and the upper wall
    vector<double> _lowerWall;          ///< lower wall tos apply semi-hamonic potential
    vector<double> _upperWall;          ///< upper wall to apply semi-hamonic potential
    vector<double> _xiBoundaryForce;    ///< current  boundary force at the xi
    vector<double> _xiBoundaryForceOld; ///< previous boundary force at the xi
    Corral         _corral;             ///< corral potential object

    /* Atomistic variables*/
    vector<double> _mass;          ///< atoms' masses
    vector<double> _positionOld;   ///< atoms' pisitions  at the previous time step
    vector<double> _velocity;      ///< atoms' velocities at the previous time step
    vector<double> _force;         ///< atoms' forces     at the previous time step


    /* CV's variables*/
    vector<double> _xiMin;         ///< minimum value of \f$\xi\f$ */
    vector<double> _xiMax;         ///< maximum value of \f$\xi\f$ */
    vector<double> _xiRange;       ///< range of \f$\xi\f$ */
    vector<double> _xiVal;         ///< \f$\xi\f$ value */
    vector<double> _xiJac;         ///< \f$\xi\f$ Jacobian */

    vector<double> _xiBin;         ///< the bin of histogram related in which \f$\xi\f$ lies*/
    vector<double> _xiMass;        ///< \f$\xi\f$ mass metrix tensors */
    vector<double> _xiMassInv;     ///< inverse of \f$\xi\f$ mass metrix tensors */
    vector<double> _xiVelocity;    ///< \f$\xi\f$* velocity */
    vector<double> _xiMomentum;    ///< \f$\xi\f$* momentum */
    vector<double> _xiMomentumOld; ///< \f$\xi\f$* old momentum */
    vector<double> _xiForce;       ///< \f$\xi\f$* force to be exerted on the system */
    vector<double> _xiTotalForce;    ///< \f$\xi\f$* force calculated from the simulation*/

    /* ABF parameters */
    int _fullSample = 200;         ///< minimum number of samples in each bin to apply full abf
    int _rampSample = 100;         ///< minimum number of samples in each bin to apply ramp abf

    /* histogram variables */
    int            _totalNumBins = 0;
    double         _binSize;
    vector<int>    _nBins;         ///< number of histogram bins in each direction
    vector<int>    _accNumBins;    ///< accumulated product of number of histogram bins in each direction
    int            _binCounts;
    vector<int>    _histCounts;
    vector<double> _histForces;

    vector<double> _binLambda;     ///< ABF force (lambda)
    vector<double> _binLambdaOld;

    //@{
    /** Variables for writing logs and histograms*/
    string   _coreName;
    string   _logName;
    string   _histName;
    ofstream _logFile;
    ofstream _histFile;
    //@}


    /** Calculate mass in CV space */
    inline void CalculateMass( );

    /** Calculate momentum in CV space */
    inline void CalculateMomentum( );

    /** Calculate ABF lambda force */
    void CalculateLambda() throw(out_of_range);

    /** Calculate boundary potential force */
    void CalculateBoundaryForces();

    /** Build ldim dimensional histograms of forces and counts*/
    void BuildHistograms();

    /** Load ldim dimensional histograms of forces and counts*/
    void LoadHistograms() throw(runtime_error);

    /** ramp function to apply ABF*/
    inline double ramp(int nSamples);

    /** Write log */
    void WriteLogs() throw(runtime_error);

    /** Check if parameters and point coordinates are properly set. */
    bool CheckSettings();

    /** Update the variables and make the class ready to use. */
    void Update() throw( invalid_argument );
private :
};
