/*=========================================================================
 * File         : Puckering.h
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
#pragma once

//C++ Standard Libraries
#include <vector>
#include <numeric>
#include <limits>
#include <exception>
#include <stdexcept>
#include <algorithm>
#include <functional>
#include <cmath>
#include <iostream>
#include <fstream>
using namespace std;

/** \brief This class provides values and gradients of puckering coordinates for
 * a six-membered ring
 *
 */
class Puckering
{
public :
    /** Default Constructor */
    Puckering();

    /** Default Destructor */
    ~Puckering() = default;

    /** Copy constructor */
    Puckering(const Puckering& inPuckering) = delete;

    /** Assignmet operator */
    void operator=( const Puckering& ) = delete;

    /** Calculate the puckering coordinates as a set of collective variables */
    vector<double> CalculateCV(const vector<double> &X);

    /** Get the Puckering's value */
    const vector<double>& GetValue() const { return _value; }

    /** Get the Puckering's Jacobian */
    const vector<double>& GetJacobian() const { return _Jacobian; }


protected :
    int    _sdim;
    int    _ldim;
    int    _hdim;
    int    _nAtoms;
    double _nAtomsInv;
    double _norm;
    double _R1R1;
    double _R1R2;
    double _R2R2;
    double _temp1;
    double _temp2;

    /** Puckering coordinates */

    double         _qx;
    double         _qy;
    double         _qz;
    double         _dqx;
    double         _dqy;
    double         _dqz;

    double         _Q;
    double         _Phi;
    double         _Theta;
    vector<double> _dQ;
    vector<double> _dTheta;
    vector<double> _dPhi;


    vector<double> _sinx;      ///< sin( pi/3*(j-1))
    vector<double> _cosx;      ///< cos( pi/3*(j-1))
    vector<double> _sin2x;     ///< sin(2pi/3*(j-1))
    vector<double> _cos2x;     ///< cos(2pi/3*(j-1))
    vector<double> _signx;     ///< (-1)^(j-1)

    vector<double> _z;         ///< Centered coordinates
    vector<double> _Delta; ///< Centering matrix
    vector<double> _R;         ///< Centered coordinates
    vector<double> _R1;        ///< R'
    vector<double> _R2;        ///< R''
    vector<double> _R1xR2;     ///< R' x R''

    vector<double> _dR1;
    vector<double> _dR2;
    vector<double> _R2xRj;
    vector<double> _RjxR1;
    vector<double> _dz;

    vector<double> _value;
    vector<double> _Jacobian;


    /** Compute the puckering variable z_j */
    void CalculateZ(const vector<double> &X);

    /** Compute the derivatives of puckering variable z_j */
    void CalculateGradZ(const vector<double> &X);

    /** Calculate the puckering coordinate Q */
    double CalculateQ();

    /** Calculate gradient of the puckering coordinate Q */
    vector<double> CalculateGradQ();

    /** Calculate the puckering coordinate $f\theta$f */
    double CalculateTheta();

    /** Calculate gradient of the puckering coordinate $f\theta$f */
    vector<double> CalculateGradTheta();

    /** Calculate the puckering coordinate $f\phi$f */
    double CalculatePhi();

    /** Calculate gradient of the puckering coordinate $f\phi$f */
    vector<double> CalculateGradPhi();


private :
    inline void outer_product(const vector<double>::const_iterator &iA, const vector<double>::const_iterator &iB,
                       const vector<double>::iterator &iC);

};
