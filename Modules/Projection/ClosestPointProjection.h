/*=========================================================================
 * File         : ClosestPointProjection.h
 * Module       : Parametrization
 * Copyright    : (C)opyright 2013
 *                See COPYRIGHT statement in top level directory.
 * Authors      : Behrooz Hashemian
 * Created      : Mon Jul 29, 2013  04:58PM
 * Last modified: Mon Jul 29, 2013  04:58PM
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
#include "LocalMaxEntropy.h"
#include "PointSetIO.h"

//Intel Math Kernel Library
#include <mkl.h>

//C++ Standard Libraries
#include <vector>
#include <memory>
#include <exception>
#include <stdexcept>
#include <algorithm>
#include <functional>
#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <string>
using namespace std;

/** \class ClosestPointProjection
 * \brief This class finds low-dimensional (d) representaion of a given sample point in high-dimensional(D) space.
 *
 * Given a sample point \f$x_s\f$ in high-dimensional space \f$\mathbb{R}^D\f$,
 * it finds the low-dimensional representation in \f$\mathbb{R}^d\f$ by minimizing
 * \f$\vert \mathcal{M}(\xi_s) - x_s \vert^2 \f$, using a combination of L-BFGS
 * and Newton-Raphon mehtod with line search.
 */
class ClosestPointProjection
{
public :
    /** Default Constructor */
    ClosestPointProjection() = default;

    /** Constructors with initialization */
    ClosestPointProjection(int lowDim, int highDim, double spacing,
                           const vector<double> &xi, const vector<double> &x,
                           const vector<double> &xiSeeds,
                           double tol0=1e-8 , double gamma=1.0)
    throw(invalid_argument);

    ClosestPointProjection(int lowDim, int highDim, double spacing,
                           const vector<double> &xi, const vector<double> &x,
                           const vector<double> &xiSeeds,
                           const vector<double> &xiMarkers,
                           double tol0=1e-8 , double gamma=1.0)
    throw(invalid_argument,runtime_error);

    /** Destructor */
    ~ClosestPointProjection();

    /** Copy constructor */
    ClosestPointProjection(const ClosestPointProjection& inCPP);

    /** Assignmet operator (explicitly deleted) */
    void operator=( const ClosestPointProjection& ) = delete;

    /** different interfaces to find the closest point projection of a given sample point*/
    //@{
    const vector<double> &FindClosestPoint(const vector<double> &xAligned) throw (runtime_error);
    const vector<double> &FindClosestPoint(const vector<double> &xAligned, bool bUseSeed) throw (runtime_error);
    const vector<double> &FindClosestPoint(vector<double>::iterator iter_begin, vector<double>::iterator iter_end) throw (runtime_error);
    //@}
    const vector<double> &GetJacobian()  const { return _xiJac; }
    const vector<double> &GetX()         const { return _xClosest;  }
    const vector<double> &GetDX()        const { return _dxClosest; }
    double GetDistance() const { return _distVal; }
    const vector<double> &GetDistanceJacobian() const { return _distJac; }

    const vector<double> &GetHighSeeds() const { return _xSeeds;  }
    const vector<double> &GetLowSeeds()  const { return _xiSeeds; }

    const vector<double> &GetHighMarkers() const { return _xMarkers;  }

    /** Calculate Jacobian of distance to the manifold with respect to its hig-dimensional representation */
    double CalculateDistance();

protected :
    bool _bTrustRegion   = false;
    bool _bNewtonRaphson = false;
    int _ldim   = 0;           ///< dimension of the low-dimensional space
    int _ldim2  = 0;           ///< the low-dimensional squared
    int _hdim   = 0;           ///< dimension of the high-dimension space
    int _nNodes = 0;           ///< number of node points

    double _sqrtEpsilon;    //!< square root of epsilon
    double _gold;           //!< Golden Ratio: 0.3819660

    double _norm_dxi;
    double _alpha;             ///< line search variable
    double _norm_res;
    vector<double> _xiLS;      ///< xi value in line search
    vector<double> _dxi;       ///< delta xi in Newton method

    double         _distVal;   ///< distance to the manifold
    vector<double> _distJac;   ///< derivatives of distance

    vector<double> _xNodes;    ///< high-dimensional node points
    vector<double> _xAligned;   ///< high-dimensional new point
    vector<double> _xDiff;     ///< diffrence between sample point and its closest projection \f$ \mathcal{M}(\xi) - x_{sample} \f$
    vector<double> _dDiff;     ///< Jacobian of xDiff
    vector<double> _hDiff;     ///< Hessian of xDiff
    vector<double> _xClosest;  ///< the closest point to the given sample point
    vector<double> _dxClosest; ///< gradient of the closest point
    vector<double> _hxClosest; ///< Hessian of the closest point
    vector<double> _xiInit;    ///< Initial point to start for minimization
    vector<double> _xiVal;     ///< low-dimnational parametrization of the closest point
    vector<double> _xiValOld;  ///< low-dimnational parametrization of the closest point
    vector<double> _xiJac;     ///< \f$ \mathbf{D}\xi = \frac{\partial \xi}{\[partial x} \f$
    vector<double> _xiSeeds;   ///< low-dimensional seed points for Newton method
    vector<double> _xSeeds;    ///< high-dimensional seed points for Newton method
    vector<double> _xMarkers;  ///< high-dimensional marker points


    /* Trust region solver variables */
    _TRNSP_HANDLE_t _handle = nullptr; ///< Trust-rigion solver handle
//    _TRNSPBC_HANDLE_t _handle = nullptr; ///< Trust-rigion solver handle (with BC)
    int    _nIters;               ///< number of iterations
    int    _maxIterTR    = 100;   ///< maximum number of iterations
    int    _maxTrialIter = 100;   ///< maximum number of iterations of calculation of trial-step
    int    _stopCriteria;         ///< stop-criterion number
    double _tolTR        = 1.0e-6;///< tolerence of Trust-Region   algorithm for all of its critera
    double _initRes;              ///< initial residuals
    double _finalRes;             ///< final residuals
    vector<double> _lowerBound;   ///< parametrization lower bounds
    vector<double> _upperBound;   ///< parametrization upper bounds
    vector<int>    _info;         ///< results of input parameter checking
    vector<double> _eps;          /**< vector of precisions for stop-criteria:
                                    * eps -[stop-criterion]
                                    * 0 - [2] the trust-region area \f$ \Delta < eps(0)\f$
                                    * 1 - [3] norm of value of the functional \f$ \|F(x)\|_2 < eps(1) \f$
                                    * 2 - [4] singularity of Jacobian matrix \f$ \|J(x)(1:H,j)\| 2 < eps(2), j = 1, ..., n \f$
                                    * 3 - [5] trial step size f$ \|s\|_2 < eps(3) \f$
                                    * 4 - [6] \f$ \|F(x)\|_2 - \|F(x) - J(x)s\|_2 < eps(4) \f$
                                    * 5 -     The trial step precision. If eps(5) = 0 => 1.0e-10
                                    */
    vector<double> _fvec;         ///< function value
    vector<double> _fjac;         ///< Jacobian matrix

    /* Newton-Raphson solver variables */
    int _maxIterNR    = 20;      ///< maximum number of Newton-Raphson itration
    int _maxIterLS    = 100;     ///< maximum number of line search    iteration
    double _tolNR     = 1.0e-6;  ///< tolerence of Newton-Raphson algorithm for all of its critera
    double _tolLS     = 1.0e-14; ///< tolerence of Newton-Raphson algorithm for all of its critera

    /* Local maximum entropy parametrization */
    int             _knn;   ///< number of nearest neighbor points
    double          _spacing;
    vector<int>     _idxNN; ///< indices of the nearest neighbor points
    vector<double>  _pa;    ///< shape functions
    vector<double>  _dpa;   ///< gradient of shape functions
    vector<double>  _hpa;   ///<< Hessian of shape functions
    LocalMaxEntropy _parametrization; ///< local maximum entropy parametrization

    /** High-dimansional searcher to find the best seed. */
    shared_ptr<ANNSearch> _highDimSearch;

    /** Find the closest point by minimizing a cost function */
    int FindClosestPoint( );

    /** Use a high dimensional seed point */
    void UseSeed();

    /** Perform a least square fitting to find closest point on the manifol */
    int MinimizeWithTrustRegion();

    /** The cost function minmized in least square sense */
    int CalculateCostValue();

    /** The gradient of the cost function minmized in least square sense */
    int CalculateCostGradient();

    /** A pointer function which point to difference functions to calculate Jacobian
      * in different dimensions: 1D, 2D or ND */
    int (ClosestPointProjection::* pCorrectJacobian)();

    /** Calculate Jacobian of low-dimensional point with respect to its hig-dimensional representation */
    int CalculateJacobian();

    /** Correct the Jacobian by considering the curvature of the manifold */
    //@{
    int CorrectJacobian1D();
    int CorrectJacobian2D();
    int CorrectJacobianND();
    //@}


    /** Check if parameters and point coordinates are properly set. */
    bool CheckSettings();

    /** Update the variables and make the class ready to use. */
    int Update() throw( invalid_argument );


private :
    bool _bSeed      = true;           /**< flag to determine if we should perform a high-dimensional
                                         * search to find a seed for trust region algorithm */
    int _nTakenSeeds = 0;              ///< keep track of how many time a high-dimensional search is performed
    int _idxSeed;                       ///< the chosen seed's index
    vector<double>::iterator seedIter; ///< iterator to keep the chosen seed position in the related vector

    vector<double> _dxdx;              ///< temporary vector to store inv(dx' dx)
    vector<double> _cor;               ///< temporary vector to store curvature correction variables
    vector<double> _tempJac;           ///< temporary vector to calculate Jacobian of distance

    int MinimizeWithNewtonRaphson(bool applyLS);
    int BrentLineSearch(int maxIter, double TolX, double &alpha);
    double LineSearchFun(double alpha);
    /** Returns the absolute value of A times the sign of B */
    inline double sign(double absA, double signB) const { return (signB < 0) ? -fabs(absA) : fabs(absA); }
};
