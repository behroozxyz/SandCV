/*=========================================================================
 * File         : LocalMaxEntropy.h
 * Module       : Parametrization
 * Copyright    : (C)opyright 2013
 *                See COPYRIGHT statement in top level directory.
 * Authors      : Behrooz Hashemian
 * Created      : Wed Jun 22, 2013  11:38AM
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
#include "ANNSearch.h"

//Intel Math Kernel Library
#include "mkl.h"

//C++ Standard Libraries
#include <vector>
#include <memory>
#include <algorithm>
#include <functional>
#include <cmath>
#include <limits>
#include <iostream>
using namespace std;

/** \class LocalMaxEntropy
 * \brief This class computes the local maximum-entropy (LME) shape functions and
 * their first and second-order spation derivatives \cite Arroyo2006.
 *
 * For a given sapmle point \f$\xi_s\f$, it calculates \f$ p_a(\xi_s), \nabla p_a(\xi_s),
 * \nabla\nabla p_a(\xi_s) \f$.
 * It employs Newtorn-Raphson method in possible combination with line search.
 */
class LocalMaxEntropy
{
public :
     /** Default Constructor */
    LocalMaxEntropy() = default;

    /** Destructor */
    ~LocalMaxEntropy() = default;

    /** Copy constructor */
    LocalMaxEntropy(const LocalMaxEntropy& inLME);

    /** Assignmet operator (explicitly deleted) */
    void operator=( const LocalMaxEntropy& ) = delete;

    LocalMaxEntropy(int dim,             const vector<double> &xNodes, double spacing, double tol0=1e-6, double gamma=1.0);
    LocalMaxEntropy(int dim, int nNodes, const vector<double> &xNodes, double spacing, double tol0=1e-6, double gamma=1.0);
    LocalMaxEntropy(int dim, int nNodes,         const double *xNodes, double spacing, double tol0=1e-6, double gamma=1.0);

    /** Compute the LME Shape functions and their first and second spatial derivatives
     * with a fixed \f$ \beta \f$ for a given sample point \f$ \xi_s \f$.
     *
     * \param [in] xNeighbors    index list of the nearest neighbors of the sample point\f$ \xi_s \f$
     * \param [in] xSample       sample point at where the LME shapefunctions are computed
     *
     * \return erro identificator
     */
    int ComputeLme(const vector<double> &xSample, const vector<int> &idxNN , bool bGrad, bool bHessian);

    /** One-dimensional LME */
    int ComputeLme1D(const vector<double> &xSample, const vector<int> &idxNN , bool withGradient, bool withHessian);

    /** One-dimensional LME */
    int ComputeLme2D(const vector<double> &xSample, const vector<int> &idxNN , bool bGrad, bool bHessian);

    /** Compute the LME Shape functions and their requested derivatives */

    int ComputeLme( int dim, double *xSample, int derivativeOrder);
    int ComputeLme(vector<double>::const_iterator iter_begin, vector<double>::const_iterator iter_end, int derivativeOrder);
    int ComputeLme(const vector<double> &xSample, int derivativeOrder);
    int ComputeLme(const vector<double> &xSample,const vector<int> &idxNN, int derivativeOrder);
    /** Compute only the LME Shape functions and return them as a vector */
    const vector<double> &ComputeLme (const vector<double> &xSample, const vector<int> &idxNN);

    /** Get the vector of indices of the nearest neighbor nodes */
    const vector<int>& GetNearestNeighbors( ) { return _idxNN; }

    /** Get the vector of calculated sahpe functions */
    const vector<double>& GetShapeFunctions( ) { return _pa; }

    /** Get the vector of calculated gradient of sahpe functions */
    const vector<double>& GetShapeFunctionsGradient( ) { return _dpa; }

    /** Get the vector of calculated Hessian of sahpe functions */
    const vector<double>& GetShapeFunctionsHessian( ) { return _hpa; }

protected:
    int    _ierr      = 0;      ///< error identifier
    int    _dim       = 0;      ///< spatial dimension
    int    _dim2      = 0;      ///< dim*dim
    int    _knn       = 0;      ///< number of nearest neighbors
    double _range     = 0.0;    ///< range for neibor search (where the shape functions decay to less than _tol0)
    double _beta      = 0.0;    ///< locality LME parameter
    double _h         = 0.0;    ///< spacing between adajenct nodes
    double _gamma     = 1.0;    ///< \f$ \frac{\beta}{h^2} \f$

    // Minimization variables
    int _maxIterNR    = 100;    ///< maximum number of Newton-Raphson itration
    int _maxIterLS    = 1000;   ///< maximum number of line search    iteration
    double _tol0      = 1e-8;   ///< Tolerence to define a cut-off for nearest neighbor search \f$ \vert x - x_a \vert \le \sqrt{\frac{\log(Tol_0)}{\beta}} \f$
    double _tolNR     = 1e-14;  ///< Tolerence in \f$\lambda\f$ for Newton-Raphson method
    double _tolLS     = 1e-10;  ///< Tolerence in \f$\lambda\f$ for Newton-Raphson method
    double _sqrtEpsilon;        ///< square root of epsilon
    double _gold;               ///< Golden Ratio: 0.3819660

    // Line search bounds
    double _lowerLimit = -0.01; ///< lower limit of the line search minimization
    double _upperLimit =  5.00; ///< upper limit of the line search minimization

    vector<double> _xNodes;     ///< spatial coordinates for the nodes
    int            _nNodes;      /// total number of node points
    vector<double> _xSample;    ///< the point in wich shape functions are calculated
    vector<double> _pa;         ///< LME shape functions at a given sample point
    vector<double> _dpa;        ///< LME gradient at a given sample point
    vector<double> _hpa;        ///< LME hessian at a given sample point

    // Nearest neighbor search
    vector<int>           _idxNN;    ///< Nearest neghbor indices
    shared_ptr<ANNSearch> _searcher; ///< Nearest neighbor search object


    /** Calculate local maximum-entropy shape functions for a given sample point. */
    int ComputeShapeFunctions();


    /** This function computes the "original" gradients for one sample point (beta is fix).
        * dp_a contains the values of the nDim spacial derivatives of the shape
        * functions corresponding to the closest neighbor nodes at the sample point x (gradient)*/
    void ComputeGradients();

    /** Computation of the "original" local hessian matrix for p_a (beta is fix).
        * hp_a contains the values of the nDim second order spacial derivatives of the shape
        * functions corresponding to the closest neighbor nodes at the sample point x (hessian)*/
    void ComputeHessians();

    /** Computate of log(Z), which is needed in the line search
     * \return log(Z)
     */
    double LogZgeneral (double _alpha);

    /**
     * Local-maximum entropy shape functions main loop
     * \param applyLS    line search flag
     * \param step       At first a checking for the step size of the JUMP is made
     * \return error identifier:
     *                       0 - successful
     *                       1 - error in inverting J
     *                       2 - maximum number of iteration reached
     *                       3 - there is NaN in the res
     */
    int LmeNewton(bool applyLS, double step);


    /** Perform line search with Brent's minimization algorithm (based on \cite Press2007)*/
    int BrentLineSearch(int maxIter, double TolX, double &alpha);

    /** Returns the absolute value of A times the sign of B */
    inline double sign(double absA, double signB) const { return (signB < 0) ? -fabs(absA) : fabs(absA); }

    /** Check if parameters and point coordinates are properly set. */
    bool CheckSettings();

    /** Update the variables and make the class ready to use. */
    int Update() throw( exception );

    //--------------------------------------------------------------------------
    // one dimensional
    //--------------------------------------------------------------------------
    /** Compute shape functions in one-dimensional space*/
    int  ComputeShapeFunctions1D();
    void ComputeGradients1D();
    void ComputeHessians1D();
    double LogZ1D(double lambda);
    int LmeNewton1D(bool applyLS, double step);
    int BrentLineSearch1D(int maxIter, double TolX, double &alpha);

    //--------------------------------------------------------------------------
    // two dimensional
    //--------------------------------------------------------------------------
    /** Compute shape functions in one-dimensional space*/
    int  ComputeShapeFunctions2D();
    void ComputeGradients2D();
    void ComputeHessians2D();
    double LogZ2D(double lambda);
    int LmeNewton2D(bool applyLS, double step);
    int BrentLineSearch2D(int maxIter, double TolX, double &alpha);

private:
    double _Z;          ///< partition function
    double _alpha;      ///< step size in Newton method
    double _norm_res;
    double _norm_dlam;
    double _detJ;
    // 1D variables
    double _lam1D;
    double _dlam1D;
    double _lamLS1D;
    double _res1D;
    double _J1D;
    // vectors
    vector<double> _dx;       ///< matrix with the difference of positions between a sample point and its nearest neighbors.
                              /**< An \f$ N \times D \f$ matrix of the difference between \f$\xi\f$ and \f$\xi_a (a=1,\dots,N)\f$,
                               * where \f$\xi\f$ is the sample point for which the shape functions are evaluated and
                               * \f$\xi_a (a=1,\dots,N) \f$ are the nears neighbor nodes to \f$\xi\f$.*/
    vector<double> _dxdx;     ///< tensor product of: (x-x_a) x (x-x_a)
    vector<double> _beta_dx2; ///< exp(-beta_i*|dx_i|^2)
    vector<double> _lam;      ///< lagrange multiplier of a given sample point
    vector<double> _dlam;     ///< increment of lambda
    vector<double> _lamLS;    ///< lambda vector (only used when line search is applied);
    vector<double> _lam_dx;   ///< inner product of lambda with dx(i,:)
    vector<double> _res;      ///< gradient of the objetive function [d x 1]
    vector<double> _R;        ///< tensor product of res [d x d]
    vector<double> _Jchol;    ///< Hessian of the objective function,_J, or its Cholesky factorization, L, (J = L'*L)
    vector<double> _Jidx;     ///< matrix of inv(J)*dx
    vector<double> _M;        ///< matrix of pa * j_a x j_a
};
