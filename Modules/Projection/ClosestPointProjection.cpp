/*=========================================================================
 * File         : ClosestPointProjection.cpp
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
#include "ClosestPointProjection.h"
static const double _epsilon = numeric_limits<double>::epsilon(); ///< Machine precision

ClosestPointProjection::ClosestPointProjection(int lowDim, int highDim,
                                               double spacing,
                                               const vector<double> &xiNodes,
                                               const vector<double> &xNodes,
                                               const vector<double> &xiSeeds,
                                               double tol0, double gamma)
throw(invalid_argument)
    : _ldim(lowDim),
      _hdim(highDim),
      _xNodes(xNodes),
      _xiSeeds(xiSeeds),
      _spacing(spacing),
      _parametrization(lowDim, xiNodes, spacing, tol0, gamma)
{
//    _bTrustRegion   = true;
    _bNewtonRaphson = true;
    _nNodes  = xNodes.size()/_hdim;

    this->Update();
}

ClosestPointProjection::ClosestPointProjection(int lowDim, int highDim,
                                               double spacing,
                                               const vector<double> &xiNodes,
                                               const vector<double> &xNodes,
                                               const vector<double> &xiSeeds,
                                               const vector<double> &xiMarkers,
                                               double tol0, double gamma)
throw(invalid_argument,runtime_error)
    : _ldim(lowDim),
      _hdim(highDim),
      _xNodes(xNodes),
      _xiSeeds(xiSeeds),
      _spacing(spacing),
      _parametrization(lowDim, xiNodes, spacing, tol0, gamma)
{
//    _bTrustRegion   = true;
    _bNewtonRaphson = true;
    _nNodes  = xNodes.size()/_hdim;

    this->Update();

    /* Calculate high-dimensional representation of the low-dimensional markers */
    int nLMarkers = xiMarkers.size()/_ldim;
    _xMarkers.reserve(nLMarkers*_hdim);
    for ( auto iter = xiMarkers.begin(); iter != xiMarkers.end(); iter += _ldim )
    {
        if(_parametrization.ComputeLme( iter, iter+_ldim, 0))
        {
            cerr<<"xiMarker = [ ";
            for (auto iter2 = iter; iter2 != iter+_ldim; ++iter2)
                cerr<<*iter2<<" ";
            cerr<<"] is failed."<<endl;
            throw runtime_error("Failed to calculate xMarkers!");
        }
        else
        {
            _pa    = _parametrization.GetShapeFunctions();
            _idxNN = _parametrization.GetNearestNeighbors();
            _knn   = _idxNN.size();
            fill(_xClosest.begin(), _xClosest.end(), 0.0);
            for ( int iN = 0; iN < _knn; ++iN )
                cblas_daxpy( _hdim, _pa.at(iN), &_xNodes.at(_idxNN.at(iN)*_hdim), 1,
                             _xClosest.data(), 1);
            _xMarkers.insert( _xMarkers.end(), _xClosest.begin(), _xClosest.end() );
        }
    }
    int nHMakers = _xMarkers.size()/_hdim;
    if ( nHMakers != nLMarkers)
        cerr<<"Number of markers in high dimension ("<<nHMakers<<") and low dimension ("<<nLMarkers<<") does not match."<<endl;
}


ClosestPointProjection::~ClosestPointProjection()
{
    cout<<_nTakenSeeds<<" high-dimensional search is performed."<<endl;
}

ClosestPointProjection::ClosestPointProjection(const ClosestPointProjection &inCPP)
    : _bTrustRegion(inCPP._bTrustRegion),
      _bNewtonRaphson(inCPP._bNewtonRaphson),
      _ldim(inCPP._ldim),
      _hdim(inCPP._hdim),
      _nNodes(inCPP._nNodes),
      _xNodes(inCPP._xNodes),
      _xiSeeds(inCPP._xiSeeds),
      _spacing(inCPP._spacing),
      _parametrization(inCPP._parametrization)
{
    this->Update();
}



int ClosestPointProjection::Update( )
throw( invalid_argument )
{
    if ( this->CheckSettings() )
    {
        cerr<<__func__<<" >> CheckSettings went wrong."<<endl;
        throw invalid_argument("One or more inputs are wrong.");
    }

    cout<<"SPACING = "<<_spacing<<endl;

    _sqrtEpsilon = sqrt( _epsilon );      // square root of epsilon
    _gold        = 0.5 * (3.0-sqrt(5.0)); // Golden Ratio: 0.3819660

    if (_ldim == 1)
        pCorrectJacobian = &ClosestPointProjection::CorrectJacobian1D;
    else if (_ldim == 2)
        pCorrectJacobian = &ClosestPointProjection::CorrectJacobian2D;
    else
        pCorrectJacobian = &ClosestPointProjection::CorrectJacobianND;

    _ldim2 = _ldim*_ldim;

    /* Newton-Raphson variables*/
    _xiLS     .resize( _ldim        );
    _dxi      .resize( _ldim        );
    _dDiff    .resize( _ldim        );
    _hDiff    .resize( _ldim2       );

    /* xi and x vectors */
    _xiInit   .resize( _ldim        );
    _xiVal    .resize( _ldim        );
    _xiValOld .resize( _ldim        );
    _xiJac    .resize( _ldim *_hdim );
    _xAligned .resize(        _hdim );
    _xDiff    .resize(        _hdim );
    _xClosest .resize(        _hdim );
    _dxClosest.resize( _ldim *_hdim );
    _hxClosest.resize( _ldim2*_hdim );

    /* Temporary variables to calculate Jacobian */
    _dxdx     .resize( _ldim2 );
    _cor      .resize( _ldim2 );

    /* distance calculation vectors */
    _distJac  .resize(_hdim      );
    _tempJac  .resize(_hdim*_hdim);

    /* Trust region solver variables */
    _fvec     .resize( _hdim );
    _fjac     .resize( _hdim*_ldim );
    _info     .resize( 6 );
    _eps      .resize( 6 );
    _eps.at(0) = sqrt(_epsilon); // [2] the trust-region area | Delta | < eps(0)
    _eps.at(1) = sqrt(_epsilon); // [3] norm of value of the functional |F(x)| < eps(1)
    _eps.at(2) = sqrt(_epsilon); // [4] singularity of Jacobian matrix |J(x)(1:H,j)| < eps(2), j = 1, ..., n
    _eps.at(3) = sqrt(_epsilon); // [5] trial step size f$ |s| < eps(3)
    _eps.at(4) = _tolTR;         // [6] |F(x)| - |F(x) - J(x)s| < eps(4)
    _eps.at(5) = sqrt(_epsilon);            //  The trial step precision. If eps(5) = 0 => 1.0e-10

    /* Calculate high-dimensional representation of the low-dimensional seeds */
    int nLSeeds = _xiSeeds.size()/_ldim;
    _xSeeds.reserve(nLSeeds*_hdim);
    int counter = 0;
    for ( auto iter = _xiSeeds.begin(); iter != _xiSeeds.end(); )
    {
        if(_parametrization.ComputeLme( iter, iter+_ldim, 0))
        {
            cerr<<"xiSeeds("<<counter<<") = [ ";
            for (auto iter2 = iter; iter2 != iter+_ldim; ++iter2)
                cerr<<*iter2<<" ";
            cerr<<"] is failed."<<endl;
            iter = _xiSeeds.erase(iter,iter+_ldim);
        }
        else
        {
            _pa    = _parametrization.GetShapeFunctions();
            _idxNN = _parametrization.GetNearestNeighbors();
            _knn   = _idxNN.size();
            fill(_xClosest.begin(), _xClosest.end(), 0.0);
            for ( int iN = 0; iN < _knn; ++iN )
                cblas_daxpy( _hdim, _pa.at(iN), &_xNodes.at(_idxNN.at(iN)*_hdim), 1,
                             _xClosest.data(), 1);
            _xSeeds.insert( _xSeeds.end(), _xClosest.begin(), _xClosest.end() );
            iter += _ldim;
        }
        ++counter;
    }
    _xSeeds.shrink_to_fit();
    int nHSeeds = _xSeeds.size()/_hdim;
    cout<<nHSeeds<<" [out of "<<nLSeeds<<"] seed points initialized correctly."<<endl;


    PointSetIO io;
    io.WritePointSet("x_seeds.txt", _xSeeds, nHSeeds, _hdim);

    /* Initialize the high-dimneional searcher */
    _highDimSearch = make_shared<ANNSearch>(1, _hdim, _xSeeds);


//    /* Define the low-dimensional bounds*/
//    _lowerBound.resize(_ldim, +numeric_limits<double>::max());
//    _upperBound.resize(_ldim, -numeric_limits<double>::max());
//    size_t iD;
//    for ( int iN = 0; iN < nHSeeds; ++iN )
//        for ( iD = 0; iD < _ldim; ++iD )
//        {
//            _lowerBound.at(iD) = min( _xiSeeds.at(iN*_ldim+iD), _lowerBound.at(iD) );
//            _upperBound.at(iD) = max( _xiSeeds.at(iN*_ldim+iD), _upperBound.at(iD) );
//        }
//    for ( iD = 0; iD < _ldim; ++iD )
//    {
//        _lowerBound.at(iD) += sqrt(_epsilon);
//        _upperBound.at(iD) -= sqrt(_epsilon);
//    }
//    cout<<"lower bound = [";
//    for( auto iter : _lowerBound )
//        cout<<iter<<"  ";
//    cout<<"]"<<endl;
//    cout<<"upper bound = [";
//    for( auto iter : _upperBound )
//        cout<<iter<<"  ";
//    cout<<"]"<<endl;

    if (_bTrustRegion  ) cout<<"TRUST REGION is set."<<endl;
    if (_bNewtonRaphson) cout<<"NEWTON RAPHSON is set."<<endl;


    return 0;
}



bool ClosestPointProjection::CheckSettings( )
{
    bool bError = false;

    if (_ldim < 1)
    {
        cerr<<__func__<<" >>  Dimension of the low-dimensional space is less than 1("<<_ldim<<")."<<endl;
        bError = true;
    }
    if (_hdim < 1)
    {
        cerr<<__func__<<" >>  Dimension of the high-dimensional space is less than 1 ("<<_hdim<<")."<<endl;
        bError = true;
    }
    if (_nNodes < 1)
    {
        cerr<<__func__<<" >>  Number of nodes either is not defined or is less than 1 ( "<<_nNodes<<" )."<<endl;
        bError = true;
    }
    if ((!_bNewtonRaphson && !_bTrustRegion))
    {
        cerr<<__func__<<" >> Either Newton-Raphson or Trust Region has to be set."<<endl;
        bError = true;
    }

    return bError;
}



const vector<double>& ClosestPointProjection::FindClosestPoint(const vector<double> &xAligned, bool bUseSeed)
throw (runtime_error)
{
    _bSeed = bUseSeed;
    _xAligned = xAligned;
    if ( this->FindClosestPoint() )
        throw runtime_error("Failed to find the closest points!");
    return _xiVal;
}



const vector<double>& ClosestPointProjection::FindClosestPoint(const vector<double> &xAligned)
throw (runtime_error)
{
    _xAligned = xAligned;
    if ( this->FindClosestPoint() )
        throw runtime_error("Failed to find the closest points!");
    return _xiVal;
}



const vector<double>& ClosestPointProjection::FindClosestPoint(vector<double>::iterator iter_begin,
                                             vector<double>::iterator iter_end)
throw (runtime_error)
{
    if (distance(iter_begin,iter_end) != _hdim)
    {
        cerr<<__func__<<"The dimension of input sample points iterator is not correct."<<endl;
        return _xiVal;
    }

    copy( iter_begin, iter_end, _xAligned.begin() );
    if ( this->FindClosestPoint() )
        throw runtime_error("Not Found!");
    return _xiVal;
}


int ClosestPointProjection::FindClosestPoint()
{
    if (_bSeed) this->UseSeed();

    int ierr = 0;
    /* Perfotm a least square fit with Trust Region or Newtom-Raphson algorithm*/
    if ( _bTrustRegion )
    {
        ierr = this->MinimizeWithTrustRegion();
        if(ierr && !_bSeed)
        {
            this->UseSeed();
            ierr = this->MinimizeWithTrustRegion();
        }
        if ( !_bSeed )
        {
            transform(_xiVal.begin(), _xiVal.end(), _xiValOld.begin(), _dxi.begin(), std::minus<double>() );
            _norm_dxi = cblas_dnrm2(_ldim, _dxi.data(), 1);
            if(_norm_dxi > 2*_spacing )
            {
                ierr = 1;
                cout<<"A jump is very probable: norm dxi = "<<_norm_dxi<<endl;
            }
        }
    }
    if ( _bNewtonRaphson )
    {
        ierr = this->MinimizeWithNewtonRaphson(false);
        if (ierr)
        {
            cout<<"Warning: Newton-Raphson failed. Line search..."<<endl;
            _xiVal = _xiInit;
            ierr = this->MinimizeWithNewtonRaphson(true);
            if (ierr && !_bSeed)
            {
//                if (_bTrustRegion)
//                {
//                    cout<<"Warning: Newton-Raphson with line search failed. Trust region..."<<endl;
//                    //                    cerr<<"Warning: Newton-Raphson with line search failed. Trust region..."<<endl;
//                    _xiVal = _xiInit;
//                    ierr = this->MinimizeWithTrustRegion();
//                }
//                else if(!_bSeed)
                cout<<"Warning: Newton-Raphson with line search failed. Next try with high-dim seed..."<<endl;
                //                    cerr<<"Warning: Newton-Raphson with line search failed. Next try with high-dim seed..."<<endl;
                this->UseSeed();
                ierr = this->MinimizeWithNewtonRaphson(true);
                if (ierr)
                    cerr<<"ERROR: Newton-Raphson with High-Dim SEED failed."<<endl;
            }
        }
    }

    if( ierr )
    {
        _distVal = cblas_dnrm2(_hdim, _xDiff.data(), 1);

        cerr<<"Failiure: xiVal = [";
        for(auto iter:_xiVal)
            cerr<<iter<<" ";
        cerr<<"] , xiInit = [";
        for(auto iter:_xiInit)
            cerr<<iter<<" ";
        cerr<<"] , distance = "<<_distVal<<endl;

        _xiVal = _xiInit;
        cerr<<__func__<<" >> The initial value is taken as the solution."<<endl;

        return ierr;
    }
    else
    {
        _bSeed    = false;
        _xiValOld = _xiVal;
        _xiInit   = _xiVal;

        return this->CalculateJacobian();
    }
}



void ClosestPointProjection::UseSeed()
{
    /* Search in high-dimensional space to find nearest seed point */
    cout<<"_xAligned size = "<<_xAligned.size()<<endl;
    _idxSeed = _highDimSearch->FindNearestNeighbors( _xAligned ).at(0);
    PointSetIO ios;
    vector<double> xSeed(_xSeeds.begin()+_idxSeed*_hdim, _xSeeds.begin()+(_idxSeed+1)*_hdim);
    cout<<"xSeed size = "<<xSeed.size()<<endl;

    ios.WritePointSet("align_"+to_string(_idxSeed)+".txt", _xAligned, 1, _hdim);
    ios.WritePointSet("seed_" +to_string(_idxSeed)+".txt", xSeed    , 1, _hdim);
    seedIter = _xiSeeds.begin()+_idxSeed*_ldim;
    copy( seedIter, seedIter+_ldim, _xiInit.begin() );
    _xiVal = _xiInit;
    ++_nTakenSeeds;
    _bSeed = true;
}



// Local max-entropy shape functions are computed by using Newton-Raphson method
// with possible Brent's line search algorithm
int ClosestPointProjection::MinimizeWithNewtonRaphson( bool applyLS )
{
    int    iN, iD, jD;
    int    ierr      = 0;
    double *xNode    = nullptr;

    /* calculate shape functions at \xi */
    if (_parametrization.ComputeLme( _xiVal, 2 ))
    {
        cerr<<__func__<<"LME failed!"<<endl;
        return 1;
    }
    _idxNN = _parametrization.GetNearestNeighbors();
    _knn   = _idxNN.size();
    //        cout<<"knn = "<<_knn<<endl;

    cout<<setprecision(16);
    for ( int iter = 0; iter < _maxIterNR; ++iter )
    {
        //        cout<<"xiVal = [";
        //        for(auto  vIter: _xiVal)
        //            cout<<vIter<<" ";
        //        cout<<"]"<<endl;

        _pa    = _parametrization.GetShapeFunctions();
        _dpa   = _parametrization.GetShapeFunctionsGradient();
        _hpa   = _parametrization.GetShapeFunctionsHessian();

        /* Initialized vectors to zero */
        fill(_xClosest .begin(), _xClosest .end(), 0.0);
        fill(_dxClosest.begin(), _dxClosest.end(), 0.0);
        fill(_hxClosest.begin(), _hxClosest.end(), 0.0);

        /* Calculate x, dx, hx */
        for ( iN = 0; iN < _knn; ++iN )
        {
            xNode = &_xNodes.at(_idxNN.at(iN)*_hdim);
            cblas_daxpy( _hdim, _pa.at(iN), xNode, 1, _xClosest.data(), 1);
            for ( iD = 0; iD < _ldim ; ++iD )
                cblas_daxpy( _hdim, _dpa.at(iN*_ldim +iD), xNode, 1, &_dxClosest.at(iD*_hdim), 1);
            for ( iD = 0; iD < _ldim2; ++iD )
                cblas_daxpy( _hdim, _hpa.at(iN*_ldim2+iD), xNode, 1, &_hxClosest.at(iD*_hdim), 1);
        }

        /* diff = x - xs */
        transform( _xClosest.begin() , _xClosest.end(), _xAligned.begin(),
                   _xDiff.begin(), std::minus<double>() );
        /* dDiff = dx.xDiff */
        cblas_dgemv( CblasRowMajor, CblasNoTrans, _ldim, _hdim, 1.0,
                     _dxClosest.data(), _hdim, _xDiff.data(), 1, 0.0 , _dDiff.data(), 1 );
        /* hDiff = hx.xDiff + dxdx */
        cblas_dgemv( CblasRowMajor, CblasNoTrans, _ldim2, _hdim, 1.0,
                     _hxClosest.data(), _hdim, _xDiff.data(), 1, 0.0 , _hDiff.data(), 1 );
        for ( iD = 0; iD < _ldim; ++iD )
            for ( jD = 0; jD < _ldim; ++jD )
                _hDiff.at(iD*_ldim+jD) += cblas_ddot(_hdim, &_dxClosest.at(iD*_hdim), 1, &_dxClosest.at(jD*_hdim), 1 );

        _norm_res = cblas_dnrm2(_ldim, _dDiff.data(), 1 );
//        cout<<"normRes = "<<_norm_res<<endl;
//        cout<<"fvalue = "<<0.5*cblas_ddot(_hdim, _xDiff.data(), 1, _xDiff.data(), 1 )<<endl;
//        cout<<"dist   = "<<cblas_dnrm2(_hdim, _xDiff.data(), 1)<<endl;

//        double det = _hDiff.at(0)*_hDiff.at(3) - _hDiff.at(1)*_hDiff.at(1) ;
//        cout<<"det = "<<det<<endl;

        if ( std::isnan( _norm_res ) )
        {
            cerr<<__func__<<" >> "<< (applyLS?"with":"without") <<" line search)"
               <<" Norm of res vector is NaN."<<endl;
            cout<<__func__<<" >> "<< (applyLS?"with":"without") <<" line search)"
               <<" Norm of res vector is NaN."<<endl;
            return 3;
        }
        else if ( _norm_res < _tolNR )
        {
            // The solution is converged.
//            cout<<"iteration = "<<iter<<endl;
            return 0;
        }
        else // Solve the system of hDiff*dxi = dDiff
        {
            // Compote Cholesky factorization of hDiff = L*L' ( hDiff is overwritten by L )
            ierr = LAPACKE_dpotrf( LAPACK_COL_MAJOR, 'L', _ldim, _hDiff.data(), _ldim );
            if (ierr == 0)
            {
                // Solve L*L'*dxi = dDiff
                ierr = LAPACKE_dpotrs( LAPACK_COL_MAJOR, 'L', _ldim, 1, _hDiff.data(), _ldim, _dDiff.data(), _ldim );
                if (ierr == 0)
                {
                    // ***update dxi***
                    _dxi = _dDiff;
                }
                else if (ierr < 0)
                {
                    cerr<<__func__<<" >> The "<<ierr<<"-th parameter of LAPACKE_dpotrs had an illegal value."<<endl;
                    cout<<__func__<<" >> The "<<ierr<<"-th parameter of LAPACKE_dpotrs had an illegal value."<<endl;
                    return 1;
                }
                else if (ierr > 0)
                {
                    cerr<<__func__<<" >> Unexpexted error in LAPACKE_dpotrs ( info ="<<ierr<<")."<<endl;
                    cout<<__func__<<" >> Unexpexted error in LAPACKE_dpotrs ( info ="<<ierr<<")."<<endl;
                    return 1;
                }
            }
            else if (ierr < 0)
            {
                cerr<<__func__<<" >> The "<<ierr<<"-th parameter of LAPACKE_dpotrf had an illegal value."<<endl;
                cout<<__func__<<" >> The "<<ierr<<"-th parameter of LAPACKE_dpotrf had an illegal value."<<endl;
                return 1;
            }
            else if (ierr > 0)
            {
                cerr<<__func__<<" >> The Hessian matrix (the leading minor of order "<<ierr<<") is not positive definite."<<endl;
                cout<<__func__<<" >> The Hessian matrix (the leading minor of order "<<ierr<<") is not positive definite."<<endl;
                return 1;
            }
        }

        _alpha = 1.0;  // Withuot Line Search
        if ( applyLS ) // Apply Line Search
                if( this->BrentLineSearch( _maxIterLS, _tolLS, _alpha ) )
                    cout<<"Warning: Line Search did not converged!"<<endl;
//        cout<<"alpha = "<< _alpha <<endl;

        // Update xiVal = xiVal - alpha*dxi
        cblas_daxpy(_ldim, -_alpha, _dxi.data(), 1, _xiVal.data(), 1);

        /* Find if there is a jump in the trajectory */
        transform(_xiInit.begin(), _xiInit.end(), _xiVal.begin(), _dxi.begin(), std::minus<double>() );
        _norm_dxi = cblas_dnrm2(_ldim, _dxi.data(), 1);
        if( _norm_dxi > 2*_spacing )
        {
            cout<<"A jump is very probable: norm dxi = "<<_norm_dxi<<endl;
            cerr<<"A jump is very probable: norm dxi = "<<_norm_dxi<<endl;
            return 4;
        }
        else
        {
            /* calculate shape functions at \xi */
            if (_parametrization.ComputeLme( _xiVal, _idxNN, 2 ))
            {
                cerr<<__func__<<"LME failed!"<<endl;
                return 1;
            }

        }
    }

    cerr<<__func__<<" >> Maximum number of iterations "<< (applyLS?"with":"without")
       <<" line search reached ("<<_maxIterNR<<")."<<endl;
    cout<<__func__<<" >> Maximum number of iterations "<< (applyLS?"with":"without")
       <<" line search reached ("<<_maxIterNR<<")."<<endl;
    return 2;
}


int ClosestPointProjection::BrentLineSearch( int maxIter, double TolX, double &alpha )
{
    //Initialize variables
    double a   = -0.01;
    double b   = 5.00;
    double xm  = 0.0;

    double u   = 0.0;
    double v   = a + _gold*(b-a);
    double w   = v;
    double x   = v;

    double e   = 0.0;  // The distance moved on the step before last
    double d  = 0.0;

    double fx  = this->LineSearchFun( x );
    double fv  = this->LineSearchFun( v );
    double fu  = 0.0;
    double fw  = fx;

    double p, q, r, tol1, tol2;

    // Start iteration until convergence!
    for ( int iter = 0; iter < maxIter; ++iter )
    {
        // Update xm, tol1, tol2
        xm    = 0.5*(a+b);
        tol1  = fabs(x)*TolX + _sqrtEpsilon;
        tol2  = 2*tol1;

        // Check the convergence criterion
        if ( fabs(x-xm) <= tol2-0.5*(b-a) )
        {
            alpha = x;
            return 0;
        }
        else
        {
            // Is a parabolic fit possible?
            if ( fabs(e) > tol1)
            {
                // Yes, fit a parabola
                r  = (x-w)*(fx-fv);   //trial parabolic fit
                q  = (x-v)*(fx-fw);
                p  = (x-v)*q-(x-w)*r;
                q  = 2.0*(q-r);

                if (q > 0.0) p = -p;
                q = fabs(q);
                r = e; // eTemp;
                e = d;

                // Is the parabola acceptable?
                if ( fabs(p) < fabs((0.5*q*r)) &&  p > q*(a-x)  &&  p < q*(b-x) )
                {
                    // Yes, take a parabolic interpolation step
                    d = p/q;
                    u = x+d;

                    // function must not be evaluated too close to a or b
                    if ( u-a < tol2 || b-u < tol2 )
                        d = this->sign( tol1, xm-x );
                }
                else
                {
                    // Not acceptable, take a golden section step
                    e = (x >= xm) ? a-x : b-x ;
                    d = _gold* e;
                }
            }
            else
            {
                // Not possible, take a golden section step
                e = (x >= xm) ? a-x : b-x ;
                d = _gold* e;
            }

            //F must not be evaluated too close to x
            u  = (fabs(d) >= tol1) ?  x+d : x+this->sign(tol1,d);
            fu = this->LineSearchFun( u ); //The only function evaluation per iteration

            // Update a,b,v,w, and x
            if (fu <= fx)
            {
                if (u >= x)
                    a = x;
                else
                    b = x;
                v  = w;
                w  = x;
                x  = u;
                fv = fw;
                fw = fx;
                fx = fu;
            }
            else // fu > fx
            {
                if ( u < x )
                    a = u;
                else
                    b = u;

                if ( fu <= fw || w == x )
                {
                    v  = w;
                    w  = u;
                    fv = fw;
                    fw = fu;
                }
                else if ( fu <= fv || v == x || v == w )
                {
                    v  = u;
                    fv = fu;
                }
            }
        }
    }
    return 1;
}



double ClosestPointProjection::LineSearchFun(double alpha)
{

    int iN;
    double *xNode = nullptr;
    // lam_ls = lam + alpha*dlam
    _xiLS = _xiVal;
    cblas_daxpy( _ldim, -alpha, _dxi.data(), 1, _xiLS.data(), 1 );

    /* calculate shape functions at \xi */
    if (_parametrization.ComputeLme( _xiLS, _idxNN, 0 ))
    {
        cerr<<__func__<<"LME failed!"<<endl;
        return numeric_limits<double>::max();
    }
    _pa    = _parametrization.GetShapeFunctions();
//    _idxNN = _parametrization.GetNearestNeighbors();
//    _knn   = _idxNN.size();

    /* Initialized vectors to zero */
    fill(_xClosest .begin(), _xClosest .end(), 0.0);

    /* Calculate x */
    for ( iN = 0; iN < _knn; ++iN )
    {
        xNode = &_xNodes.at(_idxNN.at(iN)*_hdim);
        cblas_daxpy( _hdim, _pa.at(iN), xNode, 1, _xClosest.data(), 1);
    }

    /* diff = x - xs */
    transform( _xClosest.begin() , _xClosest.end(), _xAligned.begin(),
               _xDiff.begin(), std::minus<double>() );

    /* fval = 0.5*diff^2 */
    return 0.5*cblas_ddot(_hdim, _xDiff.data(), 1, _xDiff.data(), 1 );
}



int ClosestPointProjection::MinimizeWithTrustRegion( )
{

    bool   bSuccessful = false; // controls of rci cycle
    int    RCI_Request = 0;     // reverse communication interface (RCI) parameter
    double rs          = 0.0;       // initial step bound

    fill(_info.begin(),_info.end(), 0.0); // results of input parameter checking

    CalculateCostValue();    // fvec = M(xi)-xs
    CalculateCostGradient(); // fjac = DM(xi)

    /* initialize solver (allocate mamory, set initial values) */
    if (dtrnlsp_init( &_handle, &_ldim, &_hdim, _xiVal.data(), _eps.data(), &_maxIterTR, &_maxTrialIter, &rs )
//    if (dtrnlspbc_init( &_handle, &_ldim, &_hdim, _xiVal.data(), _lowerBound.data(), _upperBound.data(), _eps.data(), &_maxIter, &_maxTrialIter, &rs )
            != TR_SUCCESS)
    {
        cerr<<__func__<<" >> error in dtrnlspbc_init."<<endl;
        return 1;
    }

    /* Checks the correctness of handle and arrays containing Jacobian matrix,
       objective function, lower and upper bounds, and stopping criteria. */
    if (dtrnlsp_check( &_handle, &_ldim, &_hdim, _fjac.data(), _fvec.data(), _eps.data(), _info.data() )
//    if (dtrnlspbc_check( &_handle, &_ldim, &_hdim, _fjac.data(), _fvec.data(), _lowerBound.data(), _upperBound.data(), _eps.data(), _info.data() )
            != TR_SUCCESS)
    {
        cerr<<__func__<<" >> error in dtrnlspbc_check."<<endl;
        return 1;
    }
    else
    {
        if (_info.at(0) != 0 || // The handle is not valid.
            _info.at(1) != 0 || // The fjac array is not valid.
            _info.at(2) != 0 || // The fvec array is not valid.
            _info.at(3) != 0 || // The LW array is not valid.
            _info.at(4) != 0 || // The UP array is not valid.
            _info.at(5) != 0    // The eps array is not valid.
           )
        {
            cerr<<__func__<<" >> Input parameters for dtrnlspbc_solve are not valid."<<endl;
            for ( int i = 0; i < 6; ++i )
                if (_info.at(i) != 0)
                    cerr<<"info["<<i<<"] = "<<_info.at(i)<<endl;
            return 1;
        }
    }

    int myCounter = 0;
    /* RCI cycle */
    while (!bSuccessful)
    {
        if (dtrnlsp_solve( &_handle, _fvec.data(), _fjac.data(), &RCI_Request ) != TR_SUCCESS)
//        if (dtrnlspbc_solve( &_handle, _fvec.data(), _fjac.data(), &RCI_Request ) != TR_SUCCESS)
        {
            cerr<<__func__<<" >> error in dtrnlspbc_solve."<<endl;
            return 1;
        }
        /* according with rci_request value, we do next step */
        if (RCI_Request == -1 || RCI_Request == -2 || RCI_Request == -3 ||
            RCI_Request == -4 || RCI_Request == -5 || RCI_Request == -6)
        {
//            cout<<"RCI_Request = "<<RCI_Request<<endl;
            bSuccessful = true; // exit rci cycle succesfully
        }
        else if (RCI_Request == 1)
        {
            /* recalculate function value */
            if ( this->CalculateCostValue() )
            {
                cerr<<"Calculating fvec has been failed."<<endl;
                fill( _fvec.begin(), _fvec.end(), numeric_limits<double>::max() );
            }
        }
        else if (RCI_Request == 2)
        {
            /* Try to recalculate Jacobian matrix with LME */
            if ( this->CalculateCostGradient() )
            {
                cerr<<"Calculating fjac has been failed."<<endl;
//                fill( _fjac.begin(), _fjac.end(), numeric_limits<double>::max() );
//                return 1;
            }
        }
        ++myCounter;
    }


    int ierr = 0;
    /* get solution statuses */
    if (dtrnlsp_get( &_handle, &_nIters, &_stopCriteria, &_initRes, &_finalRes ) != TR_SUCCESS)
//    if (dtrnlspbc_get( &_handle, &_nIters, &_stopCriteria, &_initRes, &_finalRes ) != TR_SUCCESS)
    {
        cerr<<__func__<<" >> error in dtrnlspbc_get"<<endl;
        ierr = 1;
    }
    if (dtrnlsp_delete (&_handle) != TR_SUCCESS)
//    if (dtrnlspbc_delete (&_handle) != TR_SUCCESS)
        cerr<<__func__<<" >> error in dtrnlspbc_delete"<<endl;
    mkl_free_buffers ();

//    cout<<"nIters       = "<<_nIters<<endl;
//    cout<<"Counter      = "<< myCounter<<endl;
//    cout<<"stopCriteria = "<<_stopCriteria<<endl;
//    cout<<"initRes      = "<<_initRes<<endl;
//    cout<<"finalRes     = "<<_finalRes<<endl;
    /* norm of cost function has to be less than the tolerence */
    if (_stopCriteria == 1 || _stopCriteria == 4)
        ierr = 1;


    return ierr;
}



int ClosestPointProjection::CalculateCostValue( )
{

    if (_parametrization.ComputeLme( _xiVal, 0 ))
    {
        cerr<<__func__<<" >> error in calculating shape functions."<<endl;
        cerr<<__func__<<" >> xi_requested = "<<_xiVal.at(0)<<","<<_xiVal.at(1)<<endl;
        return 1;
    }
    else
    {
        _pa    = _parametrization.GetShapeFunctions();
        _idxNN = _parametrization.GetNearestNeighbors();
        _knn   = _idxNN.size();
        fill(_xClosest.begin(), _xClosest.end(), 0.0);
        for ( int iN = 0; iN < _knn; ++iN )
            cblas_daxpy( _hdim, _pa.at(iN), &_xNodes.at(_idxNN.at(iN)*_hdim), 1,
                         _xClosest.data(), 1);
        transform( _xClosest.begin() , _xClosest.end(), _xAligned.begin(),
                   _fvec.begin(), std::minus<double>() );
        return 0;
    }
}



int ClosestPointProjection::CalculateCostGradient( )
{
    if (_parametrization.ComputeLme( _xiVal, 1 ))
    {
        cerr<<__func__<<" >> error in calculating gradient of shape functions."<<endl;
        cerr<<__func__<<" >> xi_requested = "<<_xiVal.at(0)<<","<<_xiVal.at(1)<<endl;
        return 1;
    }
    else
    {
        int iN, iD;
        double *xNode;
        _dpa   = _parametrization.GetShapeFunctionsGradient();
        _idxNN = _parametrization.GetNearestNeighbors();
        _knn   = _idxNN.size();
        fill(_fjac.begin(), _fjac.end(), 0.0);
        for ( iN = 0; iN < _knn; ++iN )
        {
            xNode = &_xNodes.at(_idxNN.at(iN)*_hdim);
            /* To be consistent with mkl solver,_fJac [HxL] is transpose of dxClosest */
            for ( iD = 0; iD < _ldim; ++iD )
                cblas_daxpy( _hdim, _dpa.at(iN*_ldim+iD), xNode, 1, &_fjac.at(iD*_hdim), 1 );
        }
        return 0;
    }
}



int ClosestPointProjection::CalculateJacobian()
{
    int     iN, iD;
    double *xNode;

    /* calculate shape functions at \xi */
    if (_parametrization.ComputeLme( _xiVal, _idxNN, 2 ))
    {
        cerr<<__func__<<"LME failed!"<<endl;
        return 1;
    }
    _pa    = _parametrization.GetShapeFunctions();
    _dpa   = _parametrization.GetShapeFunctionsGradient();
    _hpa   = _parametrization.GetShapeFunctionsHessian();
//    _idxNN = _parametrization.GetNearestNeighbors();
//    _knn   = _idxNN.size();

    /* Initialized vectors to zero */
    fill(_xClosest .begin(), _xClosest .end(), 0.0);
    fill(_dxClosest.begin(), _dxClosest.end(), 0.0);
    fill(_hxClosest.begin(), _hxClosest.end(), 0.0);

    /* Calculate x, dx, hx */
    for ( iN = 0; iN < _knn; ++iN )
    {
        xNode = &_xNodes.at(_idxNN.at(iN)*_hdim);
        cblas_daxpy( _hdim, _pa.at(iN), xNode, 1, _xClosest.data(), 1);
        for ( iD = 0; iD < _ldim ; ++iD )
            cblas_daxpy( _hdim, _dpa.at(iN*_ldim +iD), xNode, 1, &_dxClosest.at(iD*_hdim), 1);
        for ( iD = 0; iD < _ldim2; ++iD )
            cblas_daxpy( _hdim, _hpa.at(iN*_ldim2+iD), xNode, 1, &_hxClosest.at(iD*_hdim), 1);
    }

    /* diff = x - xs */
    transform( _xClosest.begin() , _xClosest.end(), _xAligned.begin(),
               _xDiff.begin(), std::minus<double>() );

    /* Calculate Jacobian by considering the curvature of manifold */
    return (this->*pCorrectJacobian)();
}



int ClosestPointProjection::CorrectJacobian1D()
{
    /* cor =  dx*dx' + hx*diff */
    _cor.at(0) = inner_product(_dxClosest.begin(), _dxClosest.end(), _dxClosest.begin(), 0.0)
               + inner_product(_xDiff    .begin(), _xDiff    .end(), _hxClosest.begin(), 0.0);

    if (_cor.at(0) < _epsilon)
    {
        cerr<<__func__<<" >> Curvature correction matrix is singular!"<<endl;
        return 1;
    }
    else
    {
        /* Jac = inv(cor) * dx */
        _xiJac = _dxClosest;
        cblas_dscal(_hdim, 1.0/_cor.at(0), _xiJac.data(), 1);
        return 0;
    }
}



int ClosestPointProjection::CorrectJacobian2D()
{
    /* cor = dx*dx' + hx*diff */
    _cor.at(0) = inner_product(_dxClosest.begin(), _dxClosest.begin()+_hdim, _dxClosest.begin()      , 0.0)
               + inner_product(_xDiff.begin()    , _xDiff.end()            , _hxClosest.begin()      , 0.0);

    _cor.at(1) = inner_product(_dxClosest.begin(), _dxClosest.begin()+_hdim, _dxClosest.begin()+_hdim, 0.0)
               + inner_product(_xDiff.begin()    , _xDiff.end()            , _hxClosest.begin()+_hdim, 0.0);

    _cor.at(2) = _cor.at(1);

    _cor.at(3) = inner_product(_dxClosest.begin()+_hdim, _dxClosest.end()  , _dxClosest.begin()+_hdim, 0.0)
               + inner_product(_xDiff.begin(), _xDiff.end(), _hxClosest.begin()+3*_hdim, 0.0);

    double det = _cor.at(0)*_cor.at(3) - _cor.at(1)*_cor.at(2);
    if (fabs(det) < _epsilon)
    {
        cerr<<__func__<<" >> error in calculation of inv(C*T), det = "<<det<<"."<<endl;
        return 1;
    }
    else
    {
        _xiJac = _dxClosest;
        cblas_daxpby(_hdim, -_cor.at(1)/det, &_dxClosest.at(_hdim), 1, _cor.at(3)/det, &_xiJac.at(0)    , 1 );
        cblas_daxpby(_hdim, -_cor.at(2)/det, &_dxClosest.at(0)    , 1, _cor.at(0)/det, &_xiJac.at(_hdim), 1 );
        return 0;
    }
}



int ClosestPointProjection::CorrectJacobianND()
{
    /* dxdx = dx*dx' */
    cblas_dgemm( CblasRowMajor, CblasNoTrans, CblasTrans, _ldim, _ldim, _hdim,
                 1.0, _dxClosest.data(), _hdim, _dxClosest.data(), _hdim,
                 0.0 , _dxdx.data(), _ldim );

    /* corr = hx*fval */
    cblas_dgemv( CblasRowMajor, CblasNoTrans, _ldim2, _hdim, 1.0,
                 _hxClosest.data(), _hdim, _xDiff.data(), 1, 0.0 , _cor.data(), 1 );

    /* corr = dxdx + hx*fval */
    cblas_daxpy(_ldim2, 1.0, _dxdx.data(), 1, _cor.data(), 1);

    /* invCor = inv(cor) */
    int ierr = LAPACKE_dpotrf( LAPACK_ROW_MAJOR, 'L', _ldim, _cor.data(), _ldim );
    if (ierr < 0)
    {
        cerr<<__func__<<" >> The "<<ierr<<"-th parameter of LAPACKE_dpotrf had an illegal value."<<endl;
        return 1;
    }
    else if (ierr > 0)
    {
        cerr<<__func__<<" >> The J matrix (the leading minor of order "<<ierr<<") is not positive definite."<<endl;
        return 1;
    }
    else
    {
        /* xiJac = invCor*dx */
        _xiJac = _dxClosest;
        ierr = LAPACKE_dpotrs( LAPACK_ROW_MAJOR, 'L', _ldim, _hdim, _cor.data(), _ldim, _xiJac.data(), _hdim );
        if (ierr < 0)
        {
            cerr<<__func__<<" >> The "<<ierr<<"-th parameter of LAPACKE_dpotrs had an illegal value."<<endl;
            return 1;
        }
        else if (ierr > 0)
        {
            cerr<<__func__<<" >> Unexpexted error in LAPACKE_dpotrs ( info ="<<ierr<<")."<<endl;
            return 1;
        }
        else
            return 0;
    }
}

double ClosestPointProjection::CalculateDistance()
{
    _distVal = cblas_dnrm2(_hdim, _xDiff.data(), 1);

    /* distJac = ( diff (dx*xiJac - I) ) / dist */
    cblas_dgemm( CblasRowMajor, CblasTrans, CblasNoTrans, _hdim, _hdim, _ldim,
                 1.0, _dxClosest.data(), _hdim, _xiJac.data(), _hdim,
                 0.0, _tempJac.data(), _hdim );
    for ( int iD = 0; iD < _hdim; ++iD )
        _tempJac.at(iD*_hdim+iD) -= 1.0;

    cblas_dgemv(CblasRowMajor, CblasNoTrans, _hdim, _hdim, 1.0, _tempJac.data(), _hdim,
                _xDiff.data(), 1, 0.0, _distJac.data(), 1 );

    cblas_dscal(_hdim, 1.0/_distVal, _distJac.data(), 1);

    return _distVal;
}



