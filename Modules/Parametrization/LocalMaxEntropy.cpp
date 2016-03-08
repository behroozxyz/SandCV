/*=========================================================================
 * File         : LocalMaxEntropy.cpp
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
#include "LocalMaxEntropy.h"
static const double _epsilon = numeric_limits<double>::epsilon(); ///< Machine precision

LocalMaxEntropy::LocalMaxEntropy(int dim, const vector<double> &xNodes , double spacing, double tol0, double gamma) :
    _dim(dim),
    _h(spacing),
    _gamma(gamma),
    _tol0(tol0),
    _xNodes(xNodes),
    _nNodes(xNodes.size()/dim)
{
    this->Update();
}


LocalMaxEntropy::LocalMaxEntropy(int dim, int nNodes, const vector<double>& xNodes , double spacing, double tol0, double gamma) :
    _dim(dim),
    _h(spacing),
    _gamma(gamma),
    _tol0(tol0),
    _nNodes(nNodes)
{
    _xNodes.assign( xNodes.begin(), xNodes.begin()+_nNodes*_dim );
    this->Update();
}


LocalMaxEntropy::LocalMaxEntropy(int dim, int nNodes, const double *xNodes , double spacing, double tol0, double gamma) :
    _dim(dim),
    _h(spacing),
    _gamma(gamma),
    _tol0(tol0),
    _nNodes(nNodes)
{
    _xNodes.assign( xNodes, xNodes+_nNodes*_dim );
    this->Update();
}


LocalMaxEntropy::LocalMaxEntropy(const LocalMaxEntropy &inLME) :
    _dim(inLME._dim),
    _h(inLME._h),
    _gamma(inLME._gamma),
    _tol0(inLME._tol0),
    _xNodes(inLME._xNodes),
    _nNodes(inLME._nNodes)
{
    this->Update();
}


int LocalMaxEntropy::Update()
throw( exception )
{
    if ( this->CheckSettings() )
    {
        cerr<<__func__<<" >> CheckSettings went wrong."<<endl;
        throw exception();
    }
    _sqrtEpsilon = sqrt( _epsilon );      // square root of epsilon
    _gold        = 0.5 * (3.0-sqrt(5.0)); // Golden Ratio: 0.3819660

    _dim2 = _dim*_dim;
    _beta = _gamma/_h/_h;
//    _range = sqrt(-log(_tol0)/_beta);
    _knn  = ceil( pow( 2.0*sqrt(-log(_tol0)/_gamma), _dim ) );
    cout<<"knn = "<<_knn<<endl;
//    _knn  = 93; //FIXME: Temporary, only for checking
    _knn  = _knn < _nNodes ? _knn : _nNodes;

    _xSample .resize( _dim );
    if (_dim != 1)
    {
        // the vectors that depends only on dim
        _lam     .resize( _dim );
        _dlam    .resize( _dim );
        _lamLS   .resize( _dim );
        _res     .resize( _dim );
        _R       .resize( _dim2 );
        _Jchol   .resize( _dim2 );
    }
    // the vectors that depends also on knn
    _lam_dx  .resize( _knn );
    _beta_dx2.resize( _knn );
    _dx      .resize( _knn*_dim );
    _Jidx    .resize( _knn*_dim );
    _dxdx    .resize( _knn*_dim2 );
    _M       .resize( _knn*_dim2 );
    // The vectors containing the shape functions and its dervatives
    _pa .resize( _knn );
    _dpa.resize( _knn*_dim );
    _hpa.resize( _knn*_dim2 );

    // Nearesr Neighbor Searcher
    _searcher = make_shared<ANNSearch>(_knn, _dim, _xNodes);
//    _searcher = make_shared<ANNSearch>(_range, _dim, _xNodes);

    return 0;
}


bool LocalMaxEntropy::CheckSettings()
{
    bool bError = false;
    if (_dim < 1)
    {
        cerr<<__func__<<" >>  Dimension is "<<_dim<<" < 1 "<<endl;
        bError = true;
    }
    if (_nNodes < 1)
    {
        cerr<<__func__<<" >>  Number of nodes either is not defined or is less than 1 ( "<<_nNodes<<" )."<<endl;
        bError = true;
    }
    if (_xNodes.empty())
    {
        cerr<<__func__<<" >>  Data point vector is empty."<<endl;
        bError = true;
    }
    if (_xNodes.size() != _nNodes*_dim)
    {
        cerr<<__func__<<" >>  The size of the node points' vector is inconsistent with the given dimension."<<endl;
        bError = true;
    }
    if (_h <= _epsilon)
    {
        cerr<<__func__<<" >>  Sapcing ( "<<_h<<" ) is smaller than machine precision ("<<_epsilon<<")."<<endl;
        bError = true;
    }
    if ( _tol0 < _epsilon )
    {
        cerr<<__func__<<" >>  Cut-off tolerence ( Tol0 = "<<_tol0<<" ) is smaller than machine precision ("<<_epsilon<<")."<<endl;
        bError = true;
    }
    else if ( _tol0 > 1e-1 )
    {
        cerr<<__func__<<" >>  Cut-off tolerence ( Tol0 = "<<_tol0<<" ) is very large."<<endl;
        bError = true;
    }
    if ( _tolNR < _epsilon )
    {
        cerr<<__func__<<" >>  Newton-Raphson tolerence ( TolNR = "<<_tolNR<<" ) is smaller than machine precision ("<<_epsilon<<")."<<endl;
        bError = true;
    }
    if (_lowerLimit > _upperLimit)
    {
        cerr<<__func__<<" >>  Line search bounds are incorrectly set ( alphaMin = "<<_lowerLimit<<" alphaMax = "<<_upperLimit<<" )."<<endl;
        bError = true;
    }

    return bError;
}


int LocalMaxEntropy::ComputeLme( vector<double>::const_iterator iter_begin,
                                 vector<double>::const_iterator iter_end,
                                 int derivativeOrder)
{
    if (distance(iter_begin,iter_end) != _dim)
    {
        cerr<<__func__<<"The dimension of input sample points iterator is not correct."<<endl;
        return 1;
    }
    copy( iter_begin, iter_end, _xSample.begin() );
    return this->ComputeLme( _xSample, derivativeOrder );
}


int LocalMaxEntropy::ComputeLme( int dim, double *xSample, int derivativeOrder)
{
    if (dim != _dim)
    {
        cerr<<__func__<<"The dimension of input sample points is not correct."<<endl;
        return 1;
    }
    copy( &xSample[0], &xSample[dim], _xSample.begin() );
    return this->ComputeLme( _xSample, derivativeOrder );
}


// Compute the LME shape functions and the requested number of derivatives
int LocalMaxEntropy::ComputeLme(const vector<double> &xSample, int derivativeOrder)
{
    _idxNN = _searcher->FindNearestNeighbors( xSample );
//    _idxNN = _searcher->FindNearestNeighborsFR( xSample );

    bool withGradient= false;
    bool withHessian = false;
    switch ( derivativeOrder )
    {
    case 0 :
        break;
    case 1 :
        withGradient = true;
        break;
    case 2 :
        withGradient = true;  withHessian = true;
        break;
    default :
        cerr<<"The derivative order has been set incorrctly (="<<derivativeOrder<<
              "). Only first and second order derivatives are implemented.";
    }

    if (_dim == 1)
        return this->ComputeLme1D( xSample, _idxNN, withGradient, withHessian);
    else if (_dim == 2)
        return this->ComputeLme2D( xSample, _idxNN, withGradient, withHessian);
    else
        return this->ComputeLme( xSample, _idxNN, withGradient, withHessian);

}


// Compute the LME shape functions and the requested number of derivatives
// with predefined neighbor list
int LocalMaxEntropy::ComputeLme(const vector<double> &xSample,
                                const vector<int>    &idxNN,
                                int derivativeOrder)
{
    bool withGradient= false;
    bool withHessian = false;
    switch ( derivativeOrder )
    {
    case 0 :
        break;
    case 1 :
        withGradient = true;
        break;
    case 2 :
        withGradient = true;  withHessian = true;
        break;
    default :
        cerr<<"The derivative order has been set incorrctly (="<<derivativeOrder<<
              "). Only first and second order derivatives are implemented.";
    }

    if (_dim == 1)
        return this->ComputeLme1D( xSample, idxNN, withGradient, withHessian);
    else if (_dim == 2)
        return this->ComputeLme2D( xSample, idxNN, withGradient, withHessian);
    else
        return this->ComputeLme( xSample, idxNN, withGradient, withHessian);

}



// Compute only the LME shape functions and return them as a vector
// with predefined neighbor list
const vector<double>& LocalMaxEntropy::ComputeLme(const vector<double> &xSample,
                                                  const vector<int>    &idxNN)
{
    if (_dim == 1)
        this->ComputeLme1D( xSample, idxNN, false, false);
    else if (_dim == 2)
        this->ComputeLme2D( xSample, idxNN, false, false);
    else
        this->ComputeLme( xSample, idxNN, false, false);
    return _pa;
}

//______________________________________________________________________________
//------------------------------------------------------------------------------
// Local maximum entropy in one-dimensional space
//------------------------------------------------------------------------------

// Compute the LME Shape functions and its first and second spatial derivatives
int LocalMaxEntropy::ComputeLme1D(const vector<double> &xSample , const vector<int> &idxNN,
                                 bool withGradient, bool withHessian )
{
    // if the number of nearest neighbors are not the same as before,
    // resize all the correspondent vectors
    if ( _knn != idxNN.size())
    {
        _knn = idxNN.size();
        if (_knn < _dim)
        {
            cerr<<__func__<<" >>  Maximum number of nearest neighbors is less than the dimension."<<endl;
            return 1;
        }
        _beta_dx2.resize( _knn );
        _dx      .resize( _knn );
        _Jidx    .resize( _knn );
        _dxdx    .resize( _knn );
        _M       .resize( _knn );
        _pa      .resize( _knn );
        if (withGradient)
        {
            _dpa.resize( _knn );
            if (withHessian)
                _hpa.resize( _knn );
        }
    }
    // dx = x - x_a, loop over the neighbors of x_sample
    for ( int iN = 0; iN < _knn; ++iN )
        _dx.at(iN) = xSample.at(0) - _xNodes.at(idxNN.at(iN));

    // Shape functions
    _ierr = this->ComputeShapeFunctions1D( );

    if ( withGradient && _ierr == 0 )
    {
        // J = sum(pa*(dx dx)) - r*r
        _J1D = cblas_ddot( _knn, _pa.data(), 1, _dxdx.data(), 1) - _res1D*_res1D;

        // Ji_dx = invJ * (x-x_a)  [NxD]
        _Jidx = _dx; // initialize to dx
        cblas_dscal( _knn, 1.0/_J1D, _Jidx.data(), 1 );

        // Gradient of the shape functions
        this->ComputeGradients1D( );

        // Hessian of the shape functions
        if ( withHessian ) this->ComputeHessians1D( );
    }
    return _ierr;
}


// Calculate the LME shape functions in one dimensional space
int LocalMaxEntropy::ComputeShapeFunctions1D( )
{
    int    ierr = 0;
    double step = 0.0;

    // step= sum(-beta*dx2) and then norm of it andnormalize with norm(dx)
    for ( int iN = 0; iN < _knn; ++iN )
        step += _beta_dx2.at(iN) = -_beta * ( _dxdx.at(iN) = _dx.at(iN)*_dx.at(iN) );

    ierr = this->LmeNewton1D( false, 0.0 ); // withouth line search

    if ( ierr )
    {
        // step = norm(-beta*dx2)/norm(dx);
        step = fabs(step) / cblas_dnrm2(_knn, _dx.data(), 1);
        ierr = this->LmeNewton1D( true, step ); // with line search
    }

    return ierr;

}


// Local max-entropy shape functions are computed by using Newton-Raphson method
// with possible Brent's line search algorithm
int LocalMaxEntropy::LmeNewton1D ( bool applyLS, double step )
{
    int    iN;
    _lam1D = 0.0;

    for ( int iter = 0; iter < _maxIterNR; ++iter )
    {

        _Z = 0.0;
        // Z =sum ( exp(-beta*dx^2+lam.dx) )
        for ( iN = 0; iN < _knn; ++iN )
            _Z += _pa.at(iN) = exp( _beta_dx2.at(iN) + _lam1D*_dx.at(iN) );

        // Compute shape functions (normalize with Z)
        cblas_dscal( _knn, 1.0/_Z, _pa.data(), 1 );

        // residual vector: r = sum(pa*dx) [Dx1]
        _res1D = cblas_ddot( _knn, _pa.data(), 1, _dx.data(), 1);

        _norm_res = fabs(_res1D);

        if (std::isnan( _norm_res ))
        {
            cerr<<__func__<<" >> ("<< (applyLS?"with":"without") <<" line search)"
               <<" Norm of res vector is NaN."<<endl;
            return 3;
        }
        else if (_norm_res < _tolNR)
        {
            // The solution is converged.
            return 0;
        }
        else // Construct J and solve the system of J*dlam = -res
        {

            // J = sum(pa*(dx dx)) - r*r
            _J1D = cblas_ddot( _knn, _pa.data(), 1, _dxdx.data(), 1) - _res1D*_res1D;

            if ( fabs(_J1D) >= _epsilon )
            {
                // ***update dlam***
                _dlam1D = -_res1D / _J1D;
            }
            else
            {
                cerr<<__func__<<" >> ("<< (applyLS?"with":"without") <<" line search)"
               <<" The J is singular."<<endl;
                return 1;
            }
        }

        _alpha = 1.0;  // withuot Line Search
        if ( applyLS ) // Apply Line Search
        {
            // Check the size of the STEP
            if ( std::isnan(_dlam1D) )
            {
                _dlam1D = 0.0001*_lam1D;
            }
            if ( _dlam1D > step )
            {
                _dlam1D = ( 0.5*step / _dlam1D ) * _dlam1D;
            }

            this->BrentLineSearch1D( _maxIterLS, _tolLS, _alpha );
        }

        // Update lambda
        _lam1D += _alpha * _dlam1D;
    }

    cerr<<__func__<<" >> Maximum number of iterations "<< (applyLS?"with":"without")
        <<" line search reached ("<<_maxIterNR<<")."<<endl;
    return 2;
}


int LocalMaxEntropy::BrentLineSearch1D( int maxIter, double TolX, double &alpha )
{
    //Initialize variables
    double a   = _lowerLimit;
    double b   = _upperLimit;
    double xm  = 0.0;

    double u   = 0.0;
    double v   = a + _gold*(b-a);
    double w   = v;
    double x   = v;

    double e   = 0.0;  // The distance moved on the step before last
    double d  = 0.0;

    double fx  = this->LogZ1D( x );
    double fv  = this->LogZ1D( v );
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
            fu = this->LogZ1D( u ); //The only function evaluation per iteration

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

    throw("Maximum iteration is reached in Brent.");
    return 1;
}


double LocalMaxEntropy::LogZ1D ( double alpha )
{
    // lam_ls = lam + t*dlam
    _lamLS1D = _lam1D + alpha*_dlam1D;

    // exp(-beta*dx^2+lam*dx)
    _Z = 0.0;
    for ( int iN = 0; iN < _knn; ++iN )
        _Z += exp( _beta_dx2.at(iN) + _lamLS1D*_dx.at(iN) );

    // return ln of the partition function
    return log(_Z);
}


// Calculat the gradient matrix of p_a. [N x D]
void LocalMaxEntropy::ComputeGradients1D( )
{
    // dp_s = -p_s*Ji*dx
    for ( int iN = 0; iN < _knn; ++iN )
        _dpa.at(iN) = -_pa.at(iN) * _Jidx.at(iN);
}


//  Calculation of the hessian matrix of p_a. [N x D^2]
void LocalMaxEntropy::ComputeHessians1D( )
{
    int    iN, jN;
    double Dij;

    // Calculate the tensor Ma = pa * j_a x j_a
    for ( iN = 0; iN < _knn; ++iN )
        _M.at(iN) = _pa.at(iN) * _Jidx.at(iN) * _Jidx.at(iN);

    // hp  = M_a - pa * sum(Mb*(1 + D_ij))
    _hpa = _M;
    for ( iN = 0; iN < _knn; ++iN )
        for ( jN = 0; jN < _knn; ++jN )
        {
            Dij = _dx.at(jN) * _Jidx.at(iN);
            _hpa.at(iN) -= _pa.at(iN)*(1.0+Dij) * _M.at(jN);
        }
}





//______________________________________________________________________________
//------------------------------------------------------------------------------
// Local maximum entropy in two-dimensional space
//------------------------------------------------------------------------------

// Compute the LME Shape functions and its first and second spatial derivatives
int LocalMaxEntropy::ComputeLme2D(const vector<double>& xSample , const vector<int> &idxNN,
                                 bool bGradient, bool bHessian )
{
    // if the number of nearest neighbors are not the same as before,
    // resize all the correspondent vectors
    if ( _knn != idxNN.size())
    {
        _knn = idxNN.size();
        if (_knn < _dim)
        {
            cerr<<__func__<<" >>  Maximum number of nearest neighbors is less than the dimension."<<endl;
            return 1;
        }
        _beta_dx2.resize( _knn );
        _dx      .resize( _knn*_dim );
        _Jidx    .resize( _knn*_dim );
        _dxdx    .resize( _knn*_dim2 );
        _M       .resize( _knn*_dim2 );
        _pa      .resize( _knn );
        if (bGradient)
        {
            _dpa.resize( _knn*_dim );
            if (bHessian)
                _hpa.resize( _knn*_dim2 );
        }
    }

    // dx = x - x_a, loop over the neighbors of x_sample
    for ( int iN = 0; iN < _knn; ++iN )
    {
        _dx.at(iN*_dim+0) = xSample.at(0) - _xNodes.at(idxNN.at(iN)*_dim+0);
        _dx.at(iN*_dim+1) = xSample.at(1) - _xNodes.at(idxNN.at(iN)*_dim+1);
    }

    // Shape functions
    _ierr = this->ComputeShapeFunctions2D( );
    if (_ierr)
    {
        cerr<<__func__<<" >> *** LME failed at xi = [";
        for (const auto &iter : xSample )
            cerr<<iter<<" ";
        cerr<<"]"<<endl;
    }

    if ( bGradient && _ierr == 0 )
    {
        // [DxD] J = sum(pa*(dx dx)) - r x r
        _Jchol.at(0) = cblas_ddot( _knn, _pa.data(), 1, &_dxdx.at(0), _dim2 ) - _res.at(0) * _res.at(0);
        _Jchol.at(1) = cblas_ddot( _knn, _pa.data(), 1, &_dxdx.at(1), _dim2 ) - _res.at(0) * _res.at(1);
        //_Jchol.at(2) = _Jchol.at(1); // is not necessary
        _Jchol.at(3) = cblas_ddot( _knn, _pa.data(), 1, &_dxdx.at(3), _dim2 ) - _res.at(1) * _res.at(1);

        _detJ = _Jchol.at(0)*_Jchol.at(3) - _Jchol.at(1)*_Jchol.at(1) ;
//        cout<<"detJ2 = "<<_detJ<<endl;


        double *pdx   = nullptr;
        // [NxD] Ji_dx = invJ * (x-x_a)
        for ( int iN = 0; iN < _knn; ++iN )
        {
            pdx   = &_dx.at(iN*_dim);
            // Jidx = Jinv * dx
            _Jidx.at(iN*_dim  ) = ( _Jchol.at(3)*pdx[0] - _Jchol.at(1)*pdx[1] ) / _detJ;
            _Jidx.at(iN*_dim+1) = ( _Jchol.at(0)*pdx[1] - _Jchol.at(1)*pdx[0] ) / _detJ;
        }

        // Gradient of the shape functions
        this->ComputeGradients2D( );

        // Hessian of the shape functions
        if ( bHessian ) this->ComputeHessians2D( );
    }
    return _ierr;
}


// Calculate the LME shape functions in two dimensional space
int LocalMaxEntropy::ComputeShapeFunctions2D( )
{
    int    ierr   = 0;
    double step   = 0.0;
    double *pdx   = nullptr;
    double *pdxdx = nullptr;

    for ( int iN = 0; iN < _knn; ++iN )
    {
        pdx   = &_dx.at(iN*_dim);
        pdxdx = &_dxdx.at(iN*_dim2);

        // [N*D2] dxdx = dx*dx (tensor produxt)
        // pdxdx[2] = pdxdx[1] is not necessary
        pdxdx[0] = pdx[0] * pdx[0];
        pdxdx[1] = pdx[0] * pdx[1];
        pdxdx[3] = pdx[1] * pdx[1];

        // step = sum(-beta*dx2) and then normalize with norm(dx)
        step += _beta_dx2.at(iN) = -_beta * ( pdxdx[0] + pdxdx[3] );
    }


    ierr = this->LmeNewton2D( false, 0.0 ); // withouth line search

    if ( ierr )
    {
        // step = norm(-beta*dx2)/norm(dx);
        step = fabs(step) / cblas_dnrm2(_knn*_dim, _dx.data(), 1);
        ierr = this->LmeNewton2D( true, step ); // with line search
    }

    return ierr;

}


// Local max-entropy shape functions are computed by using Newton-Raphson method
// with possible Brent's line search algorithm
int LocalMaxEntropy::LmeNewton2D( bool applyLS, double step )
{
    int iN;

    _lam.at(0) = 0.0;
    _lam.at(1) = 0.0;

    for ( int iter = 0; iter < _maxIterNR; ++iter )
    {
        _Z = 0.0;
        // exp(-beta*dx^2+lam.dx)
        for ( iN = 0; iN < _knn; ++iN )
            _Z += _pa.at(iN) = exp( _beta_dx2.at(iN) + _lam.at(0)*_dx.at(iN*_dim  )
                                                     + _lam.at(1)*_dx.at(iN*_dim+1) );

        // Compute shape functions (normalize with Z)
        cblas_dscal( _pa.size(), 1.0/_Z, _pa.data(), 1 );

        // residual vector: r = sum(pa*dx)
        _res.at(0) = cblas_ddot( _knn, _pa.data(), 1, &_dx.at(0), _dim);
        _res.at(1) = cblas_ddot( _knn, _pa.data(), 1, &_dx.at(1), _dim);

        _norm_res = sqrt(_res.at(0)*_res.at(0)+_res.at(1)*_res.at(1));

        if ( std::isnan( _norm_res ) )
        {
            cerr<<__func__<<" >> ("<< (applyLS?"with":"without") <<" line search)"
               <<"  Norm of res vector is NaN."<<endl;
            cout<<__func__<<" >> ("<< (applyLS?"with":"without") <<" line search)"
               <<"  Norm of res vector is NaN."<<endl;
            return 3;
        }
        else if ( _norm_res < _tolNR )
        {
            // The solution is converged.
            return 0;
        }
        else // Construct J and solve the system of J*dlam = -res
        {
            // [DxD] J = sum(pa*(dx dx)) - r x r
            //_Jchol.at(2) = _Jchol.at(1) is not necessary
            _Jchol.at(0) = cblas_ddot( _knn, _pa.data(), 1, &_dxdx.at(0), _dim2 ) - _res.at(0) * _res.at(0);
            _Jchol.at(1) = cblas_ddot( _knn, _pa.data(), 1, &_dxdx.at(1), _dim2 ) - _res.at(0) * _res.at(1);
            _Jchol.at(3) = cblas_ddot( _knn, _pa.data(), 1, &_dxdx.at(3), _dim2 ) - _res.at(1) * _res.at(1);

            _detJ = _Jchol.at(0)*_Jchol.at(3) - _Jchol.at(1)*_Jchol.at(1) ;
//            cout<<"detJ = "<<_detJ<<endl;
            if ( _detJ > _epsilon )
            {
                // ***update dlam***
                // dlam = -J^{-1}*res
                _dlam.at(0) = ( _Jchol.at(1)*_res.at(1) - _Jchol.at(3)*_res.at(0) ) / _detJ;
                _dlam.at(1) = ( _Jchol.at(1)*_res.at(0) - _Jchol.at(0)*_res.at(1) ) / _detJ;
            }
            else
            {
                cerr<<__func__<<" >> ("<< (applyLS?"with":"without") <<" line search)"
                   <<" The J is singular."<<endl;
                cout<<__func__<<" >> ("<< (applyLS?"with":"without") <<" line search)"
                   <<" The J is singular."<<endl;
                return 1;
            }

        }

        _alpha = 1.0;  // withuot Line Search
        if ( applyLS ) // Apply Line Search
        {
            // Check the size of the STEP
            _norm_dlam = sqrt( _dlam.at(0)*_dlam.at(0) + _dlam.at(1)*_dlam.at(1) );
            if ( std::isnan(_norm_dlam) )
            {
                // dlam = 0.001*lambda
                _dlam.at(0) = 0.0001*_lam.at(0);
                _dlam.at(1) = 0.0001*_lam.at(1);
            }
            if ( _norm_dlam > step )
            {
                // dlam = (0.5*step/norm(dlam)) * dlam
                _dlam.at(0)  = _dlam.at(1) =  0.5 * step /_norm_dlam ;
                _dlam.at(0) *= _lam.at(0);
                _dlam.at(1) *= _lam.at(1);
            }

            this->BrentLineSearch2D( _maxIterLS, _tolLS, _alpha );
        }

        // Update lambda = lambda + alpha*dlam
        _lam.at(0) += _alpha * _dlam.at(0);
        _lam.at(1) += _alpha * _dlam.at(1);
    }

    cerr<<__func__<<" >> Maximum number of iterations "<< (applyLS?"with":"without")
        <<" line search reached ("<<_maxIterNR<<")."<<endl;
    cout<<__func__<<" >> Maximum number of iterations "<< (applyLS?"with":"without")
        <<" line search reached ("<<_maxIterNR<<")."<<endl;
    return 2;
}


int LocalMaxEntropy::BrentLineSearch2D( int maxIter, double TolX, double &alpha )
{
    //Initialize variables
    double a   = _lowerLimit;
    double b   = _upperLimit;
    double xm  = 0.0;

    double u   = 0.0;
    double v   = a + _gold*(b-a);
    double w   = v;
    double x   = v;

    double e   = 0.0;  // The distance moved on the step before last
    double d  = 0.0;

    double fx  = this->LogZ2D( x );
    double fv  = this->LogZ2D( v );
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
            fu = this->LogZ2D( u ); //The only function evaluation per iteration

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

    cerr<<__func__<<" >> Maximum number of iterations reached ("<<maxIter<<")."<<endl;
    cout<<__func__<<" >> Maximum number of iterations reached ("<<maxIter<<")."<<endl;
    return 1;
}


double LocalMaxEntropy::LogZ2D( double alpha )
{
    // lam_ls = lam + alpha* dlam
    _lamLS.at(0) = _lam.at(0) + alpha * _dlam.at(0);
    _lamLS.at(1) = _lam.at(1) + alpha * _dlam.at(1);

    //sum( exp(-beta*dx^2+lam*dx) )
    _Z = 0.0;
    for ( int iN = 0; iN < _knn; ++iN )
        _Z += exp( _beta_dx2.at(iN)
                   + _lamLS.at(0)*_dx.at(iN*_dim  )
                   + _lamLS.at(1)*_dx.at(iN*_dim+1) );

    // return ln of the partition function
    return log(_Z);
}


// Calculat the gradient matrix of p_a. [N x D]
void LocalMaxEntropy::ComputeGradients2D( )
{
    // dp_s = -p_s*Ji*dx
    for ( int iN = 0; iN < _knn; ++iN )
    {
        _dpa.at(iN*_dim  ) = -_pa.at(iN)*_Jidx.at(iN*_dim  );
        _dpa.at(iN*_dim+1) = -_pa.at(iN)*_Jidx.at(iN*_dim+1);
    }
}


//  Calculation of the hessian matrix of p_a. [N x D^2]
void LocalMaxEntropy::ComputeHessians2D( )
{
    int    iN, jN;
    double Dij;
    double *pJidx;
    // Calculate the tensor Ma = pa * j_a x j_a
    for ( iN = 0; iN < _knn; ++iN )
    {
        pJidx = &_Jidx.at(iN*_dim);
        _M.at(iN*_dim2+0) = _pa.at(iN) * pJidx[0] * pJidx[0];
        _M.at(iN*_dim2+1) = _pa.at(iN) * pJidx[0] * pJidx[1];
        _M.at(iN*_dim2+3) = _pa.at(iN) * pJidx[1] * pJidx[1];
    }

    // hp  = M_a - pa * sum(Mb*(1 + D_ij))
    _hpa = _M;
    for ( iN = 0; iN < _knn; ++iN )
    {
        for ( jN = 0; jN < _knn; ++jN )
        {
            // D_ij = dx(j).Jidx(i)
            Dij = _dx.at(jN*_dim  ) * _Jidx.at(iN*_dim  )
                + _dx.at(jN*_dim+1) * _Jidx.at(iN*_dim+1);
            _hpa.at(iN*_dim2+0) -= _pa.at(iN)*(1.0+Dij) * _M.at(jN*_dim2+0);
            _hpa.at(iN*_dim2+1) -= _pa.at(iN)*(1.0+Dij) * _M.at(jN*_dim2+1);
            _hpa.at(iN*_dim2+3) -= _pa.at(iN)*(1.0+Dij) * _M.at(jN*_dim2+3);
        }
        _hpa.at(iN*_dim2+2) = _hpa.at(iN*_dim2+1);
    }
}







//______________________________________________________________________________
//------------------------------------------------------------------------------
// Local maximum entropy in multi-dimensional space
//------------------------------------------------------------------------------

// Compute the LME Shape functions and its first and second spatial derivatives
int LocalMaxEntropy::ComputeLme(const vector<double>& xSample , const vector<int> &idxNN,
                                 bool bGradient, bool bHessian )
{
    int ierr;
    // if the number of nearest neighbors are not the same as before,
    // resize all the correspondent vectors
    if ( _knn != idxNN.size())
    {
        _knn      = idxNN.size();
        if (_knn < _dim)
        {
            cerr<<__func__<<" >>  Maximum number of nearest neighbors is less than the dimension."<<endl;
            return 1;
        }
        _lam_dx  .resize( _knn );
        _beta_dx2.resize( _knn );
        _dx      .resize( _knn*_dim );
        _Jidx    .resize( _knn*_dim );
        _dxdx    .resize( _knn*_dim2 );
        _M       .resize( _knn*_dim2 );
        _pa      .resize( _knn );
        if (bGradient)
        {
            _dpa.resize( _knn*_dim, 0.0);
            if (bHessian)
                _hpa.resize( _knn*_dim2, 0.0 );
        }
    }

    // dx = x - x_a, loop over the neighbors of x_sample
    for ( int iN = 0; iN < _knn; ++iN )
        transform( xSample.begin(), xSample.end(),
                   _xNodes.begin()+idxNN.at(iN)*_dim,
                   _dx.begin()+iN*_dim, std::minus<double>() );

    // Shape functions
    _ierr = this->ComputeShapeFunctions( );

    if ( bGradient && _ierr == 0 )
    {
        // Compute [NxD] Jidxs = invJ * (x-x_a)

        // J[DxD] = sum(pa*(dx dx)) - R
        // ___sum(pa*(dx dx)) = dxdx' pa
        cblas_dgemv(CblasRowMajor, CblasTrans, _knn, _dim2, 1.0, _dxdx.data(), _dim2,
                    _pa.data(), 1, 0.0, _Jchol.data(), 1);
        // ___-R
        cblas_daxpy( _dim2, -1.0, _R.data(), 1, _Jchol.data(), 1 );

        // Compote Cholesky factorization of J = L*L' ( J is overwritten by L )
        ierr = LAPACKE_dpotrf( LAPACK_COL_MAJOR, 'L', _dim, _Jchol.data(), _dim );
        if (ierr == 0)
        {
            _Jidx = _dx; // initialize to dx
            // Solve J * Jidx = dx
            ierr = LAPACKE_dpotrs( LAPACK_COL_MAJOR, 'L', _dim, _knn, _Jchol.data(), _dim, _Jidx.data(), _dim );
            //TODO: Maybe a refinement with LAPACKE_dporfs
            if (ierr == 0)
            {
                // Gradient of the shape functions
                this->ComputeGradients( );

                // Hessian of the shape functions
                if ( bHessian ) this->ComputeHessians( );
            }
            else if (ierr < 0)
                cerr<<__func__<<" >> The "<<ierr<<"-th parameter of LAPACKE_dpotrs had an illegal value."<<endl;
            else if (ierr > 0)
                cerr<<__func__<<" >> Unexpexted error in LAPACKE_dpotrs ( info ="<<ierr<<")."<<endl;
        }
        else if (ierr < 0)
            cerr<<__func__<<" >> The "<<ierr<<"-th parameter of LAPACKE_dpotrf had an illegal value."<<endl;
        else if (ierr > 0)
            cerr<<__func__<<" >> The J matrix (the leading minor of order "<<ierr<<") is not positive definite."<<endl;
        _ierr = ierr ? 2 : 0;
    }
    return _ierr;
}


// Calculate the local maximum-entropy shape functions for a given sample point,
// by employing Newton method with possible Brent's line search
int LocalMaxEntropy::ComputeShapeFunctions( )
{
    int     iD, jD;
    int     ierr = 0;
    double  step = 0.0;
    double *pdx  = nullptr;

    for ( int iN = 0; iN < _knn; ++iN )
    {
        pdx = &_dx.at(iN*_dim);

        // [N*D2] dxdx = dx*dx (tensor produxt)
        for ( iD = 0; iD < _dim; ++iD )
            for ( jD = 0; jD < _dim; ++jD )
                _dxdx.at(iN*_dim2 + iD*_dim + jD) = _dx.at(iN*_dim + iD) * _dx.at(iN*_dim + jD);
//        cblas_dgemm ( CblasRowMajor, CblasTrans, CblasNoTrans, _dim, _dim, 1,
//                      1.0, pdx, _dim, pdx, _dim, 0.0 , &_dxdx.at(iN*_dim2), _dim );

        // step = sum(-beta*dx2) and then normalize with norm(dx)
        step += _beta_dx2.at(iN) = -_beta * cblas_ddot( _dim, pdx, 1, pdx, 1);
    }


    ierr = this->LmeNewton( false, 0.0 ); // withouth line search

    if ( ierr )
    {
        // step = norm(-beta*dx2)/norm(dx);
        step = fabs(step) / cblas_dnrm2(_knn*_dim, _dx.data(), 1);
        ierr = this->LmeNewton( true, step ); // with line search
    }

    return ierr;
}


// Local max-entropy shape functions are computed by using Newton-Raphson method
// with possible Brent's line search algorithm
int LocalMaxEntropy::LmeNewton ( bool applyLS, double step )
{
    int    iN, iD, jD;
    int    ierr      = 0;

    fill(_lam.begin(),_lam.end(), 0.0);

    for ( int iter = 0; iter < _maxIterNR; ++iter )
    {
        // exp(-beta*dx^2+lam.dx)
        _Z = 0.0;
        for ( iN = 0; iN < _knn; ++iN )
            _Z += _pa.at(iN) = exp( _beta_dx2.at(iN)
                                    + cblas_ddot( _dim, _lam.data(), 1, &_dx.at(iN*_dim), 1 ) );

        // Compute shape functions (normalize with Z)
        cblas_dscal( _knn, 1.0/_Z, _pa.data(), 1 );

        // residual vector: res = sum(pa*dx) (or res = dx'pa) [Dx1]
        cblas_dgemv(CblasRowMajor, CblasTrans, _knn, _dim, 1.0, _dx.data(), _dim,
                    _pa.data(), 1, 0.0, _res.data(), 1);

        _norm_res = cblas_dnrm2(_dim, _res.data(), 1 );

        // R[DxD] = res' x res
        for ( iD = 0; iD < _dim; ++iD )
            for ( jD = 0; jD < _dim; ++jD )
                _R.at(iD*_dim+jD) = _res.at(iD) * _res.at(jD);
//        cblas_dgemm ( CblasRowMajor, CblasTrans, CblasNoTrans, _dim, _dim, 1,
//                      1.0, _res.data(), _dim, _res.data(), _dim , 0.0, _R.data(), _dim );

        if ( std::isnan( _norm_res ) )
        {
            cerr<<__func__<<" >> "<< (applyLS?"with":"without") <<" line search)"
               <<" Norm of res vector is NaN."<<endl;
            return 3;
        }
        else if ( _norm_res < _tolNR )
        {
            // The solution is converged.
            return 0;
        }
        else // Construct J and solve the system of J*dlam = res
        {
            // J[DxD] = sum(pa*(dx dx)) - R
            // ___sum(pa*(dx dx)) = dxdx' pa
            cblas_dgemv(CblasRowMajor, CblasTrans, _knn, _dim2, 1.0, _dxdx.data(), _dim2,
                        _pa.data(), 1, 0.0, _Jchol.data(), 1);
            // ___-R
            cblas_daxpy( _dim2, -1.0, _R.data(), 1, _Jchol.data(), 1 );

            // Compote Cholesky factorization of J = L*L' ( J is overwritten by L )
            ierr = LAPACKE_dpotrf( LAPACK_COL_MAJOR, 'L', _dim, _Jchol.data(), _dim );
            if (ierr == 0)
            {
                // Solve L*L' dlam = res
                ierr = LAPACKE_dpotrs( LAPACK_COL_MAJOR, 'L', _dim, 1, _Jchol.data(), _dim, _res.data(), _dim );
                //TODO: Maybe a refinement with LAPACKE_dporfs
                if (ierr == 0)
                {
                    // ***update dlam***
                    cblas_daxpby(_dim, -1.0, _res.data(), 1, 0.0, _dlam.data(), 1 );
                }
                else if (ierr < 0)
                {
                    cerr<<__func__<<" >> The "<<ierr<<"-th parameter of LAPACKE_dpotrs had an illegal value."<<endl;
                    return 1;
                }
                else if (ierr > 0)
                {
                    cerr<<__func__<<" >> Unexpexted error in LAPACKE_dpotrs ( info ="<<ierr<<")."<<endl;
                    return 1;
                }
            }
            else if (ierr < 0)
            {
                cerr<<__func__<<" >> The "<<ierr<<"-th parameter of LAPACKE_dpotrf had an illegal value."<<endl;
                return 1;
            }
            else if (ierr > 0)
            {
                cerr<<__func__<<" >> The J matrix (the leading minor of order "<<ierr<<") is not positive definite."<<endl;
                return 1;
            }
        }

        _alpha = 1.0;  // withuot Line Search
        if ( applyLS ) // Apply Line Search
        {
            // Check the size of the STEP
            _norm_dlam = cblas_dnrm2(_dim, _dlam.data(), 1);
            if ( std::isnan(_norm_dlam) )
            {
                // dlam = 0.001*lambda
                _dlam = _lam;
                cblas_dscal( _dim, 0.001, _dlam.data(), 1);
            }
            if ( _norm_dlam > step )
            {
                // dlam = (0.5*step/val) * dlam
                cblas_dscal( _dim, 0.5*step/_norm_dlam, _dlam.data(), 1);
            }

            this->BrentLineSearch( _maxIterLS, _tolLS, _alpha );
        }

        // Update lambda = lambda + alpha*dlam
        cblas_daxpy(_dim, _alpha, _dlam.data(), 1, _lam.data(), 1);
    }

    cerr<<__func__<<" >> Maximum number of iterations "<< (applyLS?"with":"without")
        <<" line search reached ("<<_maxIterNR<<")."<<endl;
    return 2;
}


int LocalMaxEntropy::BrentLineSearch( int maxIter, double TolX, double &alpha )
{
    //Initialize variables
    double a   = _lowerLimit;
    double b   = _upperLimit;
    double xm  = 0.0;

    double u   = 0.0;
    double v   = a + _gold*(b-a);
    double w   = v;
    double x   = v;

    double e   = 0.0;  // The distance moved on the step before last
    double d  = 0.0;

    double fx  = this->LogZgeneral( x );
    double fv  = this->LogZgeneral( v );
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
            fu = this->LogZgeneral( u ); //The only function evaluation per iteration

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

    throw("Maximum iteration is reached in Brent.");
    return 1;
}


double LocalMaxEntropy::LogZgeneral ( double alpha )
{
    // lam_ls = lam + alpha*dlam
    _lamLS = _lam;
    cblas_daxpy( _dim, alpha, _dlam.data(), 1, _lamLS.data(), 1 );

    // lamls*dx
    cblas_dgemv(CblasRowMajor, CblasNoTrans, _knn, _dim, 1.0,
                _dx.data(), _dim, _lamLS.data(), 1, 0.0, _lam_dx.data(), 1 );

    // exp(-beta*dx^2+lam*dx)
    _Z = 0.0;
    for ( int iN = 0; iN < _knn; ++iN )
        _Z += exp( _beta_dx2.at(iN) + _lam_dx.at(iN) );

    // return ln of the partition function
    return log(_Z);
}


// Calculat the gradient matrix of p_a. [N x D]
void LocalMaxEntropy::ComputeGradients( )
{
    // dp_a = -p_a*Ji*dx
    _dpa = _Jidx;
    for ( int iN = 0; iN < _knn; ++iN )
        cblas_dscal( _dim, -_pa.at(iN), &_dpa.at(iN*_dim), 1);
}


//  Calculation of the hessian matrix of p_a. [N x D^2]
void LocalMaxEntropy::ComputeHessians( )
{
    int     iD, jD;
    int     iN, jN;
    double  Dij;
    double *pJidx = nullptr;

    // Calculate the tensor M_i = p_i * j_a x j_a
    for ( iN = 0; iN < _knn; ++iN )
    {
        pJidx = &_Jidx.at(iN*_dim);
        // M_i = j_a x j_a
        for ( iD = 0; iD < _dim; ++iD )
            for ( jD = 0; jD < _dim; ++jD )
                _M.at(iN*_dim2 + iD*_dim + jD) = _Jidx.at(iN*_dim + iD) * _Jidx.at(iN*_dim + jD);
//        cblas_dgemm ( CblasRowMajor, CblasTrans, CblasNoTrans, _dim, _dim, 1,
//                      1.0, &_Jidx.at(aN*_dim), _dim, &_Jidx.at(aN*_dim), _dim,
//                      0.0, &_M.at(aN*_dim2), _dim );
        // M_i = p_i * j_a x j_a
        cblas_dscal(_dim2, _pa.at(iN), &_M.at(iN*_dim2), 1 );
    }

    // initialize hp to M
    _hpa = _M;

    // Hessian value at each node
    for ( iN = 0; iN < _knn; ++iN )
    {
        // hp  = M_i - p_i * sum(M_j*(1 + D_ij))
        for ( jN = 0; jN < _knn; ++jN )
        {
            Dij = cblas_ddot( _dim, &_dx.at(jN*_dim), 1, &_Jidx.at(iN*_dim), 1 );
            cblas_daxpy( _dim2, -_pa.at(iN)*(1.0+Dij), &_M.at(jN*_dim2), 1, &_hpa.at(iN*_dim2), 1 );
        }
    }
}

