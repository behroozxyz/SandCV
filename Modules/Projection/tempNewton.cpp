// HEADER
    /* Newton method variable */
    double _alpha;              ///< step size in Newton method
    double _tolNR  = 1.0e-14;   ///< Tolerence of Newton-Raphson method
    double _sqrtEpsilon;        ///< square root of epsilon
    double _gold;               ///< Golden Ratio: 0.3819660
    static const double _epsilon = numeric_limits<double>::epsilon(); ///< Machine precision
    double         _norm_res;    ///< norm of gradient of cost function\f$ \vert\nabla f(\xi)\vert\f$
    double         _norm_dxi;    ///< norm of delta xi \f$ \vert\Delta\xi\vert\f$
    vector<double> _xiLS;        ///< xi value to update the line search
    vector<double> _fGrad;      ///< gradient of cost function */
    vector<double> _fHess;      ///< hessian of cost function */

    double CostValue(double alpha);

    /** Returns the absolute value of A times the sign of B */
    inline double sign(double absA, double signB) const { return (signB < 0) ? -fabs(absA) : fabs(absA); }

    /**
     * \param applyLS    line search flag
     * \param step       At first a checking for the step size of the JUMP is made
     * \return error identifier:
     *                       0 - successful
     *                       1 - error in inverting J
     *                       2 - maximum number of iteration reached
     *                       3 - there is NaN in the res
     */
    int Newton(bool applyLS);


    /** Perform line search with Brent's minimization algorithm (based on \cite Press2007)*/
    int BrentLineSearch(int maxIter, double TolX, double &alpha);






// SOURCE


double ClosestPointProjection::CostValue( double alpha )
{
    // xi_ls = xi + alpha*dxi
    _xiLS = _xiVal;
    cblas_daxpy( _ldim, alpha, _fGrad.data(), 1, _xiLS.data(), 1 );

    if (_parametrization->ComputeLme( _xiLS, 0 ))
    {
        cerr<<__func__<<" >> error in calculating shape functions."<<endl;
        return 1;
    }
    else
    {
        _pa    = _parametrization->GetShapeFunctions();
        _idxNN = _parametrization->GetNearestNeighbors();
        _knn   = _idxNN.size();
        fill(_xClosest.begin(), _xClosest.end(), 0.0);
        for ( int iN = 0; iN < _knn; ++iN )
            cblas_daxpy( _hdim, _pa.at(iN), &_xNodes.at(_idxNN.at(iN)*_hdim), 1,
                         _xClosest.data(), 1);
        transform( _xClosest.begin() , _xClosest.end(), _xSample.begin(),
                   _xDiff.begin(), std::minus<double>() );
        return cblas_dnrm2( _hdim, _xDiff.data(), 1);
    }
}


int ClosestPointProjection::Newton( bool applyLS )
{
    int    iN, iD;
    int    ierr      = 0;
    int    maxIter   = applyLS ? 1000 : 100; // Maximum number of Newton's iterations
    int    maxIterLS = 100; // Maximum number of iterations for linesearch

    for ( int iter = 0; iter < maxIter; ++iter )
    {
        if (_parametrization->ComputeLme( _xiVal, 2 ))
        {
            cerr<<__func__<<" >> error in calculating shape functions with derivatives."<<endl;
            return 1;
        }
        else
        {
            double *xNode;
            _pa    = _parametrization->GetShapeFunctions();
            _dpa   = _parametrization->GetShapeFunctionsGradient();
            _hpa   = _parametrization->GetShapeFunctionsHessian();
            _idxNN = _parametrization->GetNearestNeighbors();
            _knn   = _idxNN.size();
            fill(_xClosest .begin(), _xClosest .end(), 0.0);
            fill(_dxClosest.begin(), _dxClosest.end(), 0.0);
            fill(_hxClosest.begin(), _hxClosest.end(), 0.0);
            for ( iN = 0; iN < _knn; ++iN )
            {
                xNode = &_xNodes.at(_idxNN.at(iN)*_hdim);
                cblas_daxpy( _hdim, _pa.at(iN), xNode, 1, _xClosest.data(), 1);
                for ( iD = 0; iD < _ldim ; ++iD )
                    cblas_daxpy( _hdim, _dpa.at(iN*_ldim +iD), xNode, 1, &_dxClosest.at(iD*_hdim), 1);
                for ( iD = 0; iD < _ldim2; ++iD )
                    cblas_daxpy( _hdim, _hpa.at(iN*_ldim2+iD), xNode, 1, &_hxClosest.at(iD*_hdim), 1);
            }
        }
        /* diff = x - xSample*/
        transform( _xClosest.begin() , _xClosest.end(), _xSample.begin(),
                   _xDiff.begin(), std::minus<double>() );

        /* -df = -dx diff */
        cblas_dgemv( CblasRowMajor, CblasNoTrans, _ldim, _hdim, -1.0,
                     _dxClosest.data(),_hdim, _xDiff.data(), 1, 0.0, _fGrad.data(), 1);

        /* hf = hx diff + dx dx  */
        cblas_dgemv( CblasRowMajor, CblasNoTrans, _ldim2, _hdim, 1.0,
                     _hxClosest.data(),_hdim, _xDiff.data(), 1, 0.0, _fHess.data(), 1);
        cblas_dgemm( CblasRowMajor, CblasNoTrans, CblasTrans, _ldim, _ldim, _hdim,
                     1.0, _dxClosest.data(), _hdim, _dxClosest.data(), _hdim ,
                     1.0, _fHess.data(), _ldim );

        /* residual = nomr(df) */
        _norm_res = cblas_dnrm2( _ldim, _fGrad.data(), 1 );

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
        else /* Solve the system of hf*dxi = -df */
        {
            /* Compote Cholesky factorization of fHess = L*L' ( fHess is overwritten by L ) */
            ierr = LAPACKE_dpotrf( LAPACK_ROW_MAJOR, 'L', _ldim, _fHess.data(), _ldim );
            if (ierr == 0)
            {
                /* Solve L*L' dxi = -df ( fGrad is overwritten by dxi ) */
                ierr = LAPACKE_dpotrs( LAPACK_ROW_MAJOR, 'L', _ldim, 1, _fHess.data(), _ldim, _fGrad.data(), 1 );
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
            /* Check the size of the STEP */
            _norm_dxi = cblas_dnrm2(_ldim, _fGrad.data(), 1);
            if ( std::isnan(_norm_dxi) )
            {
                /* dlam = 0.001*lambda */
                _fGrad = _xiVal;
                cblas_dscal( _ldim, 0.001, _fGrad.data(), 1);
            }
//            double step = 2.0;
//            if ( _norm_dxi > step )
//            {
//                /* dlam = (0.5*step/val) * dlam */
//                cblas_dscal( _ldim, 0.5*step/_norm_dxi, _fGrad.data(), 1);
//            }
            this->BrentLineSearch( maxIterLS, _tolNR*1e3, _alpha );
        }

        /* Update lambda = lambda + alpha*dlam */
        cblas_daxpy(_ldim, _alpha, _fGrad.data(), 1, _xiVal.data(), 1);
    }

    cerr<<__func__<<" >> Maximum number of iterations "<< (applyLS?"with":"without")
        <<" line search reached ("<<maxIter<<")."<<endl;
    return 2;
}


int ClosestPointProjection::BrentLineSearch( int maxIter, double TolX, double &alpha )
{
    //Initialize variables
    double a   = -0.01; // lower limit
    double b   =  5.00; // upper limit
    double xm  =  0.0;

    double u   = 0.0;
    double v   = a + _gold*(b-a);
    double w   = v;
    double x   = v;

    double e   = 0.0;  // The distance moved on the step before last
    double d  = 0.0;

    double fx  = this->CostValue( x );
    double fv  = this->CostValue( v );
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
            fu = this->CostValue( u ); //The only function evaluation per iteration

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
    return 1;
}


