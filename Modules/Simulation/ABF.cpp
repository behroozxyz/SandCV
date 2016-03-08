/*=========================================================================
 * File         : ABF.cpp
 * Module       : Simulation
 * Copyright    : (C)opyright 2013
 *                See COPYRIGHT statement in top level directory.
 * Authors      : Behrooz Hashemian
 * Created      : Tue Aug 13, 2013  11:49AM
 * Last modified: Tue Aug 13, 2013  11:49AM
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
#include"ABF.h"
static const double _epsilon = numeric_limits<double>::epsilon(); ///< Machine precision

ABF::ABF(bool bApplyABF, bool bLoadHistogram, int ldim, int sdim, int nAtoms, double timeStep,
         int logFrequency, int histFrequency,
         double corralHeight, string corralFilename,
         string coreName,
         const vector<double> &xiMin, const vector<double> &xiMax,
         const vector<int> &nBins, int fullSample, int rampSample)
    :
      _bApplyABF(bApplyABF),
      _bLoadHistogram(bLoadHistogram),
      _ldim(ldim),
      _sdim(sdim),
      _nAtoms(nAtoms),
      _logFrequency(logFrequency),
      _histFrequency(histFrequency),
      _timeStepInv(1.0/timeStep),
      _corral(corralFilename, corralHeight),
      _xiMin(xiMin),
      _xiMax(xiMax),
      _fullSample(fullSample),
      _rampSample(rampSample),
      _nBins(nBins),
      _coreName(coreName)
{
    _bCorral = corralHeight > 0.0 ? true : false;
    this->Update();
}



ABF::ABF(bool bApplyABF, bool bLoadHistogram, int ldim, int sdim, int nAtoms, double timeStep,
          int logFrequency, int histFrequency,
          double kWall, const vector<double> &lowerWall, const vector<double> &upperWall,
          string coreName,
         const vector<double> &xiMin, const vector<double> &xiMax, const vector<int> &nBins, int fullSample, int rampSample)
    :
      _bApplyABF(bApplyABF),
      _bLoadHistogram(bLoadHistogram),
      _ldim(ldim),
      _sdim(sdim),
      _nAtoms(nAtoms),
      _logFrequency(logFrequency),
      _histFrequency(histFrequency),
      _timeStepInv(1.0/timeStep),
      _kWall(kWall),
      _lowerWall(lowerWall),
      _upperWall(upperWall),
      _xiMin(xiMin),
      _xiMax(xiMax),
      _fullSample(fullSample),
      _rampSample(rampSample),
      _nBins(nBins),
      _coreName(coreName)
{
    _bHarmonic = _kWall > 0.0 ? true : false;
    this->Update();
}



ABF::~ABF()
{
    WriteHistograms("_final");
    _logFile.close();
}


void ABF::Update( )
throw( invalid_argument )
{
    if ( this->CheckSettings() )
    {
        cerr<<__func__<<" >> CheckSettings went wrong."<<endl;
        throw invalid_argument("One or more inputs are wrong.");
    }
    _bBoundaryPotential = (_bCorral||_bHarmonic) ? true : false;
    _hdim = _sdim*_nAtoms;

    _velocity          .resize(_hdim,       0.0);
    _positionOld       .resize(_hdim,       0.0);
    _force             .resize(_hdim,       0.0);
    _xiRange           .resize(_ldim,       0.0);
    _xiBin             .resize(_ldim,       0  );
    _xiVelocity        .resize(_ldim,       0.0);
    _xiMomentum        .resize(_ldim,       0.0);
    _xiMomentumOld     .resize(_ldim,       0.0);
    _xiForce           .resize(_ldim,       0.0);
    _xiTotalForce      .resize(_ldim,       0.0);
    _xiMass            .resize(_ldim*_ldim, 0.0);
    _xiMassInv         .resize(_ldim*_ldim, 0.0);
    _accNumBins        .resize(_ldim,       0  );
    _binLambda         .resize(_ldim,       0.0);
    _binLambdaOld      .resize(_ldim,       0.0);
    _xiBoundaryForce   .resize(_ldim,       0.0);
    _xiBoundaryForceOld.resize(_ldim,       0.0);

    _logName   = "log_" +_coreName+".txt";
    _histName  = "hist_"+_coreName;

    /* Calculate xi range */
    transform( _xiMax.begin(), _xiMax.end(), _xiMin.begin(), _xiRange.begin(), std::minus<double>() );

    /* Build the histogram */
    this->BuildHistograms();

    /* Load histogram data from file */
    if (_bLoadHistogram) this->LoadHistograms();

    /* Create the log file */
    _logFile.open( _logName );
//    _logFile.close();
}


bool ABF::CheckSettings( )
{
    bool bError = false;

    if (_ldim < 1)
    {
        cerr<<__func__<<" >>  CV's dimension is less than 1("<<_ldim<<")."<<endl;
        bError = true;
    }
    if (_sdim < 1)
    {
        cerr<<__func__<<" >>  Molecule ambient dimension is less than 1 ("<<_sdim<<")."<<endl;
        bError = true;
    }
    if (_nAtoms < 1)
    {
        cerr<<__func__<<" >>  Number of atoms in the molecule is less than 1 ("<<_nAtoms<<")."<<endl;
        bError = true;
    }
    if (_logFrequency < 1)
    {
        cerr<<__func__<<" >>  Output frequency is less than 1 ("<<_logFrequency<<")."<<endl;
        bError = true;
    }
    if (_bCorral&&_bHarmonic)
    {
        cerr<<__func__<<" >>  It is not possible to set both corral and semi-harmonic potential."<<endl;
        bError = true;
    }
    if (_xiMax.size() != _ldim)
    {
        cerr<<__func__<<" >> xiMax dimension ("<<_xiMax.size()<<") is not the same as ldim ("<<_ldim<<")."<<endl;
        bError = true;
    }
    if (_xiMin.size() != _ldim)
    {
        cerr<<__func__<<" >> xiMin dimension ("<<_xiMin.size()<<") is not the same as ldim ("<<_ldim<<")."<<endl;
        bError = true;
    }

    return bError;
}



vector<double> &ABF::CalculateForce(const int &stepNum,
                                    const vector<double> &mass,
                                    const vector<double> &position,
                                    const vector<double> &xiVal,
                                    const vector<double> &xiJac,
                                    bool bMaster,
                                    const vector<double> &transferedLambda)
throw(runtime_error)
{
    _bMaster = bMaster;
    if(_stepNum+1 != stepNum)
        _nSteps = 0;
    _stepNum = stepNum;
    _mass    = mass;
    _xiVal   = xiVal;
    _xiJac   = xiJac;

    if(++_nSteps <= 2)
    {
        if (_nSteps == 1)
        {
            /* Update old positions */
            _positionOld = position;

            /* No force to be appied */
            fill(_force.begin(),_force.end(), 0.0);

            /* Write log */
            if ( _stepNum % _logFrequency == 0)
                this->WriteLogs( );

            return _force;
        }
        else if (_nSteps == 2)
        {
            /* Compute atom velocities v = (r_new-r_old) / dt */
            transform( position.begin(), position.end(), _positionOld.begin(), _velocity.begin(), std::minus<double>() );
            for (auto &iter:_velocity)
                iter *= _timeStepInv;

            /* Update old positions */
            _positionOld = position;

            this->CalculateMass();
            this->CalculateMomentum();

            /* Update old momentum */
            _xiMomentumOld = _xiMomentum;

            /* No force to be appied */
            fill(_force.begin(),_force.end(), 0.0);

            /* Write log */
            if ( _stepNum % _logFrequency == 0)
                this->WriteLogs( );

            if (!_bMaster)
            {
                _binLambdaOld = transferedLambda;
                fill(_xiBoundaryForceOld.begin(), _xiBoundaryForceOld.end(),0.0);
            }
            return _force;
        }
    }

    /* Compute atom velocities v = (r_new-r_old) / dt */
    transform( position.begin(), position.end(), _positionOld.begin(), _velocity.begin(), std::minus<double>() );
    for (auto &iter:_velocity)
        iter *= _timeStepInv;

    /* Update old positions */
    _positionOld = position;

    this->CalculateMass();
    this->CalculateMomentum();

    /* Compute the system force as f = dp / dt */
    transform( _xiMomentum.begin(), _xiMomentum.end(), _xiMomentumOld.begin(), _xiTotalForce.begin(), std::minus<double>() );
    for (auto &iter:_xiTotalForce)
        iter *= _timeStepInv;

    /* Update old momentum */
    _xiMomentumOld = _xiMomentum;


    try
    {
        if (!_bMaster) _binLambda = transferedLambda;
        this->CalculateLambda();
    }
    catch (const out_of_range &e)
    {
        throw runtime_error( e.what() );
    }

    /* Write log */
    if ( _stepNum % _logFrequency == 0)
    {
        this->WriteLogs( );
    }
    /* Write histograms */
//    if ( _stepNum % _histFrequency == 0)
//    {
//        this->WriteHistograms( to_string(_stepNum) );
//    }

    return _force;
}



void ABF::CalculateMass()
{
    int iD, jD, iN;
    fill( _xiMassInv.begin(), _xiMassInv.end(), 0.0 );
    /* Compute abf Z mass weighted metric tensor */
    for ( iD = 0; iD < _ldim; ++iD )
        for ( jD = 0; jD < _ldim; ++jD )
            for ( iN = 0; iN < _nAtoms; ++iN )
            {
                _xiMassInv.at(iD*_ldim+jD) +=
                        ( _xiJac.at(iD*_hdim + iN*3 + 0) * _xiJac.at(jD*_hdim + iN*3 + 0) +
                          _xiJac.at(iD*_hdim + iN*3 + 1) * _xiJac.at(jD*_hdim + iN*3 + 1) +
                          _xiJac.at(iD*_hdim + iN*3 + 2) * _xiJac.at(jD*_hdim + iN*3 + 2) )
                        / _mass.at(iN);
            }

    /* xi mass tensor is the inverse of xiMassInv */
    if (_ldim == 1)
    {
        if ( _xiMassInv.at(0) < _epsilon )
            cerr<<__func__<<" >> xiMassInv is singular!"<<endl;
        else
            _xiMass.at(0) =  1.0/_xiMassInv.at(0);
    }
    else if (_ldim == 2)
    {
        double det = _xiMassInv.at(0) * _xiMassInv.at(3) - _xiMassInv.at(1) * _xiMassInv.at(2);
        if ( fabs(det) < _epsilon )
            cerr<<__func__<<" >> xiMassInv is singular!"<<endl;
        else
        {
            _xiMass.at(0) =  _xiMassInv.at(3) / det;
            _xiMass.at(1) = -_xiMassInv.at(1) / det;
            _xiMass.at(2) = -_xiMassInv.at(2) / det;
            _xiMass.at(3) =  _xiMassInv.at(0) / det;
        }
    }
    else if (_ldim == 3)
    {
        double det = _xiMassInv.at(1)*_xiMassInv.at(5)*_xiMassInv.at(6) + _xiMassInv.at(2)*_xiMassInv.at(3)*_xiMassInv.at(7)
                   + _xiMassInv.at(0)*_xiMassInv.at(4)*_xiMassInv.at(8) - _xiMassInv.at(2)*_xiMassInv.at(4)*_xiMassInv.at(6)
                   - _xiMassInv.at(0)*_xiMassInv.at(5)*_xiMassInv.at(7) - _xiMassInv.at(1)*_xiMassInv.at(3)*_xiMassInv.at(8);

        // to avoid zero division in singular matrices
        if ( fabs(det) < _epsilon)
        {
            cerr<<__func__<<" >> xiMassInv is singular!"<<endl;
        }
        else
        {
            _xiMass.at(0) = -_xiMassInv.at(5)*_xiMassInv.at(7) + _xiMassInv.at(4)*_xiMassInv.at(8);
            _xiMass.at(1) =  _xiMassInv.at(2)*_xiMassInv.at(7) - _xiMassInv.at(1)*_xiMassInv.at(8);
            _xiMass.at(2) = -_xiMassInv.at(2)*_xiMassInv.at(4) + _xiMassInv.at(1)*_xiMassInv.at(5);
            _xiMass.at(3) =  _xiMassInv.at(5)*_xiMassInv.at(6) - _xiMassInv.at(3)*_xiMassInv.at(8);
            _xiMass.at(4) = -_xiMassInv.at(2)*_xiMassInv.at(6) + _xiMassInv.at(0)*_xiMassInv.at(8);
            _xiMass.at(5) =  _xiMassInv.at(2)*_xiMassInv.at(3) - _xiMassInv.at(0)*_xiMassInv.at(5);
            _xiMass.at(6) = -_xiMassInv.at(4)*_xiMassInv.at(6) + _xiMassInv.at(3)*_xiMassInv.at(7);
            _xiMass.at(7) =  _xiMassInv.at(1)*_xiMassInv.at(6) - _xiMassInv.at(0)*_xiMassInv.at(7);
            _xiMass.at(8) = -_xiMassInv.at(1)*_xiMassInv.at(3) + _xiMassInv.at(0)*_xiMassInv.at(4);

            det = 1.0/det;
            for(auto &item: _xiMass)
                item *= det;
        }
    }
    else
    {
        // FIXIT: complete and check.
        // Compute Cholesky factorization of xiMassInv = L*L' ( xiMassInv is overwritten by L )
        int ierr = LAPACKE_dpotrf( LAPACK_COL_MAJOR, 'L', _ldim, _xiMassInv.data(), _ldim );
        if (ierr == 0)
        {
            // xiMass = I
            fill(_xiMass.begin(), _xiMass.end(), 0.0);
            for ( iD = 0; iD < _ldim; ++iD )
                _xiMass.at(iD*_ldim+iD) = 1.0;

            // Solve L*L'*xiMass = I
            ierr = LAPACKE_dpotrs( LAPACK_COL_MAJOR, 'L', _ldim, _ldim, _xiMassInv.data(),
                                   _ldim, _xiMass.data(), _ldim );

            if (ierr < 0)
                cerr<<__func__<<" >> The "<<ierr<<"-th parameter of LAPACKE_dpotrs had an illegal value."<<endl;
            else if (ierr > 0)
                cerr<<__func__<<" >> Unexpexted error in LAPACKE_dpotrs ( info ="<<ierr<<")."<<endl;
        }
        else if (ierr < 0)
            cerr<<__func__<<" >> The "<<ierr<<"-th parameter of LAPACKE_dpotrf had an illegal value."<<endl;
        else if (ierr > 0)
            cerr<<__func__<<" >> The Mass matrix (the leading minor of order "<<ierr<<") is not positive definite."<<endl;
    }
}


void ABF::CalculateMomentum()
{
    /* Compute xi velocity */
    for (int iD = 0; iD < _ldim; ++iD)
        _xiVelocity.at(iD) = inner_product(_velocity.begin(), _velocity.end(),
                                            _xiJac.begin()+iD*_hdim, 0.0 );
    /* compute xi momentum */
    for (int iD = 0; iD < _ldim; ++iD)
        _xiMomentum.at(iD) =  inner_product(_xiVelocity.begin(), _xiVelocity.end(),
                                            _xiMass.begin()+iD*_ldim, 0.0 );

//    cblas_dgemv( CblasRowMajor, CblasNoTrans, _ldim, _hdim, 1.0, _xiJac.data(), _hdim,
//                 _velocity.data(), 1, 0.0, _xiVelocity.data(), 1 );
//    cblas_dgemv( CblasRowMajor, CblasNoTrans, _ldim, _ldim, 1.0, _xiMass.data(), _ldim,
//                 _xiVelocity.data(), 1, 0.0, _xiMomentum.data(), 1 );
}


void ABF::CalculateLambda()
throw(out_of_range)
{
    int iD, jD;

    /* Find xi bin index */
    for ( iD = 0; iD < _ldim; ++iD )
    {
        if (_xiVal.at(iD) < _xiMin.at(iD) || _xiMax.at(iD) < _xiVal.at(iD))
        {
            cerr<<__func__<<" >> CV value ( "<<_xiVal.at(iD)<<" ) is out of range."<<endl;
            throw out_of_range("CV value is out of range of the histogram!");
        }
        _xiBin.at(iD) = floor( ( (_xiVal.at(iD)-_xiMin.at(iD)) / _xiRange.at(iD) ) * (double)_nBins.at(iD) );
        if (_xiBin.at(iD) == _nBins.at(iD))
            _xiBin.at(iD) -=  1;
    }

    /* Load the histogram information of the xi bin */
    int idx = inner_product(_xiBin.begin(), _xiBin.end(), _accNumBins.begin(), 0 );

    if (_bMaster)
    {
        _binCounts = _histCounts.at( idx );

        /* Calculate the ABF force to be applied on the molecule */
        if (_binCounts >= _rampSample && _bApplyABF )
        {
            /* Runtime average of bin forces */
            for ( iD = 0; iD < _ldim; ++iD)
                _binLambda.at(iD) = -this->ramp(_binCounts) * _histForces.at(idx*_ldim+iD) /  double(_binCounts);
        }
        else
            fill(_binLambda.begin(),_binLambda.end(), 0.0);

        /* Calculate xiBoundaryForce */
        this->CalculateBoundaryForces();

        /* Calculate xiForces as the sum of ABF forces and boundary forces, f_xi = f_abf + f_bnd */
        for ( iD = 0; iD < _ldim; ++iD)
            _xiForce.at(iD) = _binLambda.at(iD) + _xiBoundaryForce.at(iD);

        /* FIXIT: This is only to temporarily disable forces on distances */
        if (_ldim > 2) _xiForce.at(_ldim) = 0.0;

        /* Trnasfer xiForces to individual atoms*/
        fill(_force.begin(), _force.end(), 0.0);
        for (iD = 0; iD < _ldim; ++iD)
            for (jD = 0; jD < _hdim; ++jD)
                _force.at(jD) += _xiForce.at(iD) * _xiJac.at(iD*_hdim+jD);
        //    cblas_dgemv( CblasRowMajor, CblasTrans, _ldim, _hdim, 1.0, _xiJac.data(), _hdim,
        //                 _xiForce.data(), 1, 0.0, _force.data(), 1);

    }
    else
    {
        fill(_xiBoundaryForce.begin(), _xiBoundaryForce.end(), 0.0);
        fill(_force.begin()          , _force.end()          , 0.0);
    }

    /* Update histogram of samples (+=1) */
    ++_histCounts.at( idx );

    /* Update force histogram */
    for (iD = 0; iD < _ldim; ++iD)
        _histForces.at(idx*_ldim+iD) += _xiTotalForce.at(iD) - _binLambdaOld.at(iD) - _xiBoundaryForceOld.at(iD);

    /* Update old forces */
    _binLambdaOld       = _binLambda;
    _xiBoundaryForceOld = _xiBoundaryForce;
}

void ABF::CalculateBoundaryForces()
{
    if(_bHarmonic)
    {
        for (int iD = 0; iD < _ldim; iD++)
        {
            _diffLower    = _lowerWall.at(iD) - _xiVal.at(iD);
            _diffUpper    = _upperWall.at(iD) - _xiVal.at(iD);
            _xiBoundaryForce.at(iD) = _kWall *
                    ( _diffLower > 0.0 ? _diffLower : _diffUpper < 0.0 ? _diffUpper : 0.0 );
        }
    }
    else if(_bCorral)
    {
        _corral.CalculateCorral( _xiVal, _xiBoundaryForce );
    }
    else
        fill( _xiBoundaryForce.begin(), _xiBoundaryForce.end(), 0.0 );
}


void ABF::BuildHistograms()
{
    int iD, jD;
    fill(_accNumBins.begin(), _accNumBins.end(), 1);
    for( iD = 0; iD < _ldim; ++iD )
        for( jD = _ldim-1; jD > iD; --jD )
            _accNumBins.at(iD) *= _nBins.at(jD);

    _totalNumBins = _accNumBins.at(0)*_nBins.at(0);
    cout<<"accNBINS = [";
    for (auto iter : _accNumBins)
        cout<<" "<<iter;
    cout<<" ]"<<endl;
    _histCounts.resize(_totalNumBins      , 0  );
    _histForces.resize(_totalNumBins*_ldim, 0.0);
}



void ABF::LoadHistograms()
throw(runtime_error)
{
    int iD = 0;
    int iN = 0;

    if(_histCounts.size() == 0 || _histForces.size() == 0)
        throw runtime_error("The histogram has to be built before loading data from a file into it!");

    /* Number of bins in each direction of histograms */
    ifstream inHistFile( "histNumBins_"+_coreName+".histi" );
    if (!inHistFile.is_open())
        throw runtime_error("Error in loading the number-of-bin histogram file (histNumBins)!");
    istream_iterator<int> inNumBins (inHistFile);
    for (iD = 0; iD < _ldim; ++iD)
    {
        if (_nBins.at(iD) != *inNumBins)
            throw runtime_error("Number of bins of provided histogram does not match with the one set in the program!");
        ++inNumBins;
    }
    inHistFile.close();

    /* Number of samples in each bin of histogram */
    inHistFile.open( "histCounts_"+_coreName+".histi" );
    if (!inHistFile.is_open())
        throw runtime_error("Error while openning the count histogram file (histCounts)!");
    istream_iterator<double> inCounts (inHistFile);
    for (iN = 0; iN < _totalNumBins; ++iN)
    {
        _histCounts.at(iN) = *inCounts;
        ++inCounts;
    }
    inHistFile.close();

    /* Accumulated forces in each bin of histogram */
    inHistFile.open("histForces_"+_coreName+".histi" );
    if (!inHistFile.is_open())
        throw runtime_error("Error in loading the force histogram file (histForces)!");
    istream_iterator<double> inForces (inHistFile);
    for (iN = 0; iN < _totalNumBins; ++iN)
        for (iD = 0; iD < _ldim; ++iD)
        {
            _histForces.at(_ldim*iN+iD) = *inForces;
            ++inForces;
        }
    inHistFile.close();

}


double ABF::ramp(int nSamples)
{
    if (nSamples >= _fullSample)
        return 1.0;
    else if (nSamples < _rampSample)
        return 0.0;
    else
        return double(nSamples-_rampSample)/double(_fullSample-_rampSample);
}


void ABF::WriteHistograms(string postfix)
throw(runtime_error)
{
    /* Number of bins in each direction of histograms */
    _histFile.open( "histNumBins_"+_coreName+"_"+postfix+".histo" );
    if (!_histFile)
        throw runtime_error("Error while openning the number-of-bin histogram file (histNumBins)!");
    ostream_iterator<int> outNumBins (_histFile,"\n");
    copy(_nBins.begin(), _nBins.end(), outNumBins);
    _histFile.close();

    /* Number of samples in each bin of histogram */
    _histFile.open( "histCounts_"+_coreName+"_"+postfix+".histo" );
    if (!_histFile)
        throw runtime_error("Error while openning the count histogram file (histCounts)!");
    ostream_iterator<int> outCounts (_histFile,"\n");
    copy(_histCounts.begin(), _histCounts.end(), outCounts);
    _histFile.close();

    /* Accumulated forces in each bin of histogram */
    _histFile.open( "histForces_"+_coreName+"_"+postfix+".histo" );
    if (!_histFile)
        throw runtime_error("Error while openning the force histogram file (histForces)!");
    _histFile<<scientific<<setprecision(15);
    ostream_iterator<double> outForce (_histFile,"\n");
    copy(_histForces.begin(), _histForces.end(), outForce);
    _histFile.close();
}


void ABF::WriteLogs()
throw(runtime_error)
{

    /* Open the log file */
//    _logFile.open( _logName, ios_base::out|ios_base::app );

    if (!_logFile)
        throw runtime_error("Error while openning the log file!");

    /* Step number */
    _logFile<<_stepNum;

    /* CV values */
    for (auto iter : _xiVal)
        _logFile<<" "<<scientific<<setprecision(15)<<iter;

    /* CV forces */
    for (auto iter : _xiForce)
        _logFile<<" "<<scientific<<setprecision(15)<<iter;

    /* CV Boundary forces */
    for (auto iter : _xiBoundaryForce)
        _logFile<<" "<<scientific<<setprecision(15)<<iter;

    /******************** TEMPORARY ********************/
//    /* Forces */
//    _histFile.open( "force.txt",ios_base::app);
//    if (!_histFile)
//        throw runtime_error("Error while openning the number-of-bin histogram file (histNumBins)!");
//    ostream_iterator<double> outForce (_histFile,"\t");
//    copy(_force.begin(), _force.end(), outForce);
//    _histFile<<endl;
//    _histFile.close();

//    /* Jacobian */
//    _histFile.open( "jac.txt",ios_base::app);
//    if (!_histFile)
//        throw runtime_error("Error while openning the number-of-bin histogram file (histNumBins)!");
//    ostream_iterator<double> outJac (_histFile,"\t");
//    copy(_xiJac.begin(), _xiJac.end(), outJac);
//    _histFile<<endl;
//    _histFile.close();
    /******************** END OF TEMPORARY ********************/

    /* Next line */
    _logFile<<endl;

    /* Close the log file */
//    _logFile.close();
}
