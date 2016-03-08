/*=========================================================================
 * File         : SMD.cpp
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
#include"SMD.h"
static const double _epsilon = numeric_limits<double>::epsilon(); ///< Machine precision


SMD::SMD(bool bApplySMD, bool bDistance, bool bWithReturn, int ldim, int sdim, int nAtoms, double timeStep, \
         double pullingRate, double pullingConstant, double pullingMaxLength,
         int logFrequency,  string coreName,
         double kWall, const vector<double> &lowerWall, const vector<double> &upperWall,
         double xiMin, double xiMax)
    :
      _bApplySMD(bApplySMD),
      _bDistance(bDistance),
      _sdim(sdim),
      _ldim(ldim),
      _nAtoms(nAtoms),
      _logFrequency(logFrequency),
      _timeStep(timeStep),
      _kWall(kWall),
      _lowerWall(lowerWall),
      _upperWall(upperWall),
      _xiMin(xiMin),
      _xiMax(xiMax),
      _bWithReturn(bWithReturn),
      _pullingRate(pullingRate),
      _pullingConstant(pullingConstant),
      _pullingMaxLength(pullingMaxLength),
      _coreName(coreName)

{
    _bHarmonic = _kWall > 0.0 ? true : false;
    /* Calculate xi range */
    _xiRange = _xiMax - _xiMin;
    this->Update();
}



SMD::~SMD()
{
    _logFile.close();
}


void SMD::Update( )
throw( invalid_argument )
{
    if ( this->CheckSettings() )
    {
        cerr<<__func__<<" >> CheckSettings went wrong."<<endl;
        throw invalid_argument("One or more inputs are wrong.");
    }
    _hdim = _sdim*_nAtoms;
    _pullingDisplacement = _pullingRate*_timeStep;
    _bdim = _bDistance ? _ldim+1 : _ldim;

    _xiForce        .resize(_bdim, 0.0);
    _xiBoundaryForce.resize(_bdim, 0.0);
    _positionOld    .resize(_hdim, 0.0);
    _displacement   .resize(_hdim, 0.0);
    _force          .resize(_hdim, 0.0);

    _logName   = "log_" +_coreName+".txt";


    /* Create the log file */
    _logFile.open( _logName );

    /* Write the header*/
    _logFile<<"StepNumber";
    for (int iD = 0; iD < _bdim; ++iD)
        _logFile<<" xiValue"<<iD;
    _logFile<<" Distance";
    for (int iD = 0; iD < _bdim; ++iD)
        _logFile<<" xiForce"<<iD;
    for (int iD = 0; iD < _bdim; ++iD)
        _logFile<<" xiBndForce"<<iD;
    _logFile<<" xiValueEnd Work WorkToral"<<endl;


}


bool SMD::CheckSettings( )
{
    bool bError = false;

    if (_ldim != 1)
    {
        cerr<<__func__<<" >>  CV's dimension has to be one (ldim = "<<_ldim<<")."<<endl;
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
    if (_xiRange < _pullingMaxLength)
    {
        cerr<<__func__<<" >>  xi range("<<_xiRange<<") is less than pulling length("<<_pullingMaxLength<<")."<<endl;
        bError = true;
    }
    return bError;
}



vector<double> &SMD::CalculateForce( const vector<double> &position,
                                     const vector<double> &xiVal,
                                     const vector<double> &xiJac,
                                     double cdist)
{
    _xiVal = xiVal;\
    _xiJac = xiJac;
    _cdist = cdist;

    cout<<"------------------------------------------------------------"<<endl;
    cout<<"nSteps = "<<_nSteps<<endl;

    if (++_nSteps == 1)
    {
        /* Attach one end of sring */
        _xiPulling = _xiVal.at(0);

        /* Update old positions */
        _positionOld = position;

        /* Write log */
        this->WriteLogs( );

        return _force;
    }

    /* Calculate xiBoundaryForce */
    this->CalculateBoundaryForces();
    cout<<"force_bnd = ";
    for (auto iter : _xiBoundaryForce)
        cout<<" "<<iter;
    cout<<endl;

    /* Compute CV's forces, F = -k*dxi  and add boundary forces  */
    _xiForce.at(0) = _pullingConstant*(_xiPulling-_xiVal.at(0)) + _xiBoundaryForce.at(0);
    if (_bDistance)
        _xiForce.at(1) = _xiBoundaryForce.at(1);
    cout<<"force = ";
    for (auto iter : _xiForce)
        cout<<" "<<iter;
    cout<<endl;

    /* Trnasfer xiForces to individual atoms*/
    cblas_dgemv( CblasRowMajor, CblasTrans, _ldim, _hdim, 1.0, _xiJac.data(), _hdim,
                 _xiForce.data(), 1, 0.0, _force.data(), 1);
//    _force = xiJac;
//    cblas_dscal(_hdim, _xiForce.at(0), _force.data(), 1);

    /* Calculate constant rate pulling work, W = F.dx */
    transform( position.begin(), position.end(), _positionOld.begin(), _displacement.begin(), std::minus<double>() );
    _work = cblas_ddot(_hdim, _force.data(), 1, _displacement.data(), 1);
    _workTotal += _work;

    /* Update pulling length */
    if (_bWithReturn)
    {
        if( (_pullingLength < _pullingMaxLength && _bForwardMode) ||
            (_pullingLength > 0.0 && _bBackwardMode) )
        {
            _pullingLength += _pullingDisplacement;
            _xiPulling     += _pullingDisplacement;
        }
        else
        {
            _pullingDisplacement *= -1.0;
            swap(_bForwardMode,_bBackwardMode);
        }
    }
    else if( _pullingLength < _pullingMaxLength )
    {
        _pullingLength += _pullingDisplacement;
        _xiPulling     += _pullingDisplacement;
    }

    /* Write log */
    if ( _nSteps % _logFrequency == 0)
        this->WriteLogs( );

    cout<<"============================================================"<<endl;
    return _force;
}


void SMD::CalculateBoundaryForces()
{
    if(_bHarmonic)
    {
        for (int iD = 0; iD < _bdim; iD++)
        {
            _diffLower    = _lowerWall.at(iD) - _xiVal.at(iD);
            _diffUpper    = _upperWall.at(iD) - _xiVal.at(iD);
            _xiBoundaryForce.at(iD) = _kWall *
                    ( _diffLower > 0.0 ? _diffLower : _diffUpper < 0.0 ? _diffUpper : 0.0 );
        }
    }
    else if(_bCorral)
        // TODO: implement corral
        fill( _xiBoundaryForce.begin(), _xiBoundaryForce.end(), 0.0 );
    else
        fill( _xiBoundaryForce.begin(), _xiBoundaryForce.end(), 0.0 );
}


void SMD::WriteLogs()
throw(runtime_error)
{

    /*IMPORTANT: do not forget to change the header in the Update() function*/

    /* Open the log file */
//    _logFile.open( _logName, ios_base::out|ios_base::app );

    if (!_logFile)
        throw runtime_error("Error while openning the log file!");

    /* Step number */
    _logFile<<_nSteps;

    /* CV values */
    for (auto iter : _xiVal)
        _logFile<<" "<<scientific<<setprecision(15)<<iter;

    /* Distance to the manifold */
    _logFile<<" "<<scientific<<setprecision(15)<<_cdist;

    /* CV forces */
    for (auto iter : _xiForce)
        _logFile<<" "<<scientific<<setprecision(15)<<iter;

    /* CV Boundary forces */
    for (auto iter : _xiBoundaryForce)
        _logFile<<" "<<scientific<<setprecision(15)<<iter;

    /* CV value of other end of the spring */
    _logFile<<" "<<scientific<<setprecision(15)<<_xiPulling;

    /* Work at eatch time step */
    _logFile<<" "<<_work<<" "<<_workTotal;

    /* Next line */
    _logFile<<endl;

    /* Close the log file */
//    _logFile.close();

    /******************** TEMPORARY ********************/
    /* Forces */
    _histFile.open( "force.txt",ios_base::app);
    if (!_histFile)
        throw runtime_error("Error while openning force.txt!");
    ostream_iterator<double> outForce (_histFile,"\t");
    copy(_force.begin(), _force.end(), outForce);
    _histFile<<endl;
    _histFile.close();

    /* Jacobian */
    _histFile.open( "jac.txt",ios_base::app);
    if (!_histFile)
        throw runtime_error("Error while openning jac.txt!");
    ostream_iterator<double> outJac (_histFile,"\t");
    copy(_xiJac.begin(), _xiJac.end(), outJac);
    _histFile<<endl;
    _histFile.close();

    /******************** END OF TEMPORARY ********************/
}
