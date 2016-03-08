/*=========================================================================
 * File         : MPSandABF.cpp
 * Module       : Swig
 * Copyright    : (C)opyright 2013
 *                See COPYRIGHT statement in top level directory.
 * Authors      : Behrooz Hashemian
 * Created      : Thu Oct 17, 2013  11:50AM
 * Last modified: Sat Feb 22, 2014  02:45PM
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
#include "MPSandABF.h"

MPSandABF::MPSandABF(int nPartitions,
                     int bApplyABF, int nAtoms, int sdim, double timeStep,
                     int logFrequency, int histFrequency, char *coreName,
                     char *xRefFilename, int nInvolvedAtoms, int *involvedAtoms,
                     int ldim, int bDistance,
                     char *xiNodesFilename, char *xNodesFilename,
                     char *xiSeedsFilename,
                     char *xiMarkersFilename,
                     int bLoadHistogram, double *xiMin, double *xiMax, int *nBins,
                     double gammaPU, double spacingPU, double tol0PU)
    : _bDistance(bDistance),
      _ldim(ldim),
      _hdim(sdim*nAtoms),
      _nAtoms(nAtoms),
      _histFrequency(histFrequency),
      _nPartitions(nPartitions),
      _vsandcv(nPartitions),
      _vabf(nPartitions)
{

    int nLines1 = 0;
    int nLines2 = 0;
    int nCols1  = 0;
    int nCols2  = 0;
    int abfDim  = _bDistance ? _ldim+1 : _ldim;

    // boundary potentials, which sould not be active in multi-patch
    double         kWall = 0;
    vector<double> lowerWall(abfDim, -numeric_limits<double>::max());
    vector<double> upperWall(abfDim, +numeric_limits<double>::max());

    //----------------
    vector<int>    vInvolvedAtoms(involvedAtoms,involvedAtoms+nInvolvedAtoms);
    vector<int>    vNumBins(nBins,nBins+abfDim);
    vector<double> xRef;
    vector<double> xiNodes;
    vector<double> xNodes;
    vector<double> xiSeeds;
    vector<double> xMarkers;
    vector<double> xiMarkers;
    PointSetIO     reader;

    //______________________________________________________________________________
    // Read xRef
    reader.ReadPointSet(xRefFilename, xRef, nLines1, nCols1);
    if (nLines1 < 1)
        cerr<<"xRefFilename("<<xRefFilename<<") is empty."<<endl;
    else if (nLines1 > 1)
        cerr<<"xRefFilename("<<xRefFilename<<") includes more than one configuration."<<endl;

    _alignment = make_shared<ProcrustesSuperimposition>(sdim,xRef,vInvolvedAtoms);

    ofstream outFile("xi_range.txt");
    cout<<"Number of partitions = "<<nPartitions<<endl;
    for ( int iP = 0; iP < nPartitions; ++iP )
    {
        cout<<"****************************************"<<endl;
        cout<<"****Initializing partition ["<<iP<<"]..."<<endl;
        xiNodes.clear();
        xNodes.clear();
        xiSeeds.clear();
        xiMarkers.clear();

        string postfix("_part"+to_string(iP)+".txt");
        //______________________________________________________________________________
        // Read xiNodes
        reader.ReadPointSet( xiNodesFilename+postfix, xiNodes , nLines1, nCols1 );
        if (nLines1 < 1)
            cerr<<"xiNodesFilename("<<xiNodesFilename<<") is empty."<<endl;
        if (nCols1 != _ldim)
            cerr<<"Dimension of nodes points does not match with ldim ("<<nCols1<<"!="<<_ldim<<")."<<endl;
        double spacing = xiNodes.back();
        xiNodes.erase(xiNodes.end()-_ldim,xiNodes.end());
        --nLines1;
        cout<<"spacing = "<<spacing<<endl;

        //______________________________________________________________________________
        // Read xNodes
        reader.ReadPointSet( xNodesFilename+postfix, xNodes  , nLines2, nCols2 ) ;
        if (nLines2 < 1)
            cerr<<"xNodesFilename("<<xNodesFilename<<") is empty."<<endl;
        if (nLines2 != nLines1)
            cerr<<"Number of low and high dimensional points does not match("
               <<nLines1<<"!="<<nLines2<<")."<<endl;

        //______________________________________________________________________________
        // Read xiSeeds
        reader.ReadPointSet( xiSeedsFilename+postfix, xiSeeds , nLines2, nCols2 );
        if (nLines2 < 1)
            cerr<<"xiSeedsFilename("<<xiSeedsFilename<<") is empty."<<endl;
        if (nCols2 != nCols1)
            cerr<<"Dimension of seeds and nodes points does not match("
               <<nCols2<<"!="<<nCols1<<")."<<endl;

        //______________________________________________________________________________
        // Read xiMarkers
        reader.ReadPointSet( xiMarkersFilename+postfix, xiMarkers , nLines2, nCols2 );
        if (nLines2 < 1)
            cerr<<"xiMarkersFilename("<<xiMarkersFilename<<") is empty."<<endl;
        if (nCols2 != nCols1)
            cerr<<"Dimension of markers and nodes points does not match("
               <<nCols2<<"!="<<nCols1<<")."<<endl;

        //______________________________________________________________________________
        // Find the min and max of xi nodes

        vector<double> vXiMin( abfDim, +numeric_limits<double>::max() );
        vector<double> vXiMax( abfDim, -numeric_limits<double>::max() );
        size_t iD;
        for ( auto iter = xiNodes.begin(); iter != xiNodes.end(); )
            for ( iD = 0; iD < _ldim; ++iD, ++iter )
            {
                vXiMin.at(iD) = min( *iter, vXiMin.at(iD) );
                vXiMax.at(iD) = max( *iter, vXiMax.at(iD) );
            }
        if (_bDistance)
        {
            vXiMin.back() = xiMin[_ldim];
            vXiMax.back() = xiMax[_ldim];
        }
        cout<<"xi min = [";
        for( auto iter : vXiMin )
        {
            cout<<iter<<"\t";
            outFile<<iter<<"\t";
        }
        cout<<"]"<<endl;
        cout<<"xi max = [";
        for( auto iter : vXiMax )
        {
            cout<<iter<<"\t";
            outFile<<iter<<"\t";
        }
        cout<<"]"<<endl;
        outFile<<endl;


        //______________________________________________________________________________
        // Initialize SandCV
        _vsandcv.at(iP) = make_shared<SandCV>( sdim, _ldim, _hdim, spacing, xRef, vInvolvedAtoms, xiNodes, xNodes, xiSeeds, xiMarkers);
        cout<<"SandCV no. ["<<iP<<"] is initialized."<<endl;
        xMarkers = _vsandcv.at(iP)->GetMarkers();
        _xMarkers.insert(_xMarkers.end(), xMarkers.begin(), xMarkers.end());
        _mark2part.resize(_mark2part.size()+xiMarkers.size()/_ldim, iP);

        //______________________________________________________________________________
        // Initialize ABF
        _vabf.at(iP) = make_shared<ABF>(bApplyABF, bLoadHistogram, abfDim, sdim, _nAtoms, timeStep, logFrequency, histFrequency,
                                        kWall, lowerWall, upperWall, "P"+to_string(iP)+"_"+coreName, vXiMin, vXiMax, vNumBins );
        cout<<"ABF no. ["<<iP<<"] is initialized."<<endl;
        cout<<"----------------------------------------------------------------------"<<endl;
    }
    outFile.close();

    cout<<"mark2part size = "<<_mark2part.size()<<endl;

    int    nMarkers  = _xMarkers.size()/_hdim;
    _betaPU          = gammaPU/spacingPU/spacingPU;
    double rangePU   = sqrt(-log(tol0PU)/_betaPU);

    ANNSearch markerNeighborSearch(rangePU, _hdim, _xMarkers);
    vector<double> aMarker(_hdim);
    _markerNeighbors.resize(nMarkers);
    for ( int iM = 0; iM < nMarkers; ++iM )
    {
        copy(_xMarkers.begin()+iM*_hdim, _xMarkers.begin()+(iM+1)*_hdim, aMarker.begin());
        _markerNeighbors.at(iM) = markerNeighborSearch.FindNearestNeighborsFR( aMarker );
    }
    cout<<"MARKER SEARCH: "<<endl<<
          "\tbeta  = "<<_betaPU<<endl<<
          "\trange = "<<rangePU<<endl;
    _closestMarkerSearch = make_shared<ANNSearch>(1, _hdim, _xMarkers);

    // ABF forces of each partition
    _forceABF.resize(_nPartitions);
    for( auto &iter : _forceABF )
        iter.resize(_hdim, 0.0);
    _diff       .resize(_hdim);
    _normDiff   .resize(_nPartitions);
    _psi        .resize(_nPartitions);
    _forceTotal .resize(_hdim);
    _lambdaSlave.resize(_ldim);
    _dPhi       .resize(_ldim*_ldim);
    _dPhiOverlap.resize(_nPartitions*_ldim*_ldim);
    _xiOverlap  .resize(_nPartitions*_ldim);
    _fileXiOverlap  .open("xi_overlap.txt");
    _filePsiOverlap .open("psi_overlap.txt");
    _fileDPhiOverlap.open("dphi_overlap.txt");
    _fileDiffOverlap.open("diff_overlap.txt");
}

MPSandABF::~MPSandABF()
{
    _fileXiOverlap.close();
}


bool compareParts(pair<int,double> P1, pair<int,double> P2)
{
    return P1.second > P2.second;
}


bool isSmall(pair<int,double> P1)
{
    return P1.second < 1e-8;
}


int MPSandABF::TransferLambda()
{
    if (_ldim == 1)
    {
        _lambdaSlave.at(0) *= _lambdaMaster.at(0)
                * inner_product(_dxMaster.begin(), _dxMaster.end(), _dxSlave.begin() , 0.0)
                / inner_product(_dxMaster.begin(), _dxMaster.end(), _dxMaster.begin(), 0.0);
        return 0;
    }
    else if (_ldim == 2)
    {
        vector<double> temp1(3);
        temp1.at(0) = inner_product(_dxMaster.begin()      , _dxMaster.begin()+_hdim, _dxMaster.begin()      , 0.0);
        temp1.at(1) = inner_product(_dxMaster.begin()      , _dxMaster.begin()+_hdim, _dxMaster.begin()+_hdim, 0.0);
        temp1.at(2) = inner_product(_dxMaster.begin()+_hdim, _dxMaster.end()        , _dxMaster.begin()+_hdim, 0.0);

        vector<double> temp2(4);
        temp2.at(0) = inner_product(_dxMaster.begin()      , _dxMaster.begin()+_hdim, _dxSlave.begin()      , 0.0);
        temp2.at(1) = inner_product(_dxMaster.begin()      , _dxMaster.begin()+_hdim, _dxSlave.begin()+_hdim, 0.0);
        temp2.at(2) = inner_product(_dxMaster.begin()+_hdim, _dxMaster.end()        , _dxSlave.begin()       , 0.0);
        temp2.at(3) = inner_product(_dxMaster.begin()+_hdim, _dxMaster.end()        , _dxSlave.begin()+_hdim, 0.0);

        double det = temp1.at(0)*temp1.at(2) - temp1.at(1)*temp1.at(1);
        if (fabs(det) < numeric_limits<double>::epsilon())
        {
            cerr<<__func__<<" >> error in transfering CV, det = "<<det<<"."<<endl;
            return 1;
        }
        else
        {
            _dPhi.at(0) = ( temp1.at(2)*temp2.at(0) - temp1.at(1)*temp2.at(2) ) / det;
            _dPhi.at(1) = ( temp1.at(2)*temp2.at(1) - temp1.at(1)*temp2.at(3) ) / det;
            _dPhi.at(2) = (-temp1.at(1)*temp2.at(0) + temp1.at(0)*temp2.at(2) ) / det;
            _dPhi.at(3) = (-temp1.at(1)*temp2.at(1) + temp1.at(0)*temp2.at(3) ) / det;

            _lambdaSlave.at(0) = _lambdaMaster.at(0)*_dPhi.at(0) + _lambdaMaster.at(1)*_dPhi.at(2);
            _lambdaSlave.at(1) = _lambdaMaster.at(0)*_dPhi.at(1) + _lambdaMaster.at(1)*_dPhi.at(3);

            return 0;
        }
    }
    else
    {
        cerr<<"Transfering CV is implemented only for 1D and 2D!"<<endl;
        return 1;
    }
}


void MPSandABF::CalculateForces( int stepNum, double *masses, double *positions, double *forces )
{
//    cout<<"%%%%["<<stepNum<<"]%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
    /* Change the masses and positions array to the STL vectors */
    vector<double> vMasses   ( masses   , masses   +_nAtoms);
    vector<double> vPositions( positions, positions+_hdim  );
    vector<double> aPositions( _alignment->ComputeAlignment(vPositions) );

    ++stepNum;

    const int closestMarker = _closestMarkerSearch->FindNearestNeighbors(aPositions).at(0);
//    cout<<"Closest marker = "<<closestMarker<<endl;

//    ANNSearch checkMarkerSearch(100, _hdim, _xMarkers);
//    vector<int> checkMarkers = checkMarkerSearch.FindNearestNeighbors(aPositions);
//    vector<int> markerParts;
//    for(auto mid:checkMarkers)
//        markerParts.push_back(_mark2part.at(mid));
//    int partCount = count(markerParts.begin(),markerParts.end(),_mark2part.at(closestMarker));
//    cout<<"partition ["<<_mark2part.at(closestMarker)<<"] counts = "<<partCount<<endl;

    const vector<int> markID = _markerNeighbors.at(closestMarker);
//    cout<<"Number of neighbors = "<<markID.size()<<endl;

    set<int> setPartID;
    for ( const auto &mid : markID)
        setPartID.insert(_mark2part.at(mid));

    vector<pair<int,double>> partInfo;
    for (auto &pid:setPartID)
    {
        partInfo.push_back(make_pair(pid,0.0));
    }

    int nParts = partInfo.size();
//    cout<<"Parts ("<<nParts<<") :";
//    for(auto &iter:partInfo)
//        cout<<iter.first<<"  ";
//    cout<<endl;

    /* Compute the CV value of the current configuration and its related Jacobian */
    try
    {
        if (nParts == 0)
            throw(runtime_error("No partition found!"));
        else if (nParts == 1)
        {
            partInfo.front().second         = 1.0;
            _psi.at(partInfo.front().first) = 1.0;
        }
        else
        {
            _Z  = 0.0;
            _wA = 0.0;
            _pid = -1;
            fill(_psi.begin(), _psi.end(), 0.0);
            for ( const auto &mid : markID )
            {
                _pid = _mark2part.at(mid);
                transform(aPositions.begin(), aPositions.end(), _xMarkers.begin()+_hdim*mid,
                          _diff.begin(), std::minus<double>() );
                _wA = exp( -_betaPU * cblas_ddot(_hdim, _diff.data(), 1, _diff.data(), 1 ) );
                _Z            += _wA;
                _psi.at(_pid) += _wA;
            }
            _Z = 1.0 / _Z; // Inverse Z
            for ( auto &iter : _psi )
            {
                iter *= _Z;
            }

            for (auto &iter:partInfo)
                iter.second = _psi.at(iter.first);

            /* Remove small psi */
            partInfo.erase(remove_if(partInfo.begin(),partInfo.end(),isSmall),partInfo.end());

            /* Sort based on psi value */
            sort(partInfo.begin(), partInfo.end(), compareParts );

            nParts = partInfo.size();
//            cout<<"Parts ("<<nParts<<") :"<<endl;
//            for(auto &iter:partInfo)
//                cout<<"["<<iter.first<<"] = "<<iter.second<<endl;
        }

        int nFailedParts = nParts;
        fill(_normDiff.begin()   , _normDiff.end()   , NAN);
        fill(_xiOverlap.begin()  , _xiOverlap.end()  , NAN);
        fill(_dPhiOverlap.begin(), _dPhiOverlap.end(), NAN);
        int pid = -1;
        try
        {
            pid = partInfo.at(0).first;
//            cout<<"******** Master part [ "<<pid<<" ] = "<<_psi.at(pid)<<"("<<partInfo.at(0).second<<")  ********"<<endl;
            _vsandcv.at(pid)->CalculateCV(vPositions, stepNum);

            /* Calculate the ABF forces and return it */
            _forceTotal = _vabf.at(pid)->CalculateForce( stepNum, vMasses, vPositions, _vsandcv.at(pid)->GetValue(), _vsandcv.at(pid)->GetJacobian(), true );

            _lambdaMaster = _vabf.at(pid)->GetLambda();
            _xMaster      = _vsandcv.at(pid)->GetX();
            _dxMaster     = _vsandcv.at(pid)->GetDX();

            copy(_vsandcv.at(pid)->GetValue().begin(), _vsandcv.at(pid)->GetValue().end(),
                         _xiOverlap.begin()+pid*_ldim);
            nFailedParts--;
            for ( int i = 1; i < nParts; ++i )
            {
                pid = partInfo.at(i).first;
//                cout<<"******** Slave part [ "<<pid<<" ] = "<<_psi.at(pid)<<"("<<partInfo.at(i).second<<")  ********"<<endl;
                try
                {
                    _vsandcv.at(pid)->CalculateCV(vPositions, stepNum);

                    /* Transform xiForce evaluated at master part and pass it to ABF */
                    _xSlave  = _vsandcv.at(pid)->GetX();
                    _dxSlave = _vsandcv.at(pid)->GetDX();
                    if (TransferLambda())
                        throw runtime_error("Error in tranforming CVs!");
                    else
                    {
                        _vabf.at(pid)->CalculateForce( stepNum, vMasses, vPositions, _vsandcv.at(pid)->GetValue(), _vsandcv.at(pid)->GetJacobian(), false, _lambdaSlave );
                    }
                    // copy the xi value to be save in overlap file
                    /* xDiff = xM - xS */
                    transform( _xMaster.begin() , _xMaster.end(), _xSlave.begin(), _diff.begin(), std::minus<double>() );
                    _normDiff.at(pid) = inner_product(_diff.begin(), _diff.end(), _diff.begin(), 0.0);
                    copy(_vsandcv.at(pid)->GetValue().begin(), _vsandcv.at(pid)->GetValue().end(),
                         _xiOverlap.begin()+pid*_ldim);
                    copy(_dPhi.begin(), _dPhi.end(), _dPhiOverlap.begin()+pid*_ldim*_ldim);
                    nFailedParts--;
                }
                catch(const runtime_error &e)
                {
                    cout<<"Error catched ("<<e.what()<<") and "<<pid<<" ("<<_psi.at(pid)<<") check if it is important!"<<endl;
                    if (_psi.at(pid) > 1e-7)
                        cerr<<"Error :: ("<<e.what()<<") in slave part ["<<pid<<"] ("<<_psi.at(pid)<<")!"<<endl;
                    else
                        cout<<"Not important!"<<endl;
                }
            }
        }
        catch(const runtime_error &e)
        {
            cerr<<"Error catched ("<<e.what()<<") and "<<pid<<" ("<<_psi.at(pid)<<") is removed!"<<endl;
        }

        if (stepNum%100 == 0 && nParts-nFailedParts > 1)
        {
            for (auto &iter : _xiOverlap)
                _fileXiOverlap<<iter<<"\t";
            _fileXiOverlap<<endl;

            for (auto &iter : _psi)
                _filePsiOverlap<<iter<<"\t";
            _filePsiOverlap<<endl;

            for (auto &iter : _dPhiOverlap)
                _fileDPhiOverlap<<iter<<"\t";
            _fileDPhiOverlap<<endl;

            for (auto &iter : _normDiff)
                _fileDiffOverlap<<iter<<"\t";
            _fileDiffOverlap<<endl;
        }

        /* Copy the calculated forces to pass to the MD engine */
        copy( _forceTotal.begin(), _forceTotal.end(), forces );

        if (stepNum % _histFrequency == 0)
            for (auto &item: _vabf)
            {
                item->WriteHistograms( "current" );
                item->WriteHistograms( to_string(stepNum) );
            }
    }
    catch (const runtime_error &e)
    {
        fill( forces, forces+_hdim, 0.0 );
        cout<<"NO FORCE is applied because "<<e.what()<<endl;
    }

}



void MPSandABF::CalculateUniqueForces( int stepNum, double *masses, double *positions, double *forces )
{
    cout<<"%%%%["<<stepNum<<"]%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
    /* Change the masses and positions array to the STL vectors */
    vector<double> vMasses   ( masses   , masses   +_nAtoms);
    vector<double> vPositions( positions, positions+_hdim  );
    vector<double> aPositions( _alignment->ComputeAlignment(vPositions) );

    //    const vector<int> markID = _markerSearch->FindNearestNeighborsFR( aPositions );
    const vector<int> markID = _closestMarkerSearch->FindNearestNeighbors( aPositions );
    cout<<"knn markers = "<<markID.size()<<endl;

    set<int> setPartID;
    for ( const auto &mid : markID)
        setPartID.insert(_mark2part.at(mid));

    vector<int> partID;
    for (auto &pid:setPartID)
        partID.push_back(pid);

    int nParts = partID.size();
    cout<<"Parts ["<<nParts<<"] :";
    for(auto &iter:partID)
        cout<<iter<<"  ";
    cout<<endl;

    int pidMax;
    /* Compute the CV value of the current configuration and its related Jacobian */
    try
    {
        if (nParts == 0)
            throw(runtime_error("No partition found!"));
        else if (nParts == 1)
        {
            pidMax = *partID.begin();
            _psi.at(pidMax) = 1.0;
        }
        else
        {
            _Z  = 0.0;
            _wA = 0.0;
            _pid = -1;
            fill(_psi.begin(), _psi.end(), 0.0);
            for ( const auto &mid : markID )
            {
                _pid = _mark2part.at(mid);
                transform(aPositions.begin(), aPositions.end(), _xMarkers.begin()+_hdim*mid,
                          _diff.begin(), std::minus<double>() );
                _wA = exp( -_betaPU * cblas_ddot(_hdim, _diff.data(), 1, _diff.data(), 1 ) );
                _Z            += _wA;
                _psi.at(_pid) += _wA;
            }
//            double sum = 0.0;
            _Z = 1.0 / _Z; // Inverse Z
            for ( auto &iter : _psi )
            {
                iter *= _Z;
//                sum += iter;
//                cout<<"iter = "<<iter<<endl;
            }
//            cout<<"Sum psi = "<<sum<<endl;

            pidMax = 0;
            for ( int pid = 1; pid < _nPartitions; ++pid)
                if (_psi.at(pid) > _psi.at(pidMax)) pidMax = pid;
        }

        cout<<"******** part [ "<<pidMax<<" ] = "<<_psi.at(pidMax)<<"  ********"<<endl;
        try
        {
            _vsandcv.at(pidMax)->CalculateCV(vPositions, stepNum);

            /* Calculate the ABF forces and return it */
            _forceABF.at(pidMax) = _vabf.at(pidMax)->CalculateForce( stepNum, vMasses, vPositions, _vsandcv.at(pidMax)->GetValue(), _vsandcv.at(pidMax)->GetJacobian() );
        }
        catch(const runtime_error &e)
        {
            cout<<"Error catched ("<<e.what()<<") and "<<pidMax<<" ("<<_psi.at(pidMax)<<") is removed!"<<endl;
            cerr<<"Error catched ("<<e.what()<<") and "<<pidMax<<" ("<<_psi.at(pidMax)<<") is removed!"<<endl;
        }

        /* Copy the calculated forces to pass to the MD engine */
        copy( _forceABF.at(pidMax).begin(), _forceABF.at(pidMax).end(), forces );
    }
    catch (const runtime_error &e)
    {
        fill( forces, forces+_hdim, 0.0 );
        cout<<"NO FORCE is applied because "<<e.what()<<endl;
    }

}
