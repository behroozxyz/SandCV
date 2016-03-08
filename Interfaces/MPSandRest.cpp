/*=========================================================================
 * File         : MPSandRest.cpp
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
#include "MPSandRest.h"

MPSandRest::MPSandRest(int nPartitions,
                       int nAtoms, int sdim, int logFrequency,
                       char *xRefFilename, int nInvolvedAtoms, int *involvedAtoms,
                       int ldim, int bDistance,
                       char *xiNodesFilename, char *xNodesFilename,
                       char *xiSeedsFilename, char *xiMarkersFilename,
                       char *xCenterFilename, double kBias, double rBias,
                       double *xiMin, double *xiMax,
                       double gammaPU, double spacingPU, double tol0PU)
    : _bDistance(bDistance),
      _nAtoms(nAtoms),
      _ldim(ldim),
      _hdim(sdim*nAtoms),
      _nPartitions(nPartitions),
      _vsandcv(nPartitions),
      _logFrequency(logFrequency),
      _kBias(kBias),
      _rBias(rBias)
{

    int nLines1 = 0;
    int nLines2 = 0;
    int nCols1  = 0;
    int nCols2  = 0;
    int abfDim  = _bDistance ? _ldim+1 : _ldim;

    vector<int>    vInvolvedAtoms(involvedAtoms,involvedAtoms+nInvolvedAtoms);
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

    //______________________________________________________________________________
    // Read the center of harmonic potential
    reader.ReadPointSet(xCenterFilename, _xCenter, nLines1, nCols1);
    if (nLines1 < 1)
        cerr<<"xCenterFilename("<<xCenterFilename<<") is empty."<<endl;
    else if (nLines1 > 1)
        cerr<<"xCenterFilename("<<xCenterFilename<<") includes more than one configuration."<<endl;
    if (nCols1 != _hdim)
        cerr<<"xCenterFilename("<<xCenterFilename<<"):"<<
              "Dimension of zone centers("<<nCols1<<") does not match with hdim ("<<_hdim<<")."<<endl;

    _xCenter = _alignment->ComputeAlignment(_xCenter);


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
    _xCenter     .resize(_hdim);
    _diff        .resize(_hdim);
    _normDiff    .resize(_nPartitions);
    _psi         .resize(_nPartitions);
    _lambdaMaster.resize(_ldim);
    _lambdaSlave .resize(_ldim);
    _dPhi        .resize(_ldim*_ldim);
    _dPhiOverlap .resize(_nPartitions*_ldim*_ldim);
    _xiOverlap   .resize(_nPartitions*_ldim);
    _fileXiOverlap  .open("xi_overlap.txt");
    _filePsiOverlap .open("psi_overlap.txt");
    _fileDPhiOverlap.open("dphi_overlap.txt");
    _fileDiffOverlap.open("diff_overlap.txt");
    _fileDihedrals  .open("dihedrals.txt");
}

MPSandRest::~MPSandRest()
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


int MPSandRest::TransferLambda()
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


void MPSandRest::CalculateForces( int stepNum, double phiAngle, double psiAngle, double *masses, double *positions, double *forces )
{
//    cout<<"%%%%["<<stepNum<<"]%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
    /* Change the masses and positions array to the STL vectors */
//    vector<double> vMasses   ( masses   , masses   +_nAtoms);
    vector<double> vPositions( positions, positions+_hdim  );
    _alignment->ComputeDerivatives(vPositions);
    vector<double> aPositions( _alignment->GetAlignedX() );

    const int closestMarker = _closestMarkerSearch->FindNearestNeighbors(aPositions).at(0);
//    cout<<"Closest marker = "<<closestMarker<<endl;

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


    fill(_psi.begin()        , _psi.end()        , 0.0);
    fill(_normDiff.begin()   , _normDiff.end()   , 0.0);
    fill(_xiOverlap.begin()  , _xiOverlap.end()  , 0.0);
    fill(_dPhiOverlap.begin(), _dPhiOverlap.end(), 0.0);
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
        int pid = -1;
        try
        {
            pid = partInfo.at(0).first;
//            cout<<"******** Master part [ "<<pid<<" ] = "<<_psi.at(pid)<<"  ********"<<endl;
            _vsandcv.at(pid)->CalculateCV(vPositions, stepNum);


            _xMaster      = _vsandcv.at(pid)->GetX();
            _dxMaster     = _vsandcv.at(pid)->GetDX();

            copy(_vsandcv.at(pid)->GetValue().begin(), _vsandcv.at(pid)->GetValue().end(), _xiOverlap.begin()+pid*_ldim);

            nFailedParts--;
            for ( int i = 1; i < nParts; ++i )
            {
                pid = partInfo.at(i).first;
//                cout<<"******** Slave part [ "<<pid<<" ] = "<<_psi.at(pid)<<"  ********"<<endl;
                try
                {
                    _vsandcv.at(pid)->CalculateCV(vPositions, stepNum);

                    /* Transform xiForce evaluated at master part and pass it to ABF */
                    _xSlave  = _vsandcv.at(pid)->GetX();
                    _dxSlave = _vsandcv.at(pid)->GetDX();
                    if (TransferLambda())
                        throw runtime_error("Error in tranforming CVs!");

                    /* copy the xi value to be save in overlap file */
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
                    if (_psi.at(pid) > 0.05)
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
    }
    catch (const runtime_error &e)
    {
        cout<<"ERROR in MPSandCV >> "<<e.what()<<endl;
    }


    /* Apply the harmonic forces */
    transform( _xCenter.begin(), _xCenter.end(), aPositions.begin(), _diff.begin(), std::minus<double>() );
    double rd = sqrt(inner_product(_diff.begin(), _diff.end(), _diff.begin(), 0.0));
    double rf = 0.0;

    if (rd > _rBias)
    {
        for ( auto &item : _diff )
            item *= _kBias*(rd-_rBias);

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 1, _hdim, _hdim, 1.0,
                    _diff.data(), _hdim, _alignment->GetdAdX().data(), _hdim, 0.0, forces, _hdim);

        rf = sqrt(inner_product(forces, forces+_hdim, forces, 0.0));
    }
    else
        fill(forces, forces+_hdim, 0.0);

    /* Write the SandCV projection information */
    if (stepNum%_logFrequency == 0)
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

        _fileDihedrals<<phiAngle<<"\t"<<psiAngle<<"\t"<<rd<<"\t"<<rf<<endl;
    }

}
