/*=========================================================================
 * File         : SandCMapABF.cpp
 * Module       : Swig
 * Copyright    : (C)opyright 2013
 *                See COPYRIGHT statement in top level directory.
 * Authors      : Behrooz Hashemian
 * Created      : Tue Aug 13, 2013  11:49AM
 * Last modified: Tue Oct 29, 2013  12:52PM
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
#include "SandCMapABF.h"

SandCMapABF::SandCMapABF( int bApplyABF, int nAtoms, int sdim, double timeStep,
                  int logFrequency, int histFrequency, char *coreName,
                  char *xRefFilename, int nInvolvedAtoms, int *involvedAtoms,
                  int ldim, int bDistance,
                  char *xiNodesFilename, char *xNodesFilename,
                  char *xiSeedsFilename,
                  double corralHeight, char *corralFilename,
                  int bLoadHistogram, double *xiMin, double *xiMax, int *nBins)
    : _bDistance(bDistance)
{
    int adim = _bDistance ? ldim+1 : ldim;
    vector<int>    vNumBins(nBins,nBins+adim);
    vector<double> vXiMin(xiMin,xiMin+adim);
    vector<double> vXiMax(xiMax,xiMax+adim);
    //______________________________________________________________________________
    // Initialize SandCV
    this->UpdateSand( nAtoms, sdim, xRefFilename, nInvolvedAtoms, involvedAtoms,
                      ldim, xiNodesFilename, xNodesFilename, xiSeedsFilename,
                      vXiMin, vXiMax );
    //______________________________________________________________________________
    // Initialize ABF
    _abf = make_shared<ABF>( bApplyABF, bLoadHistogram, adim, sdim, nAtoms, timeStep, logFrequency, histFrequency,
                            corralHeight, corralFilename, coreName, vXiMin, vXiMax, vNumBins );
//    cout<<"ABF is initialized."<<endl;
}



SandCMapABF::SandCMapABF( int bApplyABF, int nAtoms, int sdim, double timeStep,
                  int logFrequency, int histFrequency, char *coreName,
                  char *xRefFilename, int nInvolvedAtoms, int *involvedAtoms,
                  int ldim, int bDistance,
                  char *xiNodesFilename, char *xNodesFilename,
                  char *xiSeedsFilename,
                  double kWall, double *lowerWall, double *upperWall,
                  int bLoadHistogram, double *xiMin, double *xiMax, int *nBins)
    : _bDistance(bDistance)
{
    int adim = _bDistance ? ldim+1 : ldim;
    vector<int>    vNumBins(nBins,nBins+adim);
    vector<double> vXiMin(xiMin,xiMin+adim);
    vector<double> vXiMax(xiMax,xiMax+adim);
    vector<double> vLowerWall(lowerWall,lowerWall+adim);
    vector<double> vUpperWall(upperWall,upperWall+adim);
    //______________________________________________________________________________
    // Initialize SandCV
    this->UpdateSand( nAtoms, sdim, xRefFilename, nInvolvedAtoms, involvedAtoms,
                      ldim, xiNodesFilename, xNodesFilename, xiSeedsFilename,
                      vXiMin, vXiMax);

    //______________________________________________________________________________
    // Initialize ABF
    _abf = make_shared<ABF>( bApplyABF, bLoadHistogram, adim, sdim, nAtoms, timeStep,
                            logFrequency, histFrequency,
                            kWall, vLowerWall, vUpperWall,
                            coreName, vXiMin, vXiMax, vNumBins );
    cout<<"ABF is initialized."<<endl;

}



void SandCMapABF::UpdateSand(int nAtoms, int sdim,
                         char *xRefFilename, int nInvolvedAtoms, int *involvedAtoms,
                         int ldim, char *xiNodesFilename, char *xNodesFilename,
                         char *xiSeedsFilename,
                         const vector<double> &vXiMin,
                         const vector<double> &vXiMax)
{
    int nLines1 = 0;
    int nLines2 = 0;
    int nCols1  = 0;
    int nCols2  = 0;

    _nAtoms = nAtoms;
    _hdim   = sdim*nAtoms;
    vector<int>    vInvolvedAtoms(involvedAtoms,involvedAtoms+nInvolvedAtoms);
    vector<double> xRef;
    vector<double> xiNodes;
    vector<double> xNodes;
    vector<double> xiSeeds;
    PointSetIO     reader;

    //______________________________________________________________________________
    // Read xRefFilename
    reader.ReadPointSet(xRefFilename, xRef, nLines1, nCols1);
    if (nLines1 < 1)
        cerr<<"xRefFilename("<<xRefFilename<<") is empty."<<endl;
    else if (nLines1 > 1)
        cerr<<"xRefFilename("<<xRefFilename<<") includes more than one configuration."<<endl;

    //______________________________________________________________________________
    // Read xiNodesFilename
    reader.ReadPointSet(xiNodesFilename , xiNodes , nLines1, nCols1);
    if (nLines1 < 1)
        cerr<<"xiNodesFilename("<<xiNodesFilename<<") is empty."<<endl;
    if (nCols1 != ldim)
        cerr<<"Dimension of nodes points does not match with ldim ("<<nCols1<<"!="<<ldim<<")."<<endl;
    double spacing = xiNodes.back();
    xiNodes.erase(xiNodes.end()-ldim,xiNodes.end());
    --nLines1;
    cout<<"spacing = "<<spacing<<endl;

    //______________________________________________________________________________
    // Read xNodesFilename
    reader.ReadPointSet(xNodesFilename, xNodes  , nLines2, nCols2);
    if (nLines2 < 1)
        cerr<<"xNodesFilename("<<xNodesFilename<<") is empty."<<endl;
    if (nLines2 != nLines1)
        cerr<<"Number of low and high dimensional points does not match("
           <<nLines1<<"!="<<nLines2<<")."<<endl;

    //______________________________________________________________________________
    // Read xiSeedsFilename
    reader.ReadPointSet(xiSeedsFilename , xiSeeds , nLines2, nCols2);
    if (nLines2 < 1)
        cerr<<"xiSeedsFilename("<<xiSeedsFilename<<") is empty."<<endl;
    if (nCols2 != nCols1)
        cerr<<"Dimension of seeds and nodes points does not match("
           <<nCols2<<"!="<<nCols1<<")."<<endl;

    /* Remove the seed points outside of the defined xi boundaries */
    bool isOut = false;
    int  iD    = 0;
    int  iN    = 0;
    while ( iN < nLines2 )
    {
        isOut = false;
        for ( iD = 0; iD < ldim; ++iD )
            isOut = isOut || xiSeeds.at(ldim*iN+iD) < vXiMin.at(iD) || xiSeeds.at(ldim*iN+iD) > vXiMax.at(iD) ;

        if (isOut)
        {
//            cout<<"removed seed = ";
//            for ( iD = 0; iD < _ldim; ++iD )
//                cout<<xiSeeds.at(_ldim*iN+iD)<<" ";
//            cout<<endl;
            xiSeeds.erase( xiSeeds.begin()+ldim*iN, xiSeeds.begin()+ldim*(iN+1) );
            --nLines2;
        }
        else
            ++iN;
    }
    reader.WritePointSet("xi_seeds.txt", xiSeeds, nLines2, ldim);

    //______________________________________________________________________________
    // Initialize SandCV
    _sandcv = make_shared<SandCV>( sdim, ldim, _hdim, spacing, xRef, vInvolvedAtoms, xiNodes, xNodes, xiSeeds);
    cout<<"SandCV is initialized."<<endl;

}



void SandCMapABF::CalculateForces( int stepNum, double *masses, double *positions, double *forces )
{
    /* Change the masses and positions array to the STL vectors */
    vector<double> vMasses   (masses   , masses   +_nAtoms);
    vector<double> vPositions(positions, positions+_hdim  );

    /* Compute the CV value of the current configuration and its related Jacobian */
    try
    {
        _bDistance ? _sandcv->CalculateCVdist(vPositions) : _sandcv->CalculateCV(vPositions);

        /* Calculate the Abf forces and return it */
        vector<double> forceABF( _abf->CalculateForce( stepNum, vMasses, vPositions, _sandcv->GetValue(), _sandcv->GetJacobian() ) );

        /* Copy the calculated forces to pass to the MD engine */
        copy( forceABF.begin(), forceABF.end(), forces );
    }
    catch (const runtime_error &e)
    {
        fill( forces, forces+_hdim, 0.0 );
        cout<<"NO FORCE is applied because "<<e.what()<<endl;
    }

}

