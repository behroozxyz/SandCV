/*=========================================================================
 * File         : SandCMapSMD.cpp
 * Module       : Swig
 * Copyright    : (C)opyright 2013
 *                See COPYRIGHT statement in top level directory.
 * Authors      : Behrooz Hashemian
 * Created      : Tue Aug 13, 2013  11:49AM
 * Last modified: Thu Aug 22, 2013  03:02PM
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
#include "SandCMapSMD.h"
SandCMapSMD::SandCMapSMD(int bApplySMD, int nAtoms, int sdim, double timeStep,
                           int logFrequency, char *coreName,
                           double d0,
                           int ldim, int bDistance,
                           char *xiNodesFilename, char *xNodesFilename, char *xiSeedsFilename,
                           int bWithReturn,
                           double pullingRate, double pullingConstant, double pullingMaxLength,
                           double kWall, double *lowerWall, double *upperWall,
                           double xiMin, double xiMax)
    : _bDistance(bDistance)
{

    int nLines1 = 0;
    int nLines2 = 0;
    int nCols1  = 0;
    int nCols2  = 0;

    int bdim = _bDistance ? ldim+1 : ldim;
    _hdim = sdim*nAtoms;
    vector<double> vLowerWall(lowerWall,lowerWall+bdim);
    vector<double> vUpperWall(upperWall,upperWall+bdim);
    vector<double> xiNodes;
    vector<double> xNodes;
    vector<double> xiSeeds;
    PointSetIO     reader;
    //______________________________________________________________________________
    // Read xiNodesFilename
    try
    {
        reader.ReadPointSet(xiNodesFilename , xiNodes , nLines1, nCols1);
    }
    catch( runtime_error &e)
    {
            cerr<<__func__<<" >> "<<e.what()<<endl;
    }
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
    try
    {
        reader.ReadPointSet(xNodesFilename, xNodes  , nLines2, nCols2);
    }
    catch( runtime_error &e)
    {
            cerr<<__func__<<" >> "<<e.what()<<endl;
    }
    if (nLines2 < 1)
        cerr<<"xNodesFilename("<<xNodesFilename<<") is empty."<<endl;
    if (nLines2 != nLines1)
        cerr<<"Number of low and high dimensional points does not match("
           <<nLines1<<"!="<<nLines2<<")."<<endl;

    //______________________________________________________________________________
    // Read xiSeedsFilename
    try
    {
        reader.ReadPointSet(xiSeedsFilename , xiSeeds , nLines2, nCols2);
    }
    catch( runtime_error &e)
    {
            cerr<<__func__<<" >> "<<e.what()<<endl;
    }
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
        isOut = xiSeeds.at(ldim*iN+iD) < xiMin || xiSeeds.at(ldim*iN+iD) > xiMax ;

        if (isOut)
        {
            xiSeeds.erase( xiSeeds.begin()+ldim*iN, xiSeeds.begin()+ldim*(iN+1) );
            --nLines2;
        }
        else
            ++iN;
    }
    reader.WritePointSet("xi_seeds.txt", xiSeeds, nLines2, ldim);

    //______________________________________________________________________________
    // Initialize SandCVCMap
    sandcv = new SandCVCMap( sdim, ldim, nAtoms, spacing, d0, xiNodes, xNodes, xiSeeds);
    cout<<"SandCVCMap is initialized."<<endl;

    //______________________________________________________________________________
    // Initialize SMD
    smd = new SMD( bApplySMD, _bDistance, bWithReturn, ldim, sdim, nAtoms, timeStep,
                   pullingRate, pullingConstant, pullingMaxLength,
                   logFrequency,  coreName,
                   kWall,  vLowerWall, vUpperWall,
                   xiMin, xiMax);
    cout<<"SMD is initialized."<<endl;
}


SandCMapSMD::~SandCMapSMD()
{
    if (sandcv) delete sandcv;
    if (smd) delete smd;
}


void SandCMapSMD::CalculateForces( double *positions, double *forces )
{
    /* Change the masses and positions array to the STL vectors */
    vector<double> vPositions(positions, positions+_hdim  );

    /* Compute the CV value of the current configuration and its related Jacobian */
    try
    {
        _bDistance ? sandcv->CalculateCVdist(vPositions) : sandcv->CalculateCV(vPositions);

        /* Calculate the Abf forces and return it */
        vector<double> forceSMD( smd->CalculateForce( vPositions, sandcv->GetValue(), sandcv->GetJacobian(), sandcv->GetDistance() ) );

        /* Copy the calculated forces to pass to the MD engine */
        copy( forceSMD.begin(), forceSMD.end(), forces );
    }
    catch (const runtime_error &e)
    {
        fill( forces, forces+_hdim, 0.0 );
        cout<<"NO FORCE is applied because "<<e.what()<<endl;
    }

}
