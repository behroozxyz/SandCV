/*=========================================================================
 * File         : PuckerABF.cpp
 * Module       : Swig
 * Copyright    : (C)opyright 2013
 *                See COPYRIGHT statement in top level directory.
 * Authors      : Behrooz Hashemian
 * Created      : Thu Dec 19, 2013  07:44PM
 * Last modified: Thu Dec 19, 2013  07:44PM
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
#include "PuckerABF.h"

PuckerABF::PuckerABF(int bApplyABF, int nAtoms, int sdim, double timeStep,
                     int logFrequency, int histFrequency, char *coreName,
                     double corralHeight, char *corralFilename,
                     int bLoadHistogram, double *xiMin, double *xiMax, int *nBins, int fullSample, int rampSample)
{
    int ldim = 2;
    _nAtoms = nAtoms;
    _hdim   = sdim*nAtoms;
    vector<int>    vNumBins(nBins,nBins+ldim);
    vector<double> vXiMin(xiMin,xiMin+ldim);
    vector<double> vXiMax(xiMax,xiMax+ldim);
    //______________________________________________________________________________
    // Initialize Puckering
    _puckercv = make_shared<Puckering>();

    //______________________________________________________________________________
    // Initialize ABF
    _abf = make_shared<ABF>( bApplyABF, bLoadHistogram, ldim, sdim, nAtoms, timeStep, logFrequency, histFrequency,
                            corralHeight, corralFilename, coreName, vXiMin, vXiMax, vNumBins, fullSample, rampSample);
    //    cout<<"ABF is initialized."<<endl;
}

PuckerABF::PuckerABF(int bApplyABF, int nAtoms, int sdim, double timeStep,
                     int logFrequency, int histFrequency, char *coreName,
                     double kWall, double *lowerWall, double *upperWall,
                     int bLoadHistogram, double *xiMin, double *xiMax, int *nBins, int fullSample, int rampSample)
{
    int ldim = 2;
    _nAtoms  = nAtoms;
    _hdim    = nAtoms*sdim;
    vector<int>    vNumBins(nBins,nBins+ldim);
    vector<double> vXiMin(xiMin,xiMin+ldim);
    vector<double> vXiMax(xiMax,xiMax+ldim);
    vector<double> vLowerWall(lowerWall,lowerWall+ldim);
    vector<double> vUpperWall(upperWall,upperWall+ldim);
    //______________________________________________________________________________
    // Initialize Puckering
    _puckercv = make_shared<Puckering>();

    //______________________________________________________________________________
    // Initialize ABF
    _abf = make_shared<ABF>( bApplyABF, bLoadHistogram, ldim, sdim, nAtoms, timeStep,
                            logFrequency, histFrequency,
                            kWall, vLowerWall, vUpperWall,
                            coreName, vXiMin, vXiMax, vNumBins, fullSample, rampSample);
    cout<<"ABF is initialized."<<endl;

}



void PuckerABF::CalculateForces( int stepNum, double *masses, double *positions, double *forces )
{
    /* Change the masses and positions array to the STL vectors */
    vector<double> vMasses   (masses   , masses   +_nAtoms);
    vector<double> vPositions(positions, positions+_hdim  );

    /* Compute the CV value of the current configuration and its related Jacobian */
    try
    {
        _puckercv->CalculateCV(vPositions);

        /* Calculate the Abf forces and return it */
        vector<double> forceABF( _abf->CalculateForce( stepNum, vMasses, vPositions, _puckercv->GetValue(), _puckercv->GetJacobian() ) );

        /* Copy the calculated forces to pass to the MD engine */
        copy( forceABF.begin(), forceABF.end(), forces );
    }
    catch (const runtime_error &e)
    {
        fill( forces, forces+_hdim, 0.0 );
        cout<<"NO FORCE is applied because "<<e.what()<<endl;
    }

}

