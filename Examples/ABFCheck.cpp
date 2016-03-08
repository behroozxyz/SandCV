#include <algorithm>
#include <iostream>
#include <limits>
#include <typeinfo>
#include <exception>
#include <string>
#include <chrono>
#include <iomanip>
#include <stdexcept>
#include <string>
#include "PointSetIO.h"
#include "SandABF.h"
using namespace std;
using namespace chrono;


int main(int argc, char* argv[])
{
    double min = numeric_limits<double>::min();
    double max = numeric_limits<double>::max();
    double eps = numeric_limits<double>::epsilon();
    cout<<"min:"<<min<<", max:"<<max<<" and eps = "<<eps<<endl;

    const char* xNodesFilename  = "../../../works/molsim/alanine.md/abf.sandcv/common/mappings/x_ala_abf_vac1_noh_all_connected_step4Rot_dihedCut_knn12_nb10.txt";
    const char* xiNodesFilename = "../../../works/molsim/alanine.md/abf.sandcv/common/mappings/xi_ala_abf_vac1_noh_all_connected_step4Rot_dihedCut_knn12_nb10.txt";
    const char* xiSeedsFilename = "../../../works/molsim/alanine.md/abf.sandcv/common/mappings/xis_ala_abf_vac1_noh_all_connected_step4Rot_dihedCut_knn12.txt";
    const char* corralFilename  = "../../../works/molsim/alanine.md/abf.sandcv/common/corral/markers_ala_abf_vac1_noh_all_connected_step4Rot_dihedCut_knn12_w1.5_nb40.txt";
    const char* xRefFilename    = "../../../works/molsim/alanine.md/abf.sandcv/common/x_ref_ala_vac.txt";
    string xNewFilename   = "../../data/xNew_ala_abf_vac_120K.txt";

    //--------------------------------------------------------------------------
    // [1] Configuration parameters: equivalent to the TCL configuration file
    //--------------------------------------------------------------------------
    int     bApplyABF      = 1;
    int     nAtoms         = 10;
    int     sdim           = 3;
    int     logFrequency   = 10;
    int     histFrequency  = 1000;
    int     nInvolvedAtoms = 10;
    int     ldim           = 2;
    int     hdim           = sdim*nAtoms;
    int     bDistance      = 0;
    int     bLoadHistogram = 0;
    double  timeStep       = 1.0;
    double  corralHeight   = 100.0;
    const char *coreName   = "checkABF";
    vector<int> nBins      = {72,72};
    vector<int> involvedAtoms(nAtoms);
    iota(involvedAtoms.begin(),involvedAtoms.end(),0);
    vector<double> xiMin   = {-10.0,-12.0};
    vector<double> xiMax   = { 12.0, 10.0};
    vector<double> masses  = {12.01099967956543, 15.99899959564209,
                              14.00699996948242, 12.01099967956543,
                              12.01099967956543, 12.01099967956543,
                              15.99899959564209, 14.00699996948242,
                              12.01099967956543, 12.01099967956543};

    //--------------------------------------------------------------------------
    // [2] Initialize SandCV and ABF
    //--------------------------------------------------------------------------
    SandABF sandabf( bApplyABF, nAtoms, sdim, timeStep, logFrequency, histFrequency,
                     (char *)coreName, (char *)xRefFilename, nInvolvedAtoms, involvedAtoms.data(),
                     ldim, bDistance,
                     (char *)xiNodesFilename, (char *)xNodesFilename, (char *)xiSeedsFilename,
                     corralHeight, (char *)corralFilename, bLoadHistogram,
                     xiMin.data(), xiMax.data(), nBins.data());

    cout<<"SandABF is initialized"<<endl;


    //--------------------------------------------------------------------------
    // [3] Read new configuration points
    //--------------------------------------------------------------------------
    cout<<"Reading high-dimensional new points: ["<<xNewFilename<<"] ..."<<endl;

    int nNewPoints = 0;
    int nCols      = 0;
    vector<double> xNew;
    PointSetIO IO;
    IO.ReadPointSet( xNewFilename, xNew, nNewPoints, nCols );

    cout<<"High-dimensional new points are loaded."<<endl;
    cout<<"\t dimension = "<<nCols<<endl;
    cout<<"\t frames    = "<<nNewPoints<<endl;
    if(nCols != hdim)
        cerr<<__func__<<" >> Dimension of new points does not match with high-dimensional node points!"<<endl;


    //--------------------------------------------------------------------------
    // [4] Run the simulation!
    //--------------------------------------------------------------------------
    time_point<system_clock> start = system_clock::now();
    int nTests = nNewPoints;
    vector<double> force(hdim);
    cout<<"Simulation of "<<nTests<<" points is started..."<<endl;
    for ( int iN = 0; iN < nTests; ++iN)
    {
        sandabf.CalculateForces( iN, masses.data(), xNew.data()+iN*hdim, force.data() );
    }


    time_point<system_clock> end = system_clock::now();
    duration<double> elapsed = end - start;
    cout<<"Elapsed time............................."<<elapsed.count()<<"seconds"<<endl;

    return 0;
}
