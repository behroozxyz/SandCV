#include <iostream>
#include <limits>
#include <typeinfo>
#include <exception>
#include <string>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <stdexcept>
#include "PointSetIO.h"
#include "ProcrustesSuperimposition.h"
using namespace std;
using namespace chrono;

int main(int argc, char* argv[])
{
    double min = numeric_limits<double>::min();
    double max = numeric_limits<double>::max();
    double eps = numeric_limits<double>::epsilon();
    cout<<"min:"<<min<<", max:"<<max<<" and eps = "<<eps<<endl;

    int nLines = 0;
    int nCols  = 0;
    string path      = "../../data/";
    string inputXs   = path+"xSamples_ala_raw.txt";
    string inputXref = path+"xRef_ala.txt";
    PointSetIO IO;

    //__________________________________________________________________________
    // [1-1] Load reference configuraton for alignement
    //
    cout<<"Reading reference configuration: ["<<inputXref<<"] ..."<<endl;
    vector<double> xRef;
    IO.ReadPointSet( inputXref, xRef, nLines, nCols );
    cout<<"The reference configuration is loaded."<<endl;
    cout<<"\t dimension = "<<nCols<<endl;
    cout<<"\t frames    = "<<nLines<<endl;
    if(nLines != 1 )
        cerr<<__func__<<" >> More than one reference configuration is proveded!"<<endl;

    // Set active list
    int sdim = 3;
    int nAtoms = nCols/sdim;
    vector<int> activeList(nAtoms);
    iota(activeList.begin(),activeList.end(),0);
//    activeList.erase(activeList.begin()+8);
//    activeList.erase(activeList.begin()+6);
//    activeList.erase(activeList.begin()+4);
//    activeList.erase(activeList.begin()+2);
    cout<<"size activeList = "<<activeList.size()<<endl;

    //__________________________________________________________________________
    // [1-2] Load high-dimensional sample points
    //
    cout<<"Reading high-dimensional sample points: ["<<inputXs<<"] ..."<<endl;

    int hdim     = 0;
    int nSamples = 0;
    vector<double> xSamples;

    IO.ReadPointSet( inputXs, xSamples, nSamples, hdim );
    cout<<"High-dimensional sample points are loaded."<<endl;
    cout<<"\t dimension = "<<hdim<<endl;
    cout<<"\t frames    = "<<nSamples<<endl;
    if(nCols != hdim)
        cerr<<__func__<<" >> Dimension of sample points does not match with reference configuration!"<<endl;



    //--------------------------------------------------------------------------
    // [2] Alignment
    //--------------------------------------------------------------------------
    //__________________________________________________________________________
    // [2-0] Initialize Procrustes analysis
    //
    for (auto &iter : xRef)
        iter += 1.0e2;
    ProcrustesSuperimposition PS( sdim, xRef, activeList );

    //__________________________________________________________________________
    // [2-1] Check the Jacobian of Alignemntwith finite differences
    //
    int iD, jD;
    int nTests;
    if(argv[1] == NULL)
        nTests = nSamples;
    else
        nTests = atoi(argv[1]);

    double delta  = pow(eps,1.0/3.0);
    time_point<system_clock> start = system_clock::now();
    vector<double> xOrig(hdim, 0.0);
    vector<double> xTemp(hdim, 0.0);
    vector<double> xiMinus, xiPlus;
    vector<double> jacProc  (nTests*hdim*hdim, 0.0);
    vector<double> jacFinite(nTests*hdim*hdim, 0.0);

    for ( int iN = 0; iN < nTests; ++iN)
    {
        copy( xSamples.begin()+iN*hdim, xSamples.begin()+(iN+1)*hdim, xOrig.begin() );
        PS.ComputeDerivatives( xOrig );
        copy( PS.GetdAdX().begin(), PS.GetdAdX().end(), jacProc.begin()+iN*hdim*hdim );

        for( iD = 0; iD < hdim; ++iD )
        {
            copy( xOrig.begin(), xOrig.end(), xTemp.begin() );
            xTemp.at(iD) += delta;
            xiPlus = PS.ComputeAlignment( xTemp );

            copy( xOrig.begin(), xOrig.end(), xTemp.begin() );
            xTemp.at(iD) -= delta;
            xiMinus = PS.ComputeAlignment( xTemp );

            for ( jD = 0; jD < hdim; ++jD )
                jacFinite.at(iN*hdim*hdim+jD*hdim+iD) = ( xiPlus.at(jD) - xiMinus.at(jD) ) / (2*delta);
        }
    }

    cout<<"["<<nTests<<"] tests has been done."<<endl;
    time_point<system_clock> end = system_clock::now();
    duration<double> elapsed = end - start;
    cout<<"Elapsed time............................."<<elapsed.count()<<"seconds"<<endl;

    //__________________________________________________________________________
    // [2-2] Save the jacobian matrices
    //
    IO.WritePointSet(path+"dAdx_proc.txt"  , jacProc  , nTests, hdim*hdim);
    IO.WritePointSet(path+"dAdx_finite.txt", jacFinite, nTests, hdim*hdim);
    cout<<"All matrices have been written to files."<<endl;

    return 0;
}
