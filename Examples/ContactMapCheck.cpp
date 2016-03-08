#include <iostream>
#include <limits>
#include <typeinfo>
#include <exception>
#include <string>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <stdexcept>
#include <cstring>
#include "PointSetIO.h"
#include "ContactMap.h"
using namespace std;
using namespace chrono;

void Usage ( char *str );
void ParseCommandLine( int argc, char *argv[], int &nTests, double &d0, double &delta );

int main(int argc, char* argv[])
{
    double min = numeric_limits<double>::min();
    double max = numeric_limits<double>::max();
    double eps = numeric_limits<double>::epsilon();
    cout<<"min:"<<min<<", max:"<<max<<" and eps = "<<eps<<endl;

    string path      = "../../data/";
    string inputXs   = path+"xSamples_ala_raw.txt";
    PointSetIO IO;

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


    //--------------------------------------------------------------------------
    // [2] Alignment
    //--------------------------------------------------------------------------
    //__________________________________________________________________________
    // [2-0] Initialize Contact Map
    //
    int    sdim   = 3;
    int    nAtoms = hdim / sdim;
    int    cdim   = nAtoms*(nAtoms-1)/2;

    int nTests;
    double d0;
    double delta;
    ParseCommandLine(argc, argv, nTests, d0, delta);
    ContactMap CMap(sdim, nAtoms, d0);

    //__________________________________________________________________________
    // [2-1] Check the Jacobian of Alignemnt with finite differences
    //


    int iD, jD;
    if(nTests == 0)
        nTests = nSamples;
    cout<<"nTests = "<<nTests<<endl;

    if(delta == 0.0)
        delta  = pow(eps,1.0/3.0);

    time_point<system_clock> start = system_clock::now();
    vector<double> xOrig(hdim, 0.0);
    vector<double> xTemp(hdim, 0.0);
    vector<double> xiMinus, xiPlus;
    vector<double> jacCMap  (nTests*cdim*hdim, 0.0);
    vector<double> jacFinite(nTests*cdim*hdim, 0.0);

    for ( int iN = 0; iN < nTests; ++iN)
    {
        copy( xSamples.begin()+iN*hdim, xSamples.begin()+(iN+1)*hdim, xOrig.begin() );
        CMap.ComputeDerivatives( xOrig );
        copy( CMap.GetdAdX().begin(), CMap.GetdAdX().end(), jacCMap.begin()+iN*cdim*hdim );

        for( iD = 0; iD < hdim; ++iD )
        {
            copy( xOrig.begin(), xOrig.end(), xTemp.begin() );
            xTemp.at(iD) += delta;
            xiPlus = CMap.ComputeAlignment( xTemp );

            copy( xOrig.begin(), xOrig.end(), xTemp.begin() );
            xTemp.at(iD) -= delta;
            xiMinus = CMap.ComputeAlignment( xTemp );

            for ( jD = 0; jD < cdim; ++jD )
                jacFinite.at(iN*cdim*hdim+jD*hdim+iD) = ( xiPlus.at(jD) - xiMinus.at(jD) ) / (2*delta);
        }
    }

    cout<<"["<<nTests<<"] tests has been done."<<endl;
    time_point<system_clock> end = system_clock::now();
    duration<double> elapsed = end - start;
    cout<<"Elapsed time............................."<<elapsed.count()<<"seconds"<<endl;

    //__________________________________________________________________________
    // [2-2] Save the jacobian matrices
    //
    IO.WritePointSet(path+"cmap_dAdx_anlytic.txt", jacCMap  , nTests, cdim*hdim);
    IO.WritePointSet(path+"cmap_dAdx_finite.txt" , jacFinite, nTests, cdim*hdim);
    cout<<"All matrices have been written to files."<<endl;

    return 0;
}



void Usage ( char *str )
{
    cout<<"ERROR:: "<< str <<" input parameter is NULL"<<endl;
}

void ParseCommandLine( int argc, char *argv[], int &nTests, double &d0, double &delta )
{
    bool   ok;

    // The command line arguments are read
    while ( argc > 1 )
    {
        ok = false;
        /* Options: */
        if ( ok == false && !strcmp(argv[1], "-n") )
        {
            argc--;
            argv++;
            if ( argv[1] == NULL )
                Usage ( "n" );
            nTests =  atoi ( argv[1] );
            argc--;
            argv++;
            ok = true;
        }
        if ( ok == false && !strcmp(argv[1], "-d0") )
        {
            argc--;
            argv++;
            if ( argv[1] == NULL )
                Usage ( "d0" );
            d0 =  atof ( argv[1] );
            argc--;
            argv++;
            ok = true;
        }
        if ( ok == false && !strcmp(argv[1], "-dx") )
        {
            argc--;
            argv++;
            if ( argv[1] == NULL )
                Usage ( "dx" );
            delta=  atof ( argv[1] );
            argc--;
            argv++;
            ok = true;
        }
        if ( ok == false )
        {
            cout<<"ERROR:: Can not parse argument "<<argv[1]<<endl;
        }
    }
}
