#include <iostream>
#include <typeinfo>
#include <exception>
#include <string>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <stdexcept>
#include <random>
#include "LocalMaxEntropy.h"
#include "ANNSearch.h"
using namespace std;
using namespace chrono;

int main(int argc, char* argv[])
{
    int    iS, iN, jN, kN, iD, jD; // counters
    int    correct;
    int    knn;
    double one;
    double xError;
    double dxError;
    double hxError;
    double gamma;
    double spacing;

    vector<int>    idxNN;
    vector<double> xNodes;
    vector<double> xSample;
    vector<double> xLme;
    vector<double> dxLme;
    vector<double> hxLme;
    vector<double> pa;
    vector<double> dpa;
    vector<double> hpa;

    // Chronometer
    time_point<system_clock> start, stop;
    duration<double> elapsed;

    //__________________________________________________________________________
    // Check LME shape fucntions in 1D or 2D
    //
    int    dim      = 2;
    int    dim2     = dim*dim;
    double tol0     = 1e-6;
    double xTol     = 1e-13;
    double dxTol    = 1e-12;
    double hxTol    = 1e-11;
    int    nSamples = 1e2;
    int    nNodesX  = 11; // in each direction
    int    nNodes   = pow(nNodesX,dim); // in each direction
    double upper    = 1.0;
    double lower    = 0.0;

    // calculate sapcing
    if (upper < lower) swap(upper,lower);
    spacing = (upper - lower) / double(nNodesX-1);

    // Initialize random number generators
    random_device randomDevice; //random device
    minstd_rand generator(randomDevice()); //linear congruential engine
    uniform_real_distribution<double> uniformDis(lower+spacing, upper-spacing); // uniform distribution

    // Generate node points
    xNodes.resize(nNodes*dim);
    switch (dim)
    {
    case 1 :
        for ( iN = 0; iN < nNodesX; ++iN )
            xNodes.at(iN) = lower + double(iN)*spacing;
        break;
    case 2 :
        for ( iN = 0; iN < nNodesX; ++iN )
            for ( jN = 0; jN < nNodesX; ++jN )
            {
                xNodes.at(iN*nNodesX*dim+jN*dim+0) = lower + double(iN)*spacing;
                xNodes.at(iN*nNodesX*dim+jN*dim+1) = lower + double(jN)*spacing;
            }
        break;
    case 3 :
        for ( iN = 0; iN < nNodesX; ++iN )
            for ( jN = 0; jN < nNodesX; ++jN )
                for ( kN = 0; kN < nNodesX; ++kN )
            {
                xNodes.at(iN*nNodesX*nNodesX*dim+jN*nNodesX*dim+kN*dim+0) = lower + double(iN)*spacing;
                xNodes.at(iN*nNodesX*nNodesX*dim+jN*nNodesX*dim+kN*dim+1) = lower + double(jN)*spacing;
                xNodes.at(iN*nNodesX*nNodesX*dim+jN*nNodesX*dim+kN*dim+2) = lower + double(kN)*spacing;
            }
        break;
    default :
        cerr<<"The given dimension ( "<< dim << " ) is not supported."<<endl;
    }

    gamma = 1.0;

    xSample.resize(dim);
    xLme.resize(dim);
    dxLme.resize(dim2);
    hxLme.resize(dim2*dim);

    // LME object
    LocalMaxEntropy LME( dim, xNodes, spacing, tol0, gamma );


    correct = 0;
    start = system_clock::now(); //start timer
    for ( iS = 0; iS < nSamples; ++iS )
    {
        for ( auto &xs : xSample )
            xs= uniformDis(generator);

        LME.ComputeLme( xSample, 2 );
        pa    = LME.GetShapeFunctions();
        dpa   = LME.GetShapeFunctionsGradient();
        hpa   = LME.GetShapeFunctionsHessian();
        idxNN = LME.GetNearestNeighbors();
        knn   = idxNN.size();

        one = 1.0;
        fill(xLme.begin() , xLme.end() , 0.0);
        fill(dxLme.begin(), dxLme.end(), 0.0);
        fill(hxLme.begin(), hxLme.end(), 0.0);
        for ( iN = 0; iN < knn; ++iN )
        {
//            one -= pa.at(iN);
            for ( iD = 0; iD < dim; ++iD )
                xLme.at(iD) += pa.at(iN) * xNodes.at(idxNN.at(iN)*dim+iD);
            for ( iD = 0; iD < dim; ++iD )
                for ( jD = 0; jD < dim; ++jD )
                    dxLme.at(iD*dim+jD) += dpa.at(iN*dim+jD) * xNodes.at(idxNN.at(iN)*dim+iD);
            for ( iD = 0; iD < dim; ++iD )
                for ( jD = 0; jD < dim2; ++jD )
                    hxLme.at(iD*dim2+jD) += hpa.at(iN*dim2+jD) * xNodes.at(idxNN.at(iN)*dim+iD);
        }

        // Calculate error
        xError = 0.0;
        for ( iD = 0; iD < dim; ++iD )
            xError += (xSample.at(iD)-xLme.at(iD)) * (xSample.at(iD)-xLme.at(iD));
        xError = sqrt(xError);

        dxError = 0.0;
        for ( iD = 0; iD < dim; ++iD )
            dxError += (1-dxLme.at(iD*dim+iD)) * (1-dxLme.at(iD*dim+iD));
        dxError = sqrt(dxError);

        hxError = 0.0;
        for ( iD = 0; iD < dim2*dim; ++iD )
            hxError += hxLme.at(iD) * hxLme.at(iD);
        hxError = sqrt(hxError);


        if (hxError > hxTol)
            cout<<"hxError = "<<hxError<<endl;
        if (dxError > dxTol)
            cout<<"dxError = "<<dxError<<endl;
        if (xError > xTol)
            cout<<"xError  = "<<xError<<endl;
        if (xError <= xTol && dxError <= dxTol && hxError < hxTol)
            ++correct;
    }
    cout<<"LME in "<<dim<<"-dimensional space: "<<endl;
    cout<<correct<<" [out of "<<nSamples<<"] sample points calculated correctly."<<endl;
    stop = system_clock::now(); //stop timer
    elapsed = stop - start;
    cout<<"Elapsed time............................."<<elapsed.count()<<"seconds"<<endl;

    return 0;
}
