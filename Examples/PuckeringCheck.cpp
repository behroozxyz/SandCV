#include <vector>
#include <iostream>
#include <limits>
#include <typeinfo>
#include <exception>
#include <string>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <stdexcept>
#include "Puckering.h"
#include "PointSetIO.h"
using namespace std;
using namespace chrono;

int main(int argc, char* argv[])
{
    time_point<system_clock> start, end;
    duration<double> elapsed;
    double min = numeric_limits<double>::min();
    double max = numeric_limits<double>::max();
    double eps = numeric_limits<double>::epsilon();
    cout<<"min:"<<min<<", max:"<<max<<" and eps = "<<eps<<endl;
    string     prefix = "../../data/";
    string     inputX = prefix+"bgc_test.txt";
    PointSetIO IO;
    Puckering  p;

    int ldim   = 2;
    int sdim   = 3;
    int nAtoms = 6;
    int jD, iN, iD;

    cout<<"Reading configurations: ["<<inputX<<"] ..."<<endl;
    int nPoints = 0;
    int hdim    = 0;
    vector<double> X;
    IO.ReadPointSet( inputX, X, nPoints, hdim );

    if (hdim != nAtoms*sdim)
    {
        cout<<"hdim != nAtoms*sdim";
        return 1;
    }

    double deltaX = sqrt(numeric_limits<double>::epsilon());
    cout<<"detaX = "<<deltaX<<endl;
    vector<double> Xorg   (hdim);
    vector<double> Xtemp  (hdim);
    vector<double> cvplus (2);
    vector<double> cvminus(2);
    vector<double> CV     (nPoints*ldim);
    vector<double> dCV    (nPoints*ldim*hdim);
    vector<double> dCVn   (nPoints*ldim*hdim);

    start = system_clock::now();

    for( int iC = 0; iC < nPoints; ++iC )
    {
        copy(X.begin()+iC*hdim, X.begin()+(iC+1)*hdim, Xorg.begin());
        p.CalculateCV(Xorg);
        copy(p.GetValue().begin()   , p.GetValue().end()   , CV.begin() +iC*ldim);
        copy(p.GetJacobian().begin(), p.GetJacobian().end(), dCV.begin()+iC*ldim*hdim);

        for( iN = 0; iN < nAtoms; ++iN )
            for( iD = 0; iD < sdim; ++iD)
            {
                copy(Xorg.begin(), Xorg.end(), Xtemp.begin());
                Xtemp.at(iN*sdim+iD) += deltaX;
                cvplus = p.CalculateCV(Xtemp);

                copy(Xorg.begin(), Xorg.end(), Xtemp.begin());
                Xtemp.at(iN*sdim+iD) -= deltaX;
                cvminus = p.CalculateCV(Xtemp);

                for ( jD = 0; jD < ldim; ++jD)
                    dCVn.at(iC*ldim*hdim + jD*hdim + iN*sdim + iD ) = 0.5 * (cvplus.at(jD)-cvminus.at(jD)) / deltaX;
            }
    }
    IO.WritePointSet(prefix+"Pucker.txt"    , CV  , nPoints, ldim     );
    IO.WritePointSet(prefix+"dPucker.txt"   , dCV , nPoints, ldim*hdim);
    IO.WritePointSet(prefix+"dPuckerNum.txt", dCVn, nPoints, ldim*hdim);

    end = system_clock::now();
    elapsed = end - start;
    cout<<"Elapsed time for blas........."<<elapsed.count()<<"seconds"<<endl;

    return 0;
}

