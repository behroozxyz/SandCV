#include <algorithm>
#include <iostream>
#include <limits>
#include <typeinfo>
#include <exception>
#include <string>
#include <chrono>
#include <iomanip>
#include <stdexcept>
#include "PointSetIO.h"
#include "Corral.h"
using namespace std;
using namespace chrono;


int main(int argc, char* argv[])
{
    double min = numeric_limits<double>::min();
    double max = numeric_limits<double>::max();
    double eps = numeric_limits<double>::epsilon();
    cout<<"min:"<<min<<", max:"<<max<<" and eps = "<<eps<<endl;


    string prefix = "../../data/";
    ofstream fileGrad(prefix+"corralGrad.txt");
    Corral corral(prefix+"ala_corral.txt");
    int dim = 2;
    vector<double> corralGradient = {0.0,0.0};
    vector<double> xiNew          = {0.0,0.0};
    cout<<"nMarkers = "<<corral.GetMarkers().size()/dim<<endl;
    int counter = 0;
//    vector<double> markers(corral.GetMarkers());
//    for ( vector<double>::iterator iterMarker = markers.begin();
//          iterMarker != markers.end(); iterMarker += dim )
//    {
//        copy(iterMarker, iterMarker+dim, xiNew.begin());
    double start  = -10.0;
    double finish =  10.0;
    double step   =  0.1;
    int nSteps = (finish-start)/step;
    int i,j;
    for ( i = 0; i < nSteps; ++i)
        for ( j = 0; j < nSteps; ++j)
        {
            xiNew.at(0) = start+i*step;
            xiNew.at(1) = start+j*step;
            cout<<"**********[ "<<counter<<" ]**********"<<endl;
            corral.CalculateCorral(xiNew, corralGradient);
            fileGrad<<corralGradient.at(0)<<"\t"<<corralGradient.at(1)<<endl;
            ++counter;
        }
    cout<<"Counts = "<<counter<<endl;
    return 0;
}
