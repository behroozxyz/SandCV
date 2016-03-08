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
#include "SandCVCMap.h"
using namespace std;
using namespace chrono;


int main(int argc, char* argv[])
{
    double min = numeric_limits<double>::min();
    double max = numeric_limits<double>::max();
    double eps = numeric_limits<double>::epsilon();
    cout<<"min:"<<min<<", max:"<<max<<" and eps = "<<eps<<endl;

    string prefix    = "../../data/";
    string postfix   = "ala_dihed.txt";
    string inputXn   = prefix+"xNodes_"+postfix;
    string inputXin  = prefix+"xiNodes_"+postfix;
    string inputXis  = prefix+"xiSeeds_"+postfix;
    string inputXs   = prefix+"xSamples_ala_pm160.txt";//+postfix;
    string inputDihed= prefix+"dihedSamples_ala_safe.txt";//+postfix;
    string inputXref = prefix+"xRef_ala.txt";

    PointSetIO IO;

    //--------------------------------------------------------------------------
    // [1] Load point sets
    //--------------------------------------------------------------------------
    //__________________________________________________________________________
    // [1-0] Load high-dimensional node points
    //
    cout<<"Reading high-dimensional node points: ["<<inputXn<<"] ..."<<endl;

    int nFrames = 0;
    int hdim    = 0;
    vector<double> xNodes;

    IO.ReadPointSet( inputXn, xNodes, nFrames, hdim );
    cout<<"\t dimension = "<<hdim<<endl;
    cout<<"\t frames    = "<<nFrames<<endl;

    //__________________________________________________________________________
    // [1-1] Load low-dimensional node points, spacing
    //
    cout<<"Reading low-dimensional node points: ["<<inputXin<<"] ..."<<endl;

    int nLines = 0;
    int ldim   = 0;
    vector<double> xiNodes;

    IO.ReadPointSet( inputXin, xiNodes, nLines, ldim );
    double spacing = xiNodes.back();
    xiNodes.erase(xiNodes.end()-ldim,xiNodes.end());
    cout<<"Low-dimensional nodes are loaded."<<endl;
    cout<<"\t dimension = "<<ldim<<endl;
    cout<<"\t frames    = "<<nLines-1<<endl;
    cout<<"\t spacing   = "<<spacing<<endl;

    if (nFrames != nLines-1)
    {
        cerr<<__func__<<" >> Number of points in low and high-dimensional space does not match!"<<endl;
        throw(invalid_argument("xiNodes file is not valid."));
    }

    //__________________________________________________________________________
    // [1-2] Load low-dimensional seed points
    //
    cout<<"Reading low-dimensional seed points: ["<<inputXis<<"] ..."<<endl;

    int nSeeds = 0;
    int nCols  = 0;
    vector<double> xiSeeds;

    IO.ReadPointSet( inputXis, xiSeeds, nSeeds, nCols );
    cout<<"Low-dimensional seeds are loaded."<<endl;
    cout<<"\t dimension = "<<nCols<<endl;
    cout<<"\t frames    = "<<nSeeds<<endl;

    if(nCols != ldim)
        cerr<<__func__<<" >> Dimension of seed points does not match with low-dimensional node points!"<<endl;

    //__________________________________________________________________________
    // [1-3] Load high-dimensional sample points
    //
    cout<<"Reading high-dimensional sample points: ["<<inputXs<<"] ..."<<endl;

   int nSamples = 0;
    vector<double> xSamples;

    IO.ReadPointSet( inputXs, xSamples, nSamples, nCols );
    cout<<"High-dimensional sample points are loaded."<<endl;
    cout<<"\t dimension = "<<hdim<<endl;
    cout<<"\t frames    = "<<nSamples<<endl;
    if(nCols != hdim)
        cerr<<__func__<<" >> Dimension of sample points does not match with high-dimensional node points!"<<endl;

    //__________________________________________________________________________
    // [1-3-1] Load dihedral angles of the high-dimensional sample points
    //
    cout<<"Reading dihedral angles of sample points: ["<<inputDihed<<"] ..."<<endl;

    vector<double> dihedrals;

    IO.ReadPointSet( inputDihed, dihedrals, nLines, nCols );
    cout<<"Dihedral angles of sample points are loaded."<<endl;
    cout<<"\t dimension = "<<nCols<<endl;
    cout<<"\t frames    = "<<nLines<<endl;
    if(nCols != 2)
        cerr<<__func__<<" >> Dimension of dihedral angles are not equal to two!"<<endl;
    if(nLines != nSamples)
        cerr<<__func__<<" >> Number of dihedrals angles and samples points are not the same!"<<endl;

    //__________________________________________________________________________
    // [1-4] Load reference configuraton for alignement
    //
    cout<<"Reading reference configuration: ["<<inputXref<<"] ..."<<endl;

    vector<double> xRef;

    IO.ReadPointSet( inputXref, xRef, nLines, nCols );
    cout<<"The reference configuration is loaded."<<endl;
    cout<<"\t dimension = "<<nCols<<endl;
    cout<<"\t frames    = "<<nLines<<endl;
    if(nCols != hdim)
        cerr<<__func__<<" >> Dimension of reference configuration does not match with high-dimensional node points!"<<endl;
    if(nLines != 1 )
        cerr<<__func__<<" >> More than one reference configuration is proveded!"<<endl;

    int sdim = 3;
    int nAtoms = hdim/sdim;
    vector<int> activeList(nAtoms);
    iota(activeList.begin(),activeList.end(),0);
//    activeList.erase(activeList.begin()+8);
//    activeList.erase(activeList.begin()+6);
//    activeList.erase(activeList.begin()+4);
//    activeList.erase(activeList.begin()+2);
//    cout<<"size activeList = "<<activeList.size()<<endl;


    //--------------------------------------------------------------------------
    // [2] SandCVCMap
    //--------------------------------------------------------------------------
    //__________________________________________________________________________
    // [2-0] Initialize SandCVCMap
    //
    double d0 = 5.0;
    SandCVCMap SandCVCMap(sdim, ldim, nAtoms, spacing, d0, xiNodes, xNodes, xiSeeds);
//    xSamples = SandCVCMap.GetHighSeeds();
//    cout<<"size xSample ="<<xSamples.size()<<endl;
    //__________________________________________________________________________
    // [2-1] Check the Jacobian of SandCVCMap with finite differences
    //
    int iD, jD;
    int iStart = 0;
    int nTests = 10;//nSamples;
    ++ldim; // to consider distance as a CV
    double deltaX = 1e-3;//pow(eps,1.0/3.0)*1000;
    time_point<system_clock> start = system_clock::now();
    vector<double> xOrig(hdim, 0.0);
    vector<double> xTemp(hdim, 0.0);
    vector<double> xiPlus;
    vector<double> xiMinus;
    vector<double> xiProjected(nTests*ldim, 0.0);
    vector<double> jacSandCVCMap(nTests*hdim*ldim, 0.0);
    vector<double> jacFinite(nTests*hdim*ldim, 0.0);
    for ( int iN = 0; iN < nTests; ++iN)
    {
        copy( xSamples.begin()+(iStart+iN)*hdim, xSamples.begin()+(iStart+iN+1)*hdim, xOrig.begin() );
//        cout<<"xOrig = [";
//        for (auto iter:xOrig)
//            cout<<iter<<" ";
//        cout<<"]"<<endl;
        SandCVCMap.CalculateCVdist( xOrig );
        copy( SandCVCMap.GetJacobian().begin(), SandCVCMap.GetJacobian().end(), jacSandCVCMap.begin()+iN*hdim*ldim  );
        copy( SandCVCMap.GetValue().begin(), SandCVCMap.GetValue().end(), xiProjected.begin()+iN*ldim  );

        //        cout<<"dihedral [ "<<dihedrals.at((iStart+iN)*2)<<" , "<<dihedrals.at((iStart+iN)*2+1)<<"]"<<endl;
//        cout<<"cv_value [ "<<SandCVCMap.GetValue().at(0)<<" , "<<SandCVCMap.GetValue().at(1)<<"]"<<endl;
//        cout<<"________________________________________"<<endl;
//        for( iD = 0; iD < hdim; ++iD )
//        {
//            copy( xOrig.begin(), xOrig.end(), xTemp.begin() );
//            xTemp.at(iD) += deltaX;
//            SandCVCMap.CalculateCVdist( xTemp );
//            xiPlus = SandCVCMap.GetValue();
////            cout<<"+cv = [ "<<xiPlus.at(0)<<" , "<<xiPlus.at(1)<<"]"<<endl;

//            copy( xOrig.begin(), xOrig.end(), xTemp.begin() );
//            xTemp.at(iD) -= deltaX;
//            SandCVCMap.CalculateCVdist( xTemp );
//            xiMinus = SandCVCMap.GetValue();
////            cout<<"-cv = [ "<<xiMinus.at(0)<<" , "<<xiMinus.at(1)<<"]"<<endl;
////            cout<<"________________________________________"<<endl;

//            for ( jD = 0; jD < ldim; ++jD )
//                jacFinite.at(iN*hdim*ldim+jD*hdim+iD) = ( xiPlus.at(jD) - xiMinus.at(jD) ) / (2*deltaX);
//        }
    }

//    cout<<"jacs ="<<endl;
//    for ( iD = 0; iD < hdim; ++iD)
//        cout<<"("<<jacSandCVCMap.at(iD)<<","<<jacFinite.at(iD)<<") ";
//    cout<<endl;
//    for ( iD = hdim; iD < ldim*hdim; ++iD)
//        cout<<"("<<jacSandCVCMap.at(iD)<<","<<jacFinite.at(iD)<<") ";
//    cout<<endl;

    time_point<system_clock> end = system_clock::now();
    duration<double> elapsed = end - start;
    cout<<"Elapsed time............................."<<elapsed.count()<<"seconds"<<endl;

    //__________________________________________________________________________
    // [2-2] Save the jacobian matrices
    //
    IO.WritePointSet(prefix+"xiProjected.txt"    , xiProjected, nTests, ldim, 8);
    IO.WritePointSet(prefix+"jacobian_SandCVCMap.txt", jacSandCVCMap  , nTests, hdim*ldim, 8);
    IO.WritePointSet(prefix+"jacobian_finite.txt", jacFinite  , nTests, hdim*ldim, 8);


    return 0;
}
