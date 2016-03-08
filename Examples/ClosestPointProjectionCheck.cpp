#include <algorithm>
#include <iostream>
#include <limits>
#include <typeinfo>
#include <exception>
#include <string>
#include <chrono>
#include <iomanip>
#include <stdexcept>
#include "ClosestPointProjection.h"
#include "ProcrustesSuperimposition.h"
#include "mat.h"
using namespace std;
using namespace chrono;

void Row2ColMajor( int nRow, int nCol, const double *matRow, double *matCol );
void Col2RowMajor( int nRow, int nCol, const double *matCol, double *matRow );

int main(int argc, char* argv[])
{
//    bool verbose = false;
    int i;
    time_point<system_clock> start, end;
    duration<double> elapsed;
    double min = numeric_limits<double>::min();
    double max = numeric_limits<double>::max();
    double eps = numeric_limits<double>::epsilon();
    cout<<"min:"<<min<<", max:"<<max<<" and eps = "<<eps<<endl;

    MATFile *inputMat;
    MATFile *outputMat;
    mxArray *xInput, *xOutput;
    string inputXn  = "../../xNodes.mat";    // X
    string inputXin = "../../xiNodes.mat";   // U, h_bin
    string inputXis = "../../xiSeeds72.mat"; // U
    string inputXs  = "../../xSamplesAligned.mat"; // U

    int     nFrames    = 0;
    double *pInputMat  = nullptr;
    double *pOutputMat = nullptr;

    //--------------------------------------------------------------------------
    // [1] Load point sets
    //--------------------------------------------------------------------------
    //__________________________________________________________________________
    // [1-0] Load high-dimensional node points
    //
    cout<<"Reading high-dimensional node points: ["<<inputXn<<"] ..."<<endl;
    inputMat = matOpen(inputXn.data(), "r");
    if ( inputMat == nullptr ) {
        cerr<<"Error reading file "<<inputXn<<endl;
        throw invalid_argument("File does not exist!");
    }

    xInput = matGetVariable(inputMat,"X");
    if (xInput == NULL) {
      cerr<<"Error reading matrix X from "<<inputXn<<endl;
      throw runtime_error("Cannot read theStatus!");
    }
    nFrames   = mxGetDimensions(xInput)[0];
    int hdim  = mxGetDimensions(xInput)[1];
    pInputMat = mxGetPr(xInput);
    vector<double> xNodes(nFrames*hdim);
    Col2RowMajor( nFrames, hdim, pInputMat, xNodes.data() );
    cout<<"High-dimensional nodes are loaded."<<endl;
    cout<<"\t dimension = "<<hdim<<endl;
    cout<<"\t frames    = "<<nFrames<<endl;
    mxDestroyArray(xInput);
    if ( matClose(inputMat) )
    {
        cerr<<"Error closing file "<<inputXn<<endl;
        throw runtime_error("Error in closing mat file.");
    }

    //__________________________________________________________________________
    // [1-1] Load low-dimensional node points, spacing
    //
    cout<<"Reading low-dimensional node points: ["<<inputXin<<"] ..."<<endl;
    inputMat = matOpen(inputXin.data(), "r");
    if ( inputMat == nullptr ) {
        cerr<<"Error reading file "<<inputXin<<endl;
        throw invalid_argument("File does not exist!");
    }

    xInput = matGetVariable(inputMat,"U");
    if (xInput == NULL) {
      cerr<<"Error reading matrix X from "<<inputXin<<endl;
      throw runtime_error("Cannot read theStatus!");
    }
    nFrames   = mxGetDimensions(xInput)[0];
    int ldim  = mxGetDimensions(xInput)[1];
    pInputMat = mxGetPr(xInput);
    vector<double> xiNodes(nFrames*ldim);
    Col2RowMajor( nFrames, ldim, pInputMat, xiNodes.data() );
    mxDestroyArray(xInput);
    // load spacing
    xInput = matGetVariable(inputMat,"h_bin");
    pInputMat = mxGetPr(xInput);
    double spacing = pInputMat[0];
    mxDestroyArray(xInput);
    cout<<"Low-dimensional nodes are loaded."<<endl;
    cout<<"\t dimension = "<<ldim<<endl;
    cout<<"\t frames    = "<<nFrames<<endl;
    cout<<"\t spacing   = "<<spacing<<endl;
    if ( matClose(inputMat) )
    {
        cerr<<"Error closing file "<<inputXin<<endl;
        throw runtime_error("Error in closing mat file.");
    }

    //__________________________________________________________________________
    // [1-2] Load low-dimensional seed points
    //
    cout<<"Reading low-dimensional seed points: ["<<inputXis<<"] ..."<<endl;
    inputMat = matOpen(inputXis.data(), "r");
    if ( inputMat == nullptr ) {
        cerr<<"Error reading file "<<inputXis<<endl;
        throw invalid_argument("File does not exist!");
    }

    xInput = matGetVariable(inputMat,"U");
    if (xInput == NULL) {
      cerr<<"Error reading matrix X from "<<inputXis<<endl;
      throw runtime_error("Cannot read theStatus!");
    }
    nFrames = mxGetDimensions(xInput)[0];
    if (mxGetDimensions(xInput)[1] != ldim)
    {
        cerr<<"Dimension of low-dimenional seed points does not match with node points."<<endl;
        throw runtime_error("dimension mismatch!");
    }
    pInputMat = mxGetPr(xInput);
    vector<double> xiSeeds(nFrames*ldim);
    Col2RowMajor( nFrames, ldim, pInputMat, xiSeeds.data() );
    cout<<"Low-dimensional seeds are loaded."<<endl;
    cout<<"\t dimension = "<<ldim<<endl;
    cout<<"\t frames    = "<<nFrames<<endl;
    mxDestroyArray(xInput);
    if ( matClose(inputMat) )
    {
        cerr<<"Error closing file "<<inputXis<<endl;
        throw runtime_error("Error in closing mat file.");
    }

    //__________________________________________________________________________
    // [1-3] Load high-dimensional sample points
    //
    cout<<"Reading high-dimensional sample points: ["<<inputXs<<"] ..."<<endl;
    inputMat = matOpen(inputXs.data(), "r");
    if ( inputMat == nullptr ) {
        cerr<<"Error reading file "<<inputXs<<endl;
        throw invalid_argument("File does not exist!");
    }

    xInput = matGetVariable(inputMat,"X");
    if (xInput == NULL) {
      cerr<<"Error reading matrix X from "<<inputXs<<endl;
      throw runtime_error("Cannot read theStatus!");
    }
    nFrames   = mxGetDimensions(xInput)[0];
    if (mxGetDimensions(xInput)[1] != hdim)
    {
        cerr<<"Dimension of sample points does not match with node points."<<endl;
        throw runtime_error("dimension mismatch!");
    }
    pInputMat = mxGetPr(xInput);
    vector<double> xSamples(nFrames*hdim);
    Col2RowMajor( nFrames, hdim, pInputMat, xSamples.data() );
    mxDestroyArray(xInput);

    xInput = matGetVariable(inputMat,"dihedrals");
    if (xInput == NULL) {
      cerr<<"Error reading dihedrals from "<<inputXs<<endl;
      throw runtime_error("Cannot read dihedrals!");
    }
    if (mxGetDimensions(xInput)[0] != nFrames)
    {
        cerr<<"Number of dihedrals does not match with number of sample points."<<endl;
        throw runtime_error("number of points mismatch!");
    }
    if (mxGetDimensions(xInput)[1] != ldim)
    {
        cerr<<"Dimension of dihedrals does not match with node points."<<endl;
        throw runtime_error("dimension mismatch!");
    }
    pInputMat = mxGetPr(xInput);
    vector<double> dihedrals(nFrames*ldim);
    Col2RowMajor( nFrames, ldim, pInputMat, dihedrals.data() );
    mxDestroyArray(xInput);

    cout<<"High-dimensional sample points and the related dihedral values are loaded."<<endl;
    cout<<"\t dimension = "<<hdim<<endl;
    cout<<"\t frames    = "<<nFrames<<endl;
    if ( matClose(inputMat) )
    {
        cerr<<"Error closing file "<<inputXs<<endl;
        throw runtime_error("Error in closing mat file.");
    }
    pInputMat = nullptr;
    //--------------------------------------------------------------------------
    // [2] Alignment and closest point projection
    //--------------------------------------------------------------------------
    //__________________________________________________________________________
    // [2-0] Initialize closest point projection
    //
    ClosestPointProjection CPP(ldim, hdim, spacing, xiNodes, xNodes, xiSeeds);

    //__________________________________________________________________________
    // [2-1] Save seed points
    //
    /* save x_seed */
    auto x_seed  = CPP.GetHighSeeds();
    outputMat = matOpen("../../out_xSeeds.mat","w");
    int nSeeds = x_seed.size()/hdim;
    xOutput   = mxCreateDoubleMatrix(nSeeds,hdim,mxREAL);
    pOutputMat = mxGetPr(xOutput);
    Row2ColMajor( nSeeds, hdim, x_seed.data(), pOutputMat );
    int status = matPutVariable(outputMat, "X", xOutput);
    mxDestroyArray(xOutput);
    if ( matClose(outputMat) )
    {
        cerr<<"Error closing file x_seeds"<<endl;
        throw runtime_error("Error in closing mat file.");
    }


    /* save xi_seed*/
    auto xi_seed = CPP.GetLowSeeds();
    outputMat = matOpen("../../out_xiSeeds.mat","w");
    nSeeds = xi_seed.size()/ldim;
    xOutput   = mxCreateDoubleMatrix(nSeeds,ldim,mxREAL);
    pOutputMat = mxGetPr(xOutput);
    Row2ColMajor( nSeeds, ldim, xi_seed.data(), pOutputMat );
    status = matPutVariable(outputMat, "U", xOutput);
    mxDestroyArray(xOutput);
    if ( matClose(outputMat) )
    {
        cerr<<"Error closing file xi_seeds"<<endl;
        throw runtime_error("Error in closing mat file.");
    }


    //__________________________________________________________________________
    // [2-2] Find closest point projection of some aligned sample points
    //
    int correct = 0;
    int all     = 0;
//    double err;
    vector<double> xi(ldim);
    vector<double> xiProjected;
    vector<double>::iterator iter;
    start = system_clock::now();
    for ( i = 0; i < nFrames; ++i )
    {
        if ( fabs(dihedrals.at(i*ldim)) > 175.0 || fabs(dihedrals.at(i*ldim+1)) > 175.0 )
            continue;
        else
        {
            ++all;
            iter = xSamples.begin() + i*hdim;
            try
            {
                //            cout<<"xSample = [";
                //            for ( j = 0; j < hdim; ++j )
                //                cout<<*(iter+j)<<" ";
                //            cout<<"]"<<endl;
                xi = CPP.FindClosestPoint( iter, iter+hdim );
            }
            catch(const runtime_error &e)
            {

                cout<<"Expected("<<i<<") = "<<dihedrals.at(i*ldim)<<" , "<<dihedrals.at(i*ldim+1)<<""<<endl;
//                cout<<"Failiure("<<i<<") = "<<xi.at(0)<<" , "<<xi.at(1)<<""<<endl;
                cerr<<__func__<<" >> "<<e.what()<<endl;
                continue;
            }
            ++correct;
        }
//        err = sqrt ( (xi.at(0)-dihedrals.at(i*ldim)  ) * (xi.at(0)-dihedrals.at(i*ldim)  )
//            + (xi.at(1)-dihedrals.at(i*ldim+1)) * (xi.at(1)-dihedrals.at(i*ldim+1)) );
//        if ( err > 1e-8)
//            cout<<"error("<<i<<") ="<<err<<endl;
//        xiProjected.insert( xiProjected.end(), xi.begin(), xi.end() );
    }
    cout<<"correct = "<<correct<<" / "<<all<<endl;
    end = system_clock::now();
    elapsed = end - start;
    cout<<"Elapsed time............................."<<elapsed.count()<<"seconds"<<endl;

    return 0;
}

void Row2ColMajor( int nRow, int nCol, const double *matRow, double *matCol )
{
    int i, j;
    for ( i = 0; i < nRow; ++i )
        for ( j = 0; j < nCol; ++j )
            matCol[j*nRow+i] = matRow[i*nCol+j];
}


void Col2RowMajor( int nRow, int nCol, const double *matCol, double *matRow )
{
    int i, j;
    for ( i = 0; i < nCol; ++i )
        for ( j = 0; j < nRow; ++j )
            matRow[j*nCol+i] = matCol[i*nRow+j];
}
