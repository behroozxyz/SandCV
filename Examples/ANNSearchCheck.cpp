#include <iostream>
#include <limits>
#include <typeinfo>
#include <exception>
#include <string>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <stdexcept>
#include "ANNSearch.h"
#include "mat.h"
using namespace std;
using namespace chrono;

void Row2ColMajor( int nRow, int nCol, const double *matRow, double *matCol );
void Col2RowMajor( int nRow, int nCol, const double *matCol, double *matRow );

int main(int argc, char* argv[])
{
    bool verbose = false;
    int i,j;
    time_point<system_clock> start, end;
    duration<double> elapsed;
    double min = numeric_limits<double>::min();
    double max = numeric_limits<double>::max();
    double eps = numeric_limits<double>::epsilon();
    cout<<"min:"<<min<<", max:"<<max<<" and eps = "<<eps<<endl;

    MATFile *input_mat;
    MATFile *output_mat;
    mxArray *x_input;
    string inputFile = "../../ala_noh.mat";

    //__________________________________________________________________________
    // Setting the configurations data

    cout<<"Reading file "<<inputFile<<"..."<<endl;
    input_mat = matOpen(inputFile.data(), "r");
    if ( input_mat == nullptr ) {
        cerr<<"Error reading file "<<inputFile<<endl;
        throw invalid_argument("File does not exist!");
    }

    x_input = matGetVariable(input_mat,"X");
    if (x_input == NULL) {
      printf("Error reading existing matrix LocalDouble\n");
      throw runtime_error("Cannot read theStatus!");
    }
    int     nFrames   = mxGetDimensions(x_input)[0];
    int     ndim      = mxGetDimensions(x_input)[1];
    double *p_x_nodes = mxGetPr(x_input);
    vector<double> x_nodes(nFrames*ndim);
    cout<<"pa1 = "<<p_x_nodes[0]<<endl;
    Col2RowMajor( nFrames, ndim, p_x_nodes, x_nodes.data() );
    cout<<"dimension = "<<ndim<<endl;
    cout<<"frames    = "<<nFrames<<endl;
    /* clean up before exit */
    mxDestroyArray(x_input);
    if ( matClose(input_mat) )
    {
        cout<<"Error closing file "<<inputFile<<endl;
        throw runtime_error("Error in closing mat file.");
    }

    //__________________________________________________________________________
    // x_new
//    cout<<"x_nodes: ";
//    for (auto iter=x_nodes.begin(); iter != x_nodes.end(); ++iter )
//        cout << *iter << "  ";
//    cout<<endl;
    ANNSearch searcher(1, ndim, x_nodes);
    vector<double> query(ndim, 0.0);
    vector<int> idx;
    start = system_clock::now();
    for ( i = 0; i < nFrames-1; ++i )
    {
        copy( x_nodes.begin()+i*ndim, x_nodes.begin()+(i+1)*ndim, query.begin() );
        idx = searcher.FindNearestNeighbors( query );
        if ( idx.at(0) != i )
            cout<<"idx = "<<idx.at(0)<<" is not right!"<<endl;
//        for (auto iter : query )
//            cout << iter << "  ";
//        cout<<endl;

    }
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
