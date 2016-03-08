/*=========================================================================
 * File         : PointSetIO.cpp
 * Module       : IO
 * Copyright    : (C)opyright 2013
 *                See COPYRIGHT statement in top level directory.
 * Authors      : Behrooz Hashemian
 * Created      : Mon Aug 19, 2013  03:59PM
 * Last modified: Mon Aug 19, 2013  03:59PM
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
#include "PointSetIO.h"

void PointSetIO::ReadPointSet(const string &filename, vector<double> &pointSet,
                              int &nLines, int &nCols)
throw (runtime_error)
{
    int      nColsCheck = 0;
    double   elementVal = 0.0;
    string   line;
    ifstream inFile( filename );

    cout<<"Opening "<<filename<<" ..."<<endl;
if (!inFile.is_open())
        throw runtime_error("Error in opening file.");

    getline(inFile, line);
    istringstream iss(line);
    for ( nCols = 0; iss>>elementVal; ++nCols)
        pointSet.push_back(elementVal);

    for( nLines = 1; getline(inFile, line); ++nLines)
    {
        istringstream iss(line);
        for ( nColsCheck = 0; iss>>elementVal; ++nColsCheck)
            pointSet.push_back(elementVal);
        if(nColsCheck != nCols)
            throw runtime_error("Number of columns in the file is not consistent!");
    }

    cout<<"\tLoaded: nLines = "<<nLines<<", nCols = "<<nCols<<endl;
}


void PointSetIO::WritePointSet(const string &filename, const vector<double> &pointSet,
                               const int &nLines, const int &nCols, const int &precision)
throw (runtime_error,invalid_argument)
{
    ofstream outFile(filename);
    if (!outFile.is_open())
    {
        cerr<<__func__<<" >> Cannot open/create ["<<filename<<"] for writting."<<endl;
        throw runtime_error("Error in creating file.");
    }

    int nPoints = pointSet.size();

    if (nPoints == 0)
    {
        cerr<<__func__<<" >> The point set vector is empty"<<endl;
        throw invalid_argument("vector<double> pointSet is empty.");
    }

    if (nPoints < nLines*nCols)
    {
        cerr<<__func__<<" >> Number of points in the point set vector is less than"
              " the product of the given number of rows and number of columns."<<endl;
        throw invalid_argument("size of pointSet is less than nLines*nColumns");
    }
    if (nPoints > nLines*nCols)
        /* Number of points in the point set vector is more than
           the product of the given number of rows and number of columns. */
        cerr<<__func__<<" >> Partial writing!"<<endl;

    vector<double>::const_iterator pointIter;
    for ( int iLine = 0; iLine < nLines; ++iLine )
    {
        for ( pointIter = pointSet.begin()+iLine*nCols;
              pointIter != pointSet.begin()+(iLine+1)*nCols-1; ++pointIter )
            outFile<<scientific<<setprecision(precision)<<*pointIter<<"\t";
        outFile<<scientific<<setprecision(precision)<<*pointIter<<endl;
    }
}

