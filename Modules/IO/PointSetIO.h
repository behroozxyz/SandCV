/*=========================================================================
 * File         : PointSetIO.h
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
#pragma once

// C++ Standard Libraries
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <exception>
#include <stdexcept>
#include <iterator>
using namespace std;

/** \brief This class read point sets from a text file */

class PointSetIO
{
public:
    /** Default Constructor */
    PointSetIO() = default;

    /** Default Destructor */
    ~PointSetIO() = default;

    /** Copy constructor (explicitly deleted) */
    PointSetIO( const PointSetIO& ) = delete;

    /** Assignmet operator (explicitly deleted) */
    void operator=( const PointSetIO& ) = delete;

    /** Read data from a text file and store them in a double vector*/
    void ReadPointSet(const string &filename , vector<double> &pointSet,
                      int &nLines, int &nCols) throw(runtime_error);

    /** Read data from a text file and store them in a double vector*/
    void WritePointSet(const string &filename , const vector<double> &pointSet,
                       const int &nLines, const int &nCols, const int &precision = 15)
    throw (runtime_error,invalid_argument);
};
