/*=========================================================================
 * File         : Corral.h
 * Module       : Simulation
 * Copyright    : (C)opyright 2013
 *                See COPYRIGHT statement in top level directory.
 * Authors      : Behrooz Hashemian
 * Created      : Sun Aug 11, 2013  02:14PM
 * Last modified: Sun Aug 11, 2013  02:14PM
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

//Sand Libraries
#include "PointSetIO.h"
#include "ANNSearch.h"

//Intel Math Kernel Library
#include <mkl.h>

// C++ Standard Libraries
#include <vector>
#include <string>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <limits>
using namespace std;

/** \class Corral
 * \brief This class provides all necessary functions to generate a continius corral function
 * around the input embedded data in a way that its gradients (force) near the borders
 * keep the new configurations inside the corral.
 *
 */

class Corral
{
public:
    /** Default Constructor */
    Corral() = default;

    /** Constructor with initialization */
    Corral(string filename, double corralHeight = 1.0);
    Corral(const vector<double> &xiCorral, int dim, double spacing, double corralHeight = 1.0);

    /** Default Destructor */
    ~Corral() = default;

    /** Copy constructor (explicitly deleted) */
    Corral( const Corral& ) = delete;

    /** Assignmet operator (explicitly deleted) */
    void operator=( const Corral& ) = delete;

    /** Calculate corral value and gradients */
    void CalculateCorral(const vector<double> &xiNew, vector<double> &gradient);

    vector<double> GetMarkers() { return _xiMarkers; }

private:
    vector<double> _xiMarkers;
    int    _dim            = 0;
    int    _nMarkers       = 0;
    int    _nInsideMarkers = 0;
    double _spacing        = 0.0;
    double _gammaPU        = 1.0;
    double _betaPU         = 0.0;
    double _rangePU        = 0.0;
    double _tolPU          = 1.0e-12;
    double _corralHeight   = 1.0;

    double         _Z;
    double         _w;
    double         _wIn;
    vector<double> _diff;
    vector<double> _dZ;

    PointSetIO     _reader;
    shared_ptr<ANNSearch> _markerSearch;



    /** Generate the corral markers */
    void GenerateCorralMarkers();

    /** Check if parameters and point coordinates are properly set. */
    bool CheckSettings();

    /** Update the variables and make the class ready to use. */
    void Update() throw( invalid_argument );
};
