/*=========================================================================
 * File         : SandCVCMap.h
 * Module       : Simulation
 * Copyright    : (C)opyright 2013
 *                See COPYRIGHT statement in top level directory.
 * Authors      : Behrooz Hashemian
 * Created      : Sun Aug 11, 2013  02:54PM
 * Last modified: Sun Aug 11, 2013  02:54PM
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
#include "ProcrustesSuperimposition.h"
#include "ClosestPointProjection.h"
#include "ContactMap.h"

//Intel Math Kernel Library
#include <mkl.h>

//C++ Standard Libraries
#include <vector>
#include <memory>
#include <exception>
#include <stdexcept>
#include <algorithm>
#include <functional>
#include <cmath>
#include <limits>
#include <iostream>
using namespace std;

/** \brief This class provides smooth and nonlinear data-driven collective variable,
 * to be used in tandem with modecular dynamics simulations (see \cite Hashemian2013).
 *
 * Collective variables (CVs) help us rationalize molecular conformations and sample
 * complex free energy landscapes with molecular dynamics simulations.
 * Given their importance, there is need for systematic methods that effectively
 * identify CVs for complex systems. In recent years, nonlinear manifold learning
 * has shown its ability to automatically identify molecular collective behavior.
 * Unfortunately, these methods fail to provide a differentiable function mapping
 * high-dimensional configurations to their low-dimensional representation, as
 * required in enhanced sampling methods.
 * Smooth and nonlinear data-driven collective variables (SandCVCMap), overcome this
 * deficiency by calculating the missing mapping based on an ensemble representative
 * of molecular flexibility.
 */
class SandCVCMap
{
public :
    /** Default Constructor */
    SandCVCMap() = default;

    /** Constructors with Smooth Contact Map */
    SandCVCMap(int sdim, int ldim, int nAtoms, double spacing, double d0,
               const vector<double> &xiNodes, const vector<double> &xNodes,
               const vector<double> &xiSeeds);
    SandCVCMap(int sdim, int ldim, int nAtoms, double spacing, double d0,
               const vector<double> &xiNodes, const vector<double> &xNodes,
               const vector<double> &xiSeeds, const vector<double> &xiMarkers);


    /** Default Destructor */
    ~SandCVCMap() = default;

    /** Copy constructor */
    SandCVCMap(const SandCVCMap& inSandCVCMap);

    /** Assignmet operator */
    void operator=( const SandCVCMap& ) = delete;

    /** Compute the SandCVCMap value and Jacobian */
    void CalculateCV(const vector<double> &xNew) throw (runtime_error);
    void CalculateCV(const vector<double> &xNew, int stepNum) throw (runtime_error);

    /** Compute the SandCVCMap value and Jacobian and also consider the distance to
      * the manifold as an extra CV */
    void CalculateCVdist(const vector<double> &xNew) throw (runtime_error);

    /** Get the SandCVCMap's value */
    const vector<double>& GetValue() const { return _value; }

    /** Get the SandCVCMap's Jacobian */
    const vector<double>& GetJacobian() const { return _Jacobian; }

    /** Get the distance to the closest point projection */
    double GetDistance() const { return _cdist; }

    //@{
    /** Get the SandCVCMap's Mapping(X) and its derrivatives(DX) */
    const vector<double>& GetX()  const { return _projection.GetX();  }
    const vector<double>& GetDX() const { return _projection.GetDX(); }
    //@}

    //@{
    /** Some functions for debugging */
    const vector<double> &GetLowSeeds()  const { return _projection.GetLowSeeds();  }
    const vector<double> &GetHighSeeds() const { return _projection.GetHighSeeds(); }
    //@}

    /** Get the high-dimansional markers */
    const vector<double>& GetMarkers()  const { return _projection.GetHighMarkers();  }


protected :
    bool _bUseSeed = false;
    int  _stepNum  = 0;
    int  _ldim;             /// collective variable space dimension
    int  _hdim;             /// contact map space dimension
    int  _rdim;             /// non-aligned space dimension
    double _cdist;
    vector<double> _value;
    vector<double> _Jacobian;
    vector<double> _xMarkers;


    ContactMap                _alignment;
    ClosestPointProjection    _projection;
private :
};
