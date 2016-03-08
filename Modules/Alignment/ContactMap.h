/*=========================================================================
 * File         : ContactMap.h
 * Module       : Alignment
 * Copyright    : (C)opyright 2011-2013
 *                See COPYRIGHT statement in top level directory.
 * Authors      : Behrooz Hashemian
 * Created      : Fri Nov 22, 2013  01:47PM
 * Last modified: Fri Nov 22, 2013  01:47PM
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

//Intel Math Kernel Library
#include <mkl.h>

//C++ Standard Libraries
#include <stdexcept>
#include <cmath>
#include <vector>
#include <numeric>
#include <iostream>
#include <algorithm>
using namespace std;

/** \class ContactMap
 * \brief Generate a smooth contact map of a molecular configuration
 *
 * \b References: \n
 *
 */
class ContactMap
{
public :
    /** Default Constructor */
    ContactMap() = default;

    /** Destructor */
    ~ContactMap() = default;

    /** Copy constructor */
    ContactMap(const ContactMap& inCMap);

    /** Assignmet operator (explicitly deleted) */
    void operator=( const ContactMap& ) = delete;

    ContactMap(int sdim, int nAtoms, double d0);


    /** Obtain the derivatives of contact map elements after calling ComputeDerivatives */
    const vector<double> &GetdAdX( ) const
    {  return this->_dcmap;  }

    /** Obtain the aligned configuration after calling ComputeAlignment or ComputeDerivatives */
    const vector<double> &GetAlignedX( ) const
    {  return this->_cmap;  }

    /** Perform the alignment and calculate the derivatives if applicable
     * \param [in]  x_new
     * \param [out] x_aligned
     * \return error identifier */
    vector<double> ComputeAlignment( const vector<double> &x_new_full );
    vector<double> ComputeDerivatives(const vector<double> &x_new_full );
protected:
    int    _sdim   = 3;
    int    _nAtoms = 0;
    double _d0inv     = 0.0;
    double _d      = 0.0;
    vector<double> _rij;
    vector<double> _cmap;
    vector<double> _dcmap;
};
