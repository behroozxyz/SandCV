/*=========================================================================
 * File         : ProcrustesSuperimposition.h
 * Module       : Alignment
 * Copyright    : (C)opyright 2011-2013
 *                See COPYRIGHT statement in top level directory.
 * Authors      : Behrooz Hashemian
 * Created      : Wed Jun 12, 2011  05:50PM
 * Last modified: Wed Jun 12, 2013  06:51PM
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

/** \class ProcrustesSuperimposition
 * \brief Perform Procrustes Analysis in order to superimpose a configuration with a reference one.
 *
 * \b References: \n
 *
 */
class ProcrustesSuperimposition
{
public :
    /** Default Constructor */
    ProcrustesSuperimposition() = default;

    /** Constructor with reference point initialization */
    ProcrustesSuperimposition(int dim, const vector<double> &x_ref_full );

    /** Constructor with a subset of reference point initialization */
    ProcrustesSuperimposition(int dim, const vector<double> &x_ref_full, const vector<int> &activePointsList );

    /** Destructor */
    ~ProcrustesSuperimposition() = default;

    /** Copy constructor */
    ProcrustesSuperimposition(const ProcrustesSuperimposition& inPS);

    /** Assignmet operator (explicitly deleted) */
    void operator=( const ProcrustesSuperimposition& ) = delete;


    /** Spatial dimension for the reference point set*/
    int  GetDimension() const    { return _sdim; }

    /** Get number of points of the reference configuration*/
    int  GetNumberOfPoints( ) const { return _nFullAtoms; }

    /** Set the options for scaling and reflecting*/
    void SetOptions( bool scaling, bool reflecting );

    //@{
    /** Set/Get the reference configuration's coordinate*/
    int SetReferenceConfiguration(const vector<double> &x_points_reference );
    int SetReferenceConfiguration(const vector<double> &x_points_reference, const vector<int> &involvedAtomes );

    const vector<double> &GetActivePoints( )           const { return _x_ref;    }
    const vector<double> &GetReferenceConfiguration( ) const { return _x_ref_full; }
    //@}


    //@{
    /** Obtain the Alignment Components (Rotation and Translation) */
    const vector<double> &GetRotation( ) const
    {  return this->_Rot;  }
    const vector<double> &GetTranslation( ) const
    {  return this->_Trans;  }
    //@}


    /** Obtain the Derivatives of rotation with respect to covariance matrix
     * Works Only for the Rigid 3D with Dervitives*/
    const vector<double> &GetdAdX( ) const
    {  return this->_dAdr;  }

    /** Obtain the aligned configuration after calling ComputeAlignment or ComputeDerivatives */
    const vector<double> &GetAlignedX( ) const
    {  return this->_x_aligned;  }

    /** Make the points zero mean*/
    vector<double> MakeZeroMean( vector<double> &x_points );

    /** Perform the alignment and calculate the derivatives if applicable
     * \param [in]  x_new
     * \param [out] x_aligned
     * \return error identifier */
    vector<double> ComputeAlignment( const vector<double> &x_new_full );
    vector<double> ComputeDerivatives(vector<double> x_new_full );

    /** Dirac delta */
    inline double delta( int i, int j ) const {  return ( i == j ) ? 1.0 : 0.0 ;  }

protected:
    int    _ierr          = 0;
    int    _sdim          = 3;
    int    _hdim          = 0;
    int    _nFullAtoms    = 0;
    double _nFullAtomsInv = 0;
    int    _nAtoms        = 0;
    double _nAtomsInv     = 0;

    vector<int>    _involvedAtoms;
    vector<double> _x_ref_full;
    vector<double> _x_ref;
    vector<double> _x_new_full;
    vector<double> _x_new;
    vector<double> _x_aligned;

    /** Update the variables and make the class ready to use. */
    int Update( ) throw( exception );

    /** Update the active points vectors*/
    int UpdateActivePoints( );
    int UpdateActivePoints(vector<double> x_new_full );

    /** Check if parameters and point coordinates are properly set. */
    bool CheckSettings();

private:
    double _norm_x_ref = 0.0;
    double _norm_x_new = 0.0;
    vector<int>    _ipiv;
    vector<double> _mean_x_ref;
    vector<double> _mean_x_new;
    vector<double> _M, _U, _S, _Vt;
    vector<double> _Rot;
    vector<double> _Trans;
    vector<double> _Omega;
    vector<double> _dRdr;
    vector<double> _dCdr;
    vector<double> _dAdr;

    // LAPACK functions' workspace
    char _job   = 'A';
};

