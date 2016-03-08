/*=========================================================================
 * File         : DihedralAngles.h
 * Module       : Simulation
 * Copyright    : (C)opyright 2013
 *                See COPYRIGHT statement in top level directory.
 * Authors      : Behrooz Hashemian
 * Created      : Mon Dec 16, 2013  03:40PM
 * Last modified: Mon Dec 16, 2013  03:41PM
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

//C++ Standard Libraries
#include <vector>
#include <numeric>
#include <limits>
#include <exception>
#include <stdexcept>
#include <algorithm>
#include <functional>
#include <cmath>
#include <iostream>
#include <fstream>
using namespace std;

/** \brief This class provides values and gradients of dihedral angles
 *
 */
class DihedralAngles
{
public :
    /** Default Constructor */
    DihedralAngles();

    /** Default Destructor */
    ~DihedralAngles() = default;

    /** Copy constructor */
    DihedralAngles(const DihedralAngles& inDihedralAngles) = delete;

    /** Assignmet operator */
    void operator=( const DihedralAngles& ) = delete;

    /** Compute dihedral angles values and gradients */
    void ComputeDihedralAngles( const vector<double> &x_new );

    /** Calculate auxiliary vectors */
    void CalculateAuxiliaryVectors(const vector<int> &AtomIndices );

    /** Compute dihedral angles values */
    //@{
    double ComputeDihedralAngleValue( );
    double ComputeDihedralAngleValue( int *initIndex );
    //@}

    /** Compute dihedral angles gradient */
    //@{
    double *ComputeDihedralAngleGradient( );
    double *ComputeDihedralAngleGradient( int *initIndex );
    //@}

    /** Compute dihedral angles invgradient */
    //@{
    double *ComputeDihedralAngleInvGradient( );
    double *ComputeDihedralAngleInvGradient( int *initIndex );
    //@}

    void Update();

protected:
    bool    bOneSiteForce;
    int     ierr;
    int     nDim;
    int     nAtoms;
    int     _nDihedralAngles;
    vector<int   > dihedralIndices;
    vector<double> dihedralAngleValues;
    vector<double> dihedralAngleGradients;
    vector<double> dihedralAngleInvGradients;
    vector<double> gradientTemp;
    vector<double> invgradientTemp;
    vector<double> x_new;

    double  MIN_VECTOR_MODULUS;

    // Dihedral Angles Workspace

    double dot_r12_r32;
    double dot_r43_r32;
    double dot_r32_r32;
    double dot_r12_r12;
    double dot_r43_r43;
    vector<double> pos1  [3];
    vector<double> pos2  [3];
    vector<double> pos3  [3];
    vector<double> pos4  [3];
    vector<double> r12   [3];
    vector<double> r32   [3];
    vector<double> r43   [3];
    vector<double> crossM[3];
    vector<double> crossN[3];
    vector<double> projR [3];
    vector<double> projS [3];
    vector<double> Tea   [3];
    vector<double> jac1, jac2, jac3, jac4;
    vector<double> invjac1, invjac4;

};

