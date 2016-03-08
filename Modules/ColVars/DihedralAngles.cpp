/*=========================================================================

Module    : Solid Mechanics
File      : DihedralAngles.cpp
Copyright : (C)opyright 2012
            See COPYRIGHT statement in top level directory.
Authors   : B. Hashemian
Modified  :
Purpose   : Smooth Contact Map
Date      :
Version   :
Changes   :

    This software is distributed WITHOUT ANY WARRANTY; without even
    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
    PURPOSE.  See the above copyright notices for more information.
=========================================================================*/
#include "DihedralAngles.h"

DihedralAngles::DihedralAngles( )
{
    this->MIN_VECTOR_MODULUS     = 1.0E-12;
}

void DihedralAngles::Update( )
{
    if ( this->CheckSettings() )
    {
        printf("ERROR:ContactMap settings are wrong\n");
        this->ierr = 1;
        return;
    }

    this->dihedralAngleValues    = new double [this->_nDihedralAngles           ];
    this->dihedralAngleGradients = new double [this->_nDihedralAngles*this->nDim];
    this->dihedralAngleInvGradients = new double [this->_nDihedralAngles*this->nDim];
    this->gradientTemp           = new double [this->nDim];
    this->invgradientTemp        = new double [this->nDim];

    return;
}


//check the setting before to compute
int DihedralAngles::CheckSettings( )
{
    if ( this->nAtoms < 4 )
    {
        printf("ERROR::DihedralAngles::CheckSettings :\n");
        printf("Number ofatoms has not been set correctly (%d < 4)\n", this->nAtoms);
        return 1;
    }
    else
    {
        this->nDim = 3*this->nAtoms;
    }
    if ( this->_nDihedralAngles < 1 )
    {
        printf("ERROR::DihedralAngles::CheckSettings :\n");
        printf("Number of dihedral angles has not been set correctly (%d < 1)\n", this->_nDihedralAngles);
        return 1;
    }
    if ( this->dihedralIndices == NULL )
    {
        printf("ERROR::DihedralAngles::CheckSettings :\n");
        printf("Dihedral angles' index vector has not been set.\n");
        return 1;
    }
    else
    {
        int maxIndex = smArray::Max(4*this->_nDihedralAngles, this->dihedralIndices);
        int minIndex = smArray::Min(4*this->_nDihedralAngles, this->dihedralIndices);
        if (maxIndex > this->nAtoms-1)
        {
            printf("ERROR::DihedralAngles::CheckSettings :\n");
            printf("Indices of dihedral angles is out of range [(maxIndex)%d > (nAtoms-1)%d]\n", maxIndex, this->nAtoms-1);
            return 1;
        }
        if (minIndex < 0)
        {
            printf("ERROR::DihedralAngles::CheckSettings :\n");
            printf("Indices of dihedral angles is out of range [(minIndex)%d < 0]\n", minIndex);
            return 1;
        }
    }

    this->Clear();

    return 0;
}



void DihedralAngles::ComputeDihedralAngles(const vector<double> &x_new )
{
    int *initIndex;
    this->x_new = x_new;
    for ( int i = 0; i < _nDihedralAngles; ++i)
    {
        initIndex = &this->dihedralIndices[4*i];
        this->CalculateAuxiliaryVectors( initIndex );
        this->dihedralAngleValues[i] = this->ComputeDihedralAngleValue();
        smArray::Copy(this->nDim, this->ComputeDihedralAngleGradient(),
                      &this->dihedralAngleGradients[i*this->nDim]);
        smArray::Copy(this->nDim, this->ComputeDihedralAngleInvGradient(),
                      &this->dihedralAngleInvGradients[i*this->nDim]);
    }
}



void DihedralAngles::CalculateAuxiliaryVectors(const vector<int> &AtomIndices)
{
    double *pos1, *pos2, *pos3, *pos4;
    int atom1, atom2, atom3, atom4;

    /* torsion angle atoms */
    atom1 = AtomIndices[0];
    atom2 = AtomIndices[1];
    atom3 = AtomIndices[2];
    atom4 = AtomIndices[3];

    pos1 = &this->x_new[3 * atom1];
    pos2 = &this->x_new[3 * atom2];
    pos3 = &this->x_new[3 * atom3];
    pos4 = &this->x_new[3 * atom4];


    smArray::Fill(this->nDim, 0.0, this->gradientTemp);
    this->jac1 = &this->gradientTemp[3 * atom1];
    this->jac2 = &this->gradientTemp[3 * atom2];
    this->jac3 = &this->gradientTemp[3 * atom3];
    this->jac4 = &this->gradientTemp[3 * atom4];

    smArray::Fill(this->nDim, 0.0, this->invgradientTemp);
    this->invjac1 = &this->invgradientTemp[3 * atom1];
    this->invjac4 = &this->invgradientTemp[3 * atom4];

    /* torsion angle vectors */
    this->r12[0] = pos1[0] - pos2[0];
    this->r12[1] = pos1[1] - pos2[1];
    this->r12[2] = pos1[2] - pos2[2];

    this->r32[0] = pos3[0] - pos2[0];
    this->r32[1] = pos3[1] - pos2[1];
    this->r32[2] = pos3[2] - pos2[2];

    this->r43[0] = pos4[0] - pos3[0];
    this->r43[1] = pos4[1] - pos3[1];
    this->r43[2] = pos4[2] - pos3[2];


    /* cross product vectors between vectors this->r43 and this->r32 (N vector) */
    this->crossN[0] = this->r43[1] * this->r32[2] - this->r43[2] * this->r32[1];
    this->crossN[1] = this->r43[2] * this->r32[0] - this->r43[0] * this->r32[2];
    this->crossN[2] = this->r43[0] * this->r32[1] - this->r43[1] * this->r32[0];

    /* vector dot products and factors */
    this->dot_r12_r32 = this->r12[0] * this->r32[0] + this->r12[1] * this->r32[1] + this->r12[2] * this->r32[2];
    this->dot_r43_r32 = this->r43[0] * this->r32[0] + this->r43[1] * this->r32[1] + this->r43[2] * this->r32[2];
    this->dot_r32_r32 = this->r32[0] * this->r32[0] + this->r32[1] * this->r32[1] + this->r32[2] * this->r32[2];

}



double DihedralAngles::ComputeDihedralAngleValue(int *initIndex)
{
    this->CalculateAuxiliaryVectors(initIndex);
    return this->ComputeDihedralAngleValue();
}



double DihedralAngles::ComputeDihedralAngleValue()
{
    double cos_angle, sign_angle;

    /* R and S projection vectors perpendicular
         * to line 32 of lines 12 and 43 respectively */
    this->projR[0] = this->r12[0] - this->r32[0] * this->dot_r12_r32 / this->dot_r32_r32;
    this->projR[1] = this->r12[1] - this->r32[1] * this->dot_r12_r32 / this->dot_r32_r32;
    this->projR[2] = this->r12[2] - this->r32[2] * this->dot_r12_r32 / this->dot_r32_r32;

    this->projS[0] = this->r43[0] - this->r32[0] * this->dot_r43_r32 / this->dot_r32_r32;
    this->projS[1] = this->r43[1] - this->r32[1] * this->dot_r43_r32 / this->dot_r32_r32;
    this->projS[2] = this->r43[2] - this->r32[2] * this->dot_r43_r32 / this->dot_r32_r32;

    /* torsion angle */
    sign_angle = SIGN( 1, this->r12[0] * this->crossN[0] + this->r12[1] * this->crossN[1] + this->r12[2] * this->crossN[2] );

    cos_angle  = ( this->projR[0] * this->projR[0] + this->projR[1] * this->projR[1] + this->projR[2] * this->projR[2] ) *
            ( this->projS[0] * this->projS[0] + this->projS[1] * this->projS[1] + this->projS[2] * this->projS[2] );
    cos_angle  = ( this->projR[0] * this->projS[0] + this->projR[1] * this->projS[1] + this->projR[2] * this->projS[2] ) / sqrt(cos_angle);
    cos_angle  = MIN( 1.0, MAX( -1.0, cos_angle ) );

    return sign_angle * acos( cos_angle );
}



double *DihedralAngles::ComputeDihedralAngleGradient(int *initIndex)
{
    this->CalculateAuxiliaryVectors(initIndex);
    return this->ComputeDihedralAngleGradient();
}



double *DihedralAngles::ComputeDihedralAngleGradient()
{
    double mod_M, mod_N;

    /* cross product vectors between vectors this->r12 and this->r32 (M vector)
     */
    this->crossM[0] = this->r12[1] * this->r32[2] - this->r12[2] * this->r32[1];
    this->crossM[1] = this->r12[2] * this->r32[0] - this->r12[0] * this->r32[2];
    this->crossM[2] = this->r12[0] * this->r32[1] - this->r12[1] * this->r32[0];
    /* norm of vector this->crossM */
    mod_M = sqrt( this->crossM[0] * this->crossM[0] + this->crossM[1] * this->crossM[1] + this->crossM[2] * this->crossM[2] );
    mod_M = MAX( MIN_VECTOR_MODULUS, mod_M );
    /* norm of vector this->crossN (computed in ComputeDihedralAngleVectors) */
    mod_N = sqrt( this->crossN[0] * this->crossN[0] + this->crossN[1] * this->crossN[1] + this->crossN[2] * this->crossN[2] );
    mod_N = MAX( MIN_VECTOR_MODULUS, mod_N );


    /* torsion jacobians this->jac1 and this->jac4... */
    this->jac1[0] =  this->crossM[0] * sqrt( this->dot_r32_r32 ) / ( mod_M * mod_M );
    this->jac1[1] =  this->crossM[1] * sqrt( this->dot_r32_r32 ) / ( mod_M * mod_M );
    this->jac1[2] =  this->crossM[2] * sqrt( this->dot_r32_r32 ) / ( mod_M * mod_M );

    this->jac4[0] = -this->crossN[0] * sqrt( this->dot_r32_r32 ) / ( mod_N * mod_N );
    this->jac4[1] = -this->crossN[1] * sqrt( this->dot_r32_r32 ) / ( mod_N * mod_N );
    this->jac4[2] = -this->crossN[2] * sqrt( this->dot_r32_r32 ) / ( mod_N * mod_N );

    /* T vector to compute the central gradients */
    this->Tea[0] = this->jac1[0] * this->dot_r12_r32 / this->dot_r32_r32 +
            this->jac4[0] * this->dot_r43_r32 / this->dot_r32_r32;

    this->Tea[1] = this->jac1[1] * this->dot_r12_r32 / this->dot_r32_r32 +
            this->jac4[1] * this->dot_r43_r32 / this->dot_r32_r32;

    this->Tea[2] = this->jac1[2] * this->dot_r12_r32 / this->dot_r32_r32 +
            this->jac4[2] * this->dot_r43_r32 / this->dot_r32_r32;


    /* ...Central jacobians this->jac2 and this->jac3 */
    this->jac2[0] = -this->jac1[0] + this->Tea[0];
    this->jac2[1] = -this->jac1[1] + this->Tea[1];
    this->jac2[2] = -this->jac1[2] + this->Tea[2];

    this->jac3[0] = -this->jac4[0] - this->Tea[0];
    this->jac3[1] = -this->jac4[1] - this->Tea[1];
    this->jac3[2] = -this->jac4[2] - this->Tea[2];


    return this->gradientTemp;
}



double *DihedralAngles::ComputeDihedralAngleInvGradient(int *initIndex)
{
    this->CalculateAuxiliaryVectors(initIndex);
    return this->ComputeDihedralAngleInvGradient();
}



double *DihedralAngles::ComputeDihedralAngleInvGradient()
{
    double mod_M, mod_N;
    double fact1, fact4;

    this->dot_r12_r12 = this->r12[0] * this->r12[0] + this->r12[1] * this->r12[1] + this->r12[2] * this->r12[2];
    this->dot_r43_r43 = this->r43[0] * this->r43[0] + this->r43[1] * this->r43[1] + this->r43[2] * this->r43[2];

    /* cross product vectors between vectors this->r12 and this->r32 (M vector)
     */
    this->crossM[0] = this->r12[1] * this->r32[2] - this->r12[2] * this->r32[1];
    this->crossM[1] = this->r12[2] * this->r32[0] - this->r12[0] * this->r32[2];
    this->crossM[2] = this->r12[0] * this->r32[1] - this->r12[1] * this->r32[0];
    /* norm of vector this->crossM */
    mod_M = sqrt( this->crossM[0] * this->crossM[0] + this->crossM[1] * this->crossM[1] + this->crossM[2] * this->crossM[2] );
    mod_M = MAX( MIN_VECTOR_MODULUS, mod_M );
    /* norm of vector this->crossN (computed in ComputeDihedralAngleVectors) */
    mod_N = sqrt( this->crossN[0] * this->crossN[0] + this->crossN[1] * this->crossN[1] + this->crossN[2] * this->crossN[2] );
    mod_N = MAX( MIN_VECTOR_MODULUS, mod_N );

    fact1 = std::sqrt( this->dot_r12_r12 - this->dot_r12_r12 * this->dot_r12_r32 * this->dot_r12_r32
                      / (this->dot_r32_r32 * this->dot_r12_r12) );
    fact4 = std::sqrt( this->dot_r43_r43 - this->dot_r43_r43 * this->dot_r43_r32 * this->dot_r43_r32
                      / (this->dot_r32_r32 * this->dot_r43_r43) );

    /* torsion jacobians this->jac1 and this->jac4... */
    this->invjac1[0] =  fact1 * this->crossM[0] / mod_M;
    this->invjac1[1] =  fact1 * this->crossM[1] / mod_M;
    this->invjac1[2] =  fact1 * this->crossM[2] / mod_M;

    if (!this->bOneSiteForce)
    {
        this->invjac4[0] = fact4 * this->crossN[0] * mod_N;
        this->invjac4[1] = fact4 * this->crossN[1] * mod_N;
        this->invjac4[2] = fact4 * this->crossN[2] * mod_N;
    }

    return this->invgradientTemp;
}


