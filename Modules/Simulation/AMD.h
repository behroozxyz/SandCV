/*=========================================================================

Module    : Solid Mechanics
File      : smAcceleratedMD.h
Copyright : (C)opyright 2011++
            See COPYRIGHT statement in top level directory.
Authors   : B. Hashemian
Modified  :
Purpose   : Computing the forces of accelerated molecular dynamics
Date      :
Version   :
Changes   :

    This software is distributed WITHOUT ANY WARRANTY; without even
    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
    PURPOSE.  See the above copyright notices for more information.
=========================================================================*/

#ifndef __smAcceleratedMD_h
#define __smAcceleratedMD_h

// C standard lib
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <time.h>
#include <sys/timeb.h>

// SolMec lib
#include <smPointSetIO.h>
#include <smPointSet.h>
#include <smInputOutput.h>
#include <smMath.h>
#include <smDimensionReduction.h>
#include <smContactMap.h>
#include <smProcrustesAnalysis.h>
#include <smClosestPointProjectionLS.h>
#include <smContinuousCorral.h>
#include <smDihedralAngles.h>

/** \class smAcceleratedMD
 * \brief The smAcceleratedMD class
 *
 * \b References: \n
 *
 * \ingroup
 */
class SOLMEC_EXPORT smAcceleratedMD : public smDimensionReduction
{
public :
    // ----------------------------------------------------
    // Typedefs
    //-----------------------------------------------------
    typedef smAcceleratedMD Self;
    typedef smSmartPointer<Self> Pointer;

    typedef smEuclidianDistance Euclidian;
    typedef smNNAdjacencyStructure< Euclidian > NNAdjacencyStructure;
    // ----------------------------------------------------
    // Methods
    smNewMacro(Self) //Build-> :: New()

    /** Return the name of this class as a string. */
    const char *GetNameOfClass() const
    {	return "smAcceleratedMD";               }


    /** This function returns the error identifier
      * \return ierr 0: all computations finished succesfully \n
      *   1: an input parameter either has not been set or has been set bad \n
      */
    int GetError() const
    {   return this->ierr;                      }


    /** Set verbose ON or OFF (default is [OFF]) */
    void SetVerboseON( )
    {    this->bVerbose = true;                 }
    void SetVerboseOFF( )
    {    this->bVerbose = false;                }


    /** Set dihedral angles as the collective variables (default is [OFF]) */
    //@{
    void SetDihedralAnglesON( )
    {   this->bUseDihedralAngles = true;        }
    void SetDihedralAnglesOFF( )
    {   this->bUseDihedralAngles = false;       }
    //@}

    /** Set/Get numebr of atoms */
    //@{
    void SetNumberOfAtoms( int nAtoms )
    {   this->nAtoms = nAtoms;                  }
    int  GetNumberOfAtoms( ) const
    {   return this->nAtoms;                    }
    //@}

    /** Set the force calculation only for first atom ON/OFF */
    //@{
    void SetOneSiteForceON( )
    {   this->bOneSiteForce = true;     }
    void SetOneSiteForceOFF( )
    {   this->bOneSiteForce = false;    }
    //@}

    /** Set/Get dihedral angle indices */
    //@{
    void SetDihedralIndices( int *dihedralIndices )
    {   this->dihedralIndices = dihedralIndices;}
    int *GetDihedralIndices( ) const
    {   return this->dihedralIndices;           }
    //@}

    /** Set periodicity on the embedded data */
    void SetPeriodicity( double *xi_min, double *xi_max )
    {
        this->bPeriodic = true;
        this->xi_min = xi_min;
        this->xi_max = xi_max;
    }

    /** Set Usage of Procrustes Analysis ON or OFF (default is [OFF]) */
    //@{
    void SetProcrustesAnalysisON( )
    {   this->bUseProcrustesAnalysis = true;    }
    void SetProcrustesAnalysisOFF( )
    {   this->bUseProcrustesAnalysis = false;   }
    //@}


    /** Set active atoms for aligning with procrustes analysis */
    void SetActiveReferenceAtoms( int nActiveAtoms, int *activeAtomsIndex )
    {
        this->nActiveAtoms     = nActiveAtoms;
        this->activeAtomsIndex = activeAtomsIndex;
    }


    /** Set Usage of Contact Map ON or OFF (default is [OFF]) */
    //@{
    void SetContactMapON( )
    {   this->bUseContactMap = true;            }
    void SetContactMapOFF( )
    {   this->bUseContactMap = false;           }
    //@}

    /** Set Usage of Contact Map ON or OFF (default is [OFF]) */
    //@{
    void SetAlignedCPPON( )
    {   this->bUseAlignedCPP = true;            }
    void SetAlignedCPPOFF( )
    {   this->bUseAlignedCPP = false;           }
    //@}

    /** Set Usage of Contact Map ON or OFF (default is [ON]) */
    //@{
    void SetCurvatureCorrectionON( )
    {   this->bUseCurvatureCorrection = true;   }
    void SetCurvatureCorrectionOFF( )
    {   this->bUseCurvatureCorrection = false;  }
    //@}

    /** Define the required parameter for adding weights to some dimensions in calculating
      * closest point projection
      */
    void SetWeightedCPP( int nWeights, int *weightIndex, double *weightValue)
    {
        this->bUseWeightedCPP = true;
        this->nWeights        = nWeights;
        this->weightIndex     = weightIndex;
        this->weightValue     = weightValue;
    }


    /** Set Usage of Corral ON or OFF (default is [OFF]) */
    //@{
    void SetCorralON( )
    {   this->bUseCorral = true;                }
    void SetCorralOFF( )
    {   this->bUseCorral = false;               }
    //@}


    /** Set/Get Corral Scale */
    //@{
    void   SetCorralScale( double scale )
    {   this->corralScale = scale;              }
    double GetCorralScale( ) const
    {   return this->corralScale;               }
    //@}

    /** Set numebr of bins for corral
    */
    void SetNumberOfCorralBins( int nBins )
    {   this->nCorralBins = nBins;              }
    int  GetNumberOfCorralBins( ) const
    {   return this->nCorralBins;               }


    /** Set a way that the corral markers get loaded from a file */
    void SetLoadCorralMarkers( const char *markersFilename )
    {
        this->bLoadCorralMarkers = true;
        this->markersFilename    = markersFilename;
    }


    void SetContactFunction( char *contactFunctionName  )
    {   this->contactFunction = contactFunctionName; }

    /** Set Contact Map Parameters
    * kappa  : cut-off length (Angstrom)
    * epsilon: window size
    */
    void SetContactMapParameters( double kappa, double epsilon )
    {
        this->kappa   = kappa;
        this->epsilon = epsilon;
    }

    /** Set Contact Map Parameters with special consideration for end-to-end distance
    * kappa   : cut-off length (Angstrom)
    * epsilon : window size
    * kappaE  : cut-off length for end-to-end contact
    * epsilonE: window size    for end-to-end contact
    */
    void SetContactMapParameters( double kappa, double epsilon, double kappaE, double epsilonE )
    {
        this->bUseEndToEndContact = true;
        this->kappa               = kappa;
        this->epsilon             = epsilon;
        this->kappaE              = kappaE;
        this->epsilonE            = epsilonE;
    }

    /** Set a factor which is multiplied by end-to-end contact force
    * factorE : force factor for end-to-end contact
    */
    void SetEndToEndFactor( double factorE )
    {
        this->factorE  = factorE;
    }

    /** This function sets the gamma parameter used to computed the thermalization
      * parameter for each node, \f$ \beta_a = \gamma h_a^2 \f$
      *
      * \param gamma_ns  nodal spacing      [1.0, 4.0], default value 0.0
      * \param gamma_lme local max-entropy  [0.4, 2.0], default value 0.0
      *
      * Overloaded function for setting the gamma parameters.
      * \param gamma[] [gamma_ns, gamma_lme]
      */
    //@{
    void SetGammaLME( double gamma_lme )
    {   this->gamma_lme = gamma_lme;            }
    void SetGammaNS( double gamma_ns )
    {   this->gamma_ns = gamma_ns;              }
    void SetGammaPU( double gamma_pu )
    {   this->gamma_pu = gamma_pu;              }
    double GetGammaNS() const
    {   return this->gamma_ns;                  }
    double GetGammaLME() const
    {   return this->gamma_lme;                 }
    double GetGammaPU() const
    {   return this->gamma_pu;                  }
    //@}


    /** The node points spacing is fixed to a constant value */
    //@{
    void SetNodeSpacing( double nodalSpacing )
    {   this->nodalSpacing = nodalSpacing;      }
    double GetNodeSpacing( ) const
    {   return this->nodalSpacing;              }
    //@}


    /** This function defines the number of nearest neighbors requested for computing
      * the spacing between node points. Default value NNR=12
      */
    //@{
    void SetNumberOfNeighbours(int NNR)
    {   this->nNearestNeighbors = NNR;          }
    int GetNumberOfNeighbors( ) const
    {   return this->nNearestNeighbors;         }
    //@}


    /** The high-dimensional node points filename */
    //@{
    void SetNodePointsFilename( const char* xNodeFilename )
    {   this->xNodesFilename = xNodeFilename;   }
    const char *GetNodePointsFilename( ) const
    {   return this->xNodesFilename;            }
    //@}


    /** The low-dimensional node points filename */
    //@{
    void SetEmbeddedPointsFilename( const char* xiNodeFilename )
    {   this->xiNodesFilename = xiNodeFilename; }
    const char *GetEmbeddedPointsFilename( ) const
    {   return this->xiNodesFilename;           }
    //@}


    /** The high-dimensional reference node point filename for procrustes analysis alignment*/
    //@{
    void SetReferencePointFilename( const char* xReferenceFilename )
    {   this->xReferenceFilename = xReferenceFilename; }
    const char *GetReferencePointFilename( ) const
    {   return this->xReferenceFilename;               }
    //@}


    /** The low-dimensional seed points filename */
    //@{
    void SetSeedPointsFilename( const char* xiSeedsFilename )
    {   this->xiSeedsFilename = xiSeedsFilename; }
    const char *GetSeedPointsFilename( ) const
    {   return this->xiSeedsFilename;           }
    //@}

    /** Get distance between the new points and its projection
      */
    double GetDistance( )
    {   return this->distance;                  }


    /** Compute forces
      * \return error identifier */
    int GetForces( double *x_new, double **x_projected, double **xi_projected, double **dxi_dx);

    /** Compute forces and corral gradints
      * \return error identifier */
    int GetForces( double *x_new, double **x_projected, double **xi_projected, double **dxi_dx,
                   double **dC_dxi );


    /** Compute the closest point projection without calculating the forces
      * \return error identifier */
    //@{
    void GetProjection( double *x_new );
    int  GetProjection( double *x_new, double **x_projected, double **xi_projected );
    //@}

    /** Apply Curvature Correction
   * \return error identifier */
    int ComputeJacobian( double *CorrectedJacobian );




    /** Bring the algorithm's outputs up-to-date. */
    void Update();


protected:
    bool        bVerbose;
    bool        bPeriodic;
    bool        bOneSiteForce;
    bool        bUseCorral;
    bool        bUseProcrustesAnalysis;
    bool        bUseContactMap;
    bool        bUseAlignedCPP;
    bool        bUseWeightedCPP;
    bool        bUseEndToEndContact;
    bool        bLoadCorralMarkers;
    bool        bUseCurvatureCorrection;
    bool        bUseDihedralAngles;
    const char *markersFilename;
    const char *xNodesFilename;
    const char *xiNodesFilename;
    const char *xiSeedsFilename;
    int         ierr;
    int         rDim;
    int         nCorralBins;
    int        *dihedralIndices;
    double      corralScale;
    double      nodalSpacing;
    double     *x_ref;

    double     *xi_min;
    double     *xi_max;
    
    // Localization variables
    int         nNearestNeighbors;
    double      gamma_ns;
    double      gamma_lme;
    double      gamma_pu;

    // General Alignment valiables
    int         nAtoms;
    double     *x_aligned;
    double     *dxi_dx_dealigned;

    // Procrustes Analysis variables
    const char *xReferenceFilename;
    int         nActiveAtoms;
    int        *activeAtomsIndex;

    // Contact Map valriables
    char       *contactFunction;
    double      kappa;
    double      epsilon;
    double      kappaE;
    double      epsilonE;
    double      factorE;

    // Weighted closest point projection
    int         nWeights;
    int        *weightIndex;
    double     *weightValue;

    //Curvature Correction Workspace
    double      distance;
    double     *xi_projected;
    double     *x_projected;
    double     *dx_projected;
    double     *ddx_projected;
    double     *Jacobian;
    double     *CorrectedJacobian;
    double     *Gtemp;
    double     *Ginv;
    double     *u ;
    double     *ddu;
    double     *du_xi;
    double     *correction;
    double     *correctionInv;


    smPointSet::Pointer                 nodes1;
    smPointSet::Pointer                 nodes2;
    smPointSet::Pointer                 nodes3;
    smPointSet::Pointer                 nodes4;
    smPointSetIO::Pointer           reader1;
    smPointSetIO::Pointer           reader2;
    smPointSetIO::Pointer           reader3;
    smPointSetIO::Pointer           reader4;
    smContactMap::Pointer               ContactMap;
    smProcrustesAnalysis::Pointer       ProcrustesAnalysis;
    smClosestPointProjectionLS::Pointer cpp;
    smContinuousCorral::Pointer         corral;
    smDihedralAngles::Pointer           dha;

    /** Constructor */
    smAcceleratedMD(void);

    /** Destructor */
    ~smAcceleratedMD();

    /** Check if parameters and node/sample point coordinates are propely set. */
    int CheckSettings();

    /** Method to delete the preallocated member variables */
    void Clear();

    /**  In this function are initialized the objects and variables needed to compute
         * the parametrical samples coordinates in the embedded low-dimensional space.
         */
    int Initialize( );

private:

    smAcceleratedMD(const Self&);	//purposely not implemented
    void operator=(const Self&);	//purposely not implemented
};

#endif //__smAcceleratedMD_h

