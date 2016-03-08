/*=========================================================================

Module    : Solid Mechanics
File      : smAcceleratedMD.cpp
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
#include "smAcceleratedMD.h"

smAcceleratedMD::smAcceleratedMD( )
{
d    //____________________________________________________________________________________
    // Inputs
    this->bPeriodic              = false;
    this->bUseContactMap         = false;
    this->bUseProcrustesAnalysis = false;
    this->bUseCorral             = false;
    this->bUseAlignedCPP         = false;
    this->bUseWeightedCPP        = false;
    this->bUseEndToEndContact    = false;
    this->bLoadCorralMarkers     = false;
    this->bUseCurvatureCorrection = true;
    this->bOneSiteForce          = false;
    this->bUseDihedralAngles     = false;
    this->corralScale            = 1.0;
    this->tDim                   = 0;
    this->nDim                   = 0;
    this->nPts                   = 0;
    this->nodalSpacing           = 0.0;
    this->xNodesFilename         = NULL;
    this->xiNodesFilename        = NULL;
    this->xReferenceFilename     = NULL;
    this->x_nodes                = NULL;
    this->xi_nodes               = NULL;

    this->xi_min                 = NULL;
    this->xi_max                 = NULL;
    
    // gamma
    this->gamma_pu               = 1.0;
    this->gamma_ns               = 1.0;
    this->gamma_lme              = 1.0;

    // Procrustes Analysis
    this->nActiveAtoms           = 0;
    this->activeAtomsIndex       = NULL;
    this->x_ref                  = NULL;

    // Contact Map
    this->kappa                  = 0.0;
    this->epsilon                = 0.0;
    this->kappaE                 = 0.0;
    this->epsilonE               = 0.0;
    this->factorE                = 1.0;

    // Weighted closest point projection
    this->nWeights               = 0;
    this->weightIndex            = NULL;
    this->weightValue            = NULL;

    // Workspace variables
    this->x_aligned              = NULL;
    this->dxi_dx_dealigned       = NULL;
    this->xi_projected           = NULL;
    this->x_projected            = NULL;
    this->dx_projected           = NULL;
    this->ddx_projected          = NULL;
    this->u                      = NULL;
    this->ddu                    = NULL;
    this->du_xi                  = NULL;
    this->Gtemp                  = NULL;
    this->Ginv                   = NULL;
    this->correction             = NULL;
    this->correctionInv          = NULL;
    this->Jacobian               = NULL;
    this->CorrectedJacobian      = NULL;


    // Class Pointers
    this->nodes1                 = false;
    this->nodes2                 = false;
    this->nodes3                 = false;
    this->reader1                = false;
    this->reader2                = false;
    this->reader3                = false;
    this->ProcrustesAnalysis     = false;
    this->ContactMap             = false;
    this->cpp                    = false;
    this->dha                    = false;

    //____________________________________________________________________________________
    // Outputs
    this->ierr                   = 0;
}




smAcceleratedMD::~smAcceleratedMD( )
{
    this->Clear();
}




//cleaning preallocated member variables
void smAcceleratedMD::Clear( )
{
    if ( this->x_aligned         != NULL ) delete[] this->x_aligned         ;
    if ( this->dxi_dx_dealigned  != NULL ) delete[] this->dxi_dx_dealigned  ;
    if ( this->u                 != NULL ) delete[] this->u                 ;
    if ( this->ddu               != NULL ) delete[] this->ddu               ;
    if ( this->du_xi             != NULL ) delete[] this->du_xi             ;
    if ( this->Gtemp             != NULL ) delete[] this->Gtemp             ;
    if ( this->Ginv              != NULL ) delete[] this->Ginv              ;
    if ( this->correction        != NULL ) delete[] this->correction        ;
    if ( this->correctionInv     != NULL ) delete[] this->correctionInv     ;
    if ( this->Jacobian          != NULL ) delete[] this->Jacobian          ;
    if ( this->CorrectedJacobian != NULL ) delete[] this->CorrectedJacobian ;

    this->xi_nodes          = NULL;
    this->x_nodes           = NULL;
    this->x_ref             = NULL;
    this->x_aligned         = NULL;
    this->dxi_dx_dealigned  = NULL;
    this->xi_projected      = NULL;
    this->x_projected       = NULL;
    this->dx_projected      = NULL;
    this->ddx_projected     = NULL;
    this->u                 = NULL;
    this->ddu               = NULL;
    this->du_xi             = NULL;
    this->Gtemp             = NULL;
    this->Ginv              = NULL;
    this->correction        = NULL;
    this->correctionInv     = NULL;
    this->Jacobian          = NULL;
    this->CorrectedJacobian = NULL;
    
    this->ierr              = 0;

}




void smAcceleratedMD::Update( )
{
    if (this->CheckSettings())
    {
        printf("ERROR:smAcceleratedMD input settings are wrong\n");
        this->ierr = 1;
        return;
    }

    // member variables and objects are initialized
    this->ierr = this->Initialize();
    if (this->ierr != 0)
    {
        printf("ERROR:smAcceleratedMD Initialization failed with Error = %d\n",this->ierr);
    }

    return;
}



//check the setting before to compute
int smAcceleratedMD::CheckSettings()
{
    if(!this->bUseDihedralAngles)
    {
        if (this->bUseProcrustesAnalysis && this->bUseContactMap)
        {
            printf("ERROR::smAcceleratedMD::CheckSettings :\n");
            printf("Both contact map and procrustes analysis have been chosen!\n");
            return 1;
        }
        else if (!this->bUseProcrustesAnalysis && !this->bUseContactMap)
        {
            printf("ERROR::smAcceleratedMD::CheckSettings :\n");
            printf("Neither contact map nor procrustes analysis has been chosen!\n");
            return 1;
        }

        if (this->xNodesFilename == NULL)
        {
            printf("ERROR::smAcceleratedMD::CheckSettings :\n");
            printf("High-dimensional node points filename has not been set!\n");
            return 1;
        }
        if (this->xiNodesFilename == NULL)
        {
            printf("ERROR::smAcceleratedMD::CheckSettings :\n");
            printf("Embedded node points filename has not been set!\n");
            return 1;
        }
        if (this->bUseProcrustesAnalysis)
        {
            if (this->xReferenceFilename == NULL)
            {
                printf("ERROR::smAcceleratedMD::CheckSettings :\n");
                printf("Reference filename for procrustes analysis has not been set!\n");
                return 1;
            }
        }
        if (this->nActiveAtoms > 0)
        {
            if (this->activeAtomsIndex == NULL)
            {
                printf("ERROR::smAcceleratedMD::CheckSettings :\n");
                printf("Procrustes analysis active atoms indices vector has not been set!\n");
                return 1;
            }
        }
    }
    if (this->bLoadCorralMarkers)
    {
        if (this->markersFilename == NULL)
        {
            printf("ERROR::smAcceleratedMD::CheckSettings :\n");
            printf("Markers filename for generating corral has not been set!\n");
            return 1;
        }
    }
    //cleaning preallocated member variables
    this->Clear();

    return 0;
}




int smAcceleratedMD::Initialize()
{
    if(this->bUseDihedralAngles)
    {
        this->SetEmbeddedDimension(2);
        this->SetNumberOfNodePoints(1); // will be overwritten
        this->dha = smDihedralAngles::New();
        if (this->bOneSiteForce)
        {
            this->dha->SetOneSiteForceON();
            cout<<"One site force is ON!"<<endl;
        }
        this->dha->SetNumberOfAtoms( this->nAtoms );
        this->dha->SetDihedralAngles( this->tDim, this->dihedralIndices );
        this->dha->Update();
        cout<<"tDim   ===== "<<this->tDim<<endl;
        cout<<"nAtoms ===== "<<this->nAtoms<<endl;
    }
    else
    {
        //NOTE [1]
        /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        *  Reading input data sets
        */
        //____________________________________________________________________________________
        // [1.1] Load node points in HIGH dimension
        this->reader1 = smPointSetIO::New();
        this->reader1->SetPointsFileName( this->xNodesFilename );
        this->reader1->Update();
        this->ierr    = this->reader1->GetError();
        if (this->ierr)
        {
            printf("ERROR::smAcceleratedMD::Initialize :\n");
            printf("Reader of the point coordinates file have error = %d\n", this->ierr);
            return 1;
        }
        this->nodes1  = this->reader1->GetOutput();
        this->SetNodePoints        ( this->nodes1->GetPoints()         );
        this->SetDimension         ( this->nodes1->GetDimension()      );
        this->SetNumberOfNodePoints( this->nodes1->GetNumberOfPoints() );
        printf("[%s] is loaded. \n", this->xNodesFilename);
        //____________________________________________________________________________________
        // [1.2] Load node points in LOW dimension
        this->reader2 = smPointSetIO::New();
        this->reader2->SetPointsFileName( this->xiNodesFilename );
        this->reader2->Update();
        this->ierr    = this->reader2->GetError();
        if (this->ierr)
        {
            printf("ERROR::smAcceleratedMD::Initialize :\n");
            printf("Reader of the point coordinates file have error = %d\n", this->ierr);
            return 1;
        }
        this->nodes2    = this->reader2->GetOutput();
        this->SetEmbeddedNodePoints( (double *) this->nodes2->GetPoints() );
        this->SetEmbeddedDimension ( this->nodes2->GetDimension() );
        if (this->tDim > 4)
        {
            printf("ERROR::smAcceleratedMD::Initialize :\n");
            printf("Embedded dimension more than 3 is not supported yet!\n");
            return 1;
        }
        if (this->nodes2->GetNumberOfPoints() != this->nPts)
        {
            if(!this->nodalSpacing) {
                if (this->nodes2->GetNumberOfPoints()-1 != this->nPts) {
                    printf("ERROR::smAcceleratedMD::Initialize :\n");
                    printf("Number of sample points in Low and High dimension doesn't match!\n");
                    return 1;
                }
            }
        }
        printf("[%s] is loaded. \n", this->xiNodesFilename);
        //____________________________________________________________________________________
        // [1.3] Load reference node point in HIGH dimension for procrustes analysis
        if (bUseProcrustesAnalysis)
        {
            this->reader3 = smPointSetIO::New();
            this->reader3->SetPointsFileName( this->xReferenceFilename );
            this->reader3->Update();
            this->ierr = this->reader3->GetError();
            if (this->ierr)
            {
                printf("ERROR::smAcceleratedMD::Initialize :\n");
                printf("Reader of the point coordinates file have error = %d\n", this->ierr);
                return 1;
            }
            this->nodes3 = this->reader3->GetOutput();
            this->x_ref  = (double *)this->nodes3->GetPoints();
            if (this->nDim != this->nodes3->GetDimension())
            {
                printf("ERROR::smAcceleratedMD::Initialize :\n");
                printf("Dimension of reference point and node points doesn't match!\n");
                return 1;
            }
            if (this->nodes3->GetNumberOfPoints() != 1)
            {
                printf("ERROR::smAcceleratedMD::Initialize :\n");
                printf("Number of reference points has to be ONE!\n");
                return 1;
            }
            printf("[%s] is loaded. \n", this->xReferenceFilename);
        }
        //____________________________________________________________________________________
        // [1.4] Load seed points in LOW dimension
        this->reader4 = smPointSetIO::New();
        this->reader4->SetPointsFileName( this->xiSeedsFilename );
        this->reader4->Update();
        this->ierr    = this->reader4->GetError();
        if (this->ierr)
        {
            printf("ERROR::smAcceleratedMD::Initialize :\n");
            printf("Reader of the point coordinates file have error = %d\n", this->ierr);
            return 1;
        }
        this->nodes4    = this->reader4->GetOutput();
        printf("[%s] is loaded. \n", this->xiSeedsFilename);
        //____________________________________________________________________________________
        // [1.5] Print Loaded Data Information
        if (this->bVerbose)
        {
            printf("\nInputs:\n" );
            printf("\tNumber of node points = %d\n", this->nPts           );
            printf("\tEmbedded  dimension   = %d\n", this->tDim           );
            printf("\tHigh      dimension   = %d\n", this->nDim           );
            printf("\tx_nodes  filename     : %s\n", this->xNodesFilename );
            printf("\txi_nodes filename     : %s\n", this->xiNodesFilename);
            printf("\txi_seeds filename     : %s\n", this->xiSeedsFilename);
            if (bUseProcrustesAnalysis)
                printf("\tx_ref    filename     : %s\n", this->xReferenceFilename);
        }
    }
    //NOTE [2]
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    *  Initialize Corral
    */
    if(this->bUseCorral)
    {
        if (this->bVerbose)
            printf("THE CORRAL IS ON   Corral height = %6.4f\n", this->corralScale);
        this->corral = smContinuousCorral::New();
        this->corral->SetEmbeddedDimension(this->tDim);
        this->corral->SetNumberOfNodePoints(this->nPts);
        this->corral->SetPartitionOfUnityOptions(this->gamma_pu);
        this->corral->SetCorralScale(this->corralScale);
        if (this->bLoadCorralMarkers)
        {
            this->corral->SetLoadMarkersON();
            this->corral->SetMarkersFilename(this->markersFilename);
        }
        else
        {
            this->corral->SetEmbeddedNodePoints(this->xi_nodes);
            this->corral->SetNumberOfBins(this->nCorralBins);
        }
        this->corral->Update();
    }
    if(!this->bUseDihedralAngles)
    {
        //NOTE [3]
        /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        *  Initialize Closest Point Projection
        */
        double Tol0  = 1e-10;    //Tolerence of LME
        double TolNR = 1e-14;    //Tolerence of Newton-Raphson

        this->cpp = smClosestPointProjectionLS::New();
        if (this->bVerbose)
            this->cpp->SetVerboseON();
        else
            this->cpp->SetVerboseOFF();
        this->cpp->SetGamma( this->gamma_ns, this->gamma_lme );
        this->cpp->SetTargetZero( Tol0 );
        this->cpp->SetNewtonTolerance( TolNR );
        if(!this->nodalSpacing) {
            this->nodalSpacing  = this->xi_nodes[this->tDim*this->nPts];
            cout<<"NodalSpacing = "<<this->nodalSpacing<<endl;
        }
        this->cpp->SetNodeSpacing( this->nodalSpacing );
        this->cpp->SetNumberOfNeighbours( this->nNearestNeighbors );
        this->cpp->SetDimension( this->nDim );
        this->cpp->SetNodePoints( this->x_nodes );
        this->cpp->SetEmbeddedDimension( this->tDim );
        this->cpp->SetEmbeddedNodePoints( this->xi_nodes );
        this->cpp->SetNumberOfNodePoints( this->nPts );
        this->cpp->SetUseLastPointON(); // use last converged point as the new seed
        // Define seed points
        this->cpp->SetNumberOfSeedPoints( this->nodes4->GetNumberOfPoints() );
        this->cpp->SetSeedPoints( (double *)this->nodes4->GetPoints() );
        if(!this->bUseCorral) {
            this->cpp->SetPeriodicity(xi_min,xi_max);
        }
        // _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
//         if(this->bUseCorral)
//         {
//             // Define corral markers as seed points
//             this->cpp->SetNumberOfSeedPoints( this->corral->GetNumberOfInsideMarkerPoints() );
//             this->cpp->SetSeedPoints( this->corral->GetMarkerPoints() );
//         }
        // _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
        if (this->bUseAlignedCPP && this->bUseProcrustesAnalysis)
            this->cpp->SetAlignmentON();
        if (this->bUseWeightedCPP)
            this->cpp->SetWeights( this->nWeights, this->weightIndex, this->weightValue );
        this->cpp->Update();

        if (this->cpp->GetError( ))
        {
            printf("ERROR::smAcceleratedMD::Initialize :\n");
            printf("Initialization Closest Point Projection has error = %d\n",
                   this->cpp->GetError( ));
            return 1;
        }
        //NOTE [4]
        /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        *  Initialize Procrustes Analysis or Contact Map
        */
        //____________________________________________________________________________________
        // [3.1] Initializing Procrustes Analysis
        if (this->bUseProcrustesAnalysis)
        {
            this->nAtoms = this->nDim/3;
            this->rDim   = this->nDim;

            this->ProcrustesAnalysis = smProcrustesAnalysis::New();
            this->ProcrustesAnalysis->SetDimension( 3 );
            this->ProcrustesAnalysis->SetNumberOfPoints( this->nAtoms );
            if (this->nActiveAtoms > 0) {
                cout<<"active atoms = "<<this->nActiveAtoms<<endl;
                smIO::PrintVector(this->nActiveAtoms,this->activeAtomsIndex);
                this->ProcrustesAnalysis->SetReferencePoints( this->nActiveAtoms, this->activeAtomsIndex, this->x_ref );
            }
            else
                this->ProcrustesAnalysis->SetReferencePoints( this->x_ref );
            this->ProcrustesAnalysis->SetOptions( false, false );
            this->ProcrustesAnalysis->Update();
            this->ierr = this->ProcrustesAnalysis->GetError();
            if (this->ierr)
            {
                printf("ERROR::smAcceleratedMD::Initialize :\n");
                printf("Procrustes-analysis initialization failed with error = %d\n", this->ierr);
                return 2;
            }
        }
        //____________________________________________________________________________________
        // [4.1] Initializing Contact Map
        else if (this->bUseContactMap)
        {
            this->nAtoms = int(0.5*(sqrt(this->nDim*8.0+1.0)+1.0));
            this->rDim   = 3*this->nAtoms;

            this->ContactMap = smContactMap::New();
            this->ContactMap->SetSpatialDimension( 3 );
            this->ContactMap->SetNumberOfAtoms( this->nAtoms );
            this->ContactMap->SetContactFunction(this->contactFunction);
            cout<<"Contact Function is " << this->contactFunction<<endl;
            if (this->bUseEndToEndContact)
            {
                this->ContactMap->SetContactMapParameters( this->kappa, this->epsilon, this->kappaE, this->epsilonE);
                this->ContactMap->SetEndToEndFactor( this->factorE );
                cout<<"Factor = "<<factorE<<endl;
            }
            else
                this->ContactMap->SetContactMapParameters( this->kappa, this->epsilon );
            this->ContactMap->Update();
            this->ierr = this->ContactMap->GetError();
            if (this->ierr)
            {
                printf("ERROR::smAcceleratedMD::Initialize :\n");
                printf("Contact-map initialization failed with error = %d\n", this->ierr);
                return 2;
            }
        }
        //____________________________________________________________________________________
        // [4.2] Allocating memory for alignment variables
        this->x_aligned        = new double [this->nDim];
        this->dxi_dx_dealigned = new double [this->rDim*this->tDim];

        //NOTE [5]
        /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        *  Initialize the Curvature Correction Workspace
        */
        //____________________________________________________________________________________
        // [5.1] Allocating memory for curvature correction variables
        this->u                 = new double [this->nDim];
        this->ddu               = new double [this->tDim*this->tDim];
        this->du_xi             = new double [this->nDim*this->tDim];
        this->Gtemp             = new double [this->tDim*this->tDim];
        this->Ginv              = new double [this->tDim*this->tDim];
        this->correction        = new double [this->tDim*this->tDim];
        this->correctionInv     = new double [this->tDim*this->tDim];
        this->Jacobian          = new double [this->nDim*this->tDim];
        this->CorrectedJacobian = new double [this->nDim*this->tDim];
    }
    return ierr;
}




int smAcceleratedMD::GetForces( double *x_new, double **x_projected, double **xi_projected, double **dxi_dx,
                                double **dC_dxi)
{
    if(this->bUseDihedralAngles) {
        this->dha->ComputeDihedralAngles( x_new );
        *x_projected  = x_new;
        *xi_projected = this->dha->GetDihedralAnglesValues();
        *dxi_dx       = this->dha->GetDihedralAnglesGradients();
        *dC_dxi       = this->dha->GetDihedralAnglesInvGradients();
    } else {
        //____________________________________________________________________________________
        // [1] Find the closest point projection of the new configuration
        this->GetProjection( x_new );
//        smIO::WriteMatrix(1,this->nDim,this->x_projected ,"xProjected.txt" );
//        smIO::WriteMatrix(1,this->tDim,this->xi_projected,"xiProjected.txt");
        *x_projected  = this->x_projected;
        *xi_projected = this->xi_projected;

        //____________________________________________________________________________________
        // [2] Apply Curvature Correction
        this->ComputeJacobian( this->CorrectedJacobian );

        //____________________________________________________________________________________
        // [3] Apply Dealignment
        if (this->bUseContactMap) {
            this->ContactMap->DeContactMap( this->tDim, x_new, this->x_aligned, this->CorrectedJacobian, this->dxi_dx_dealigned );
        }
        else if (bUseProcrustesAnalysis) {
            this->ProcrustesAnalysis->Dealignment( this->tDim, x_new, this->CorrectedJacobian, this->dxi_dx_dealigned );
        }
        else {
            printf("ERROR::smAcceleratedMD::GetForces :\n");
            printf("No dealignment method has been set.\n");
            return 1;
        }

        *dxi_dx = this->dxi_dx_dealigned;


        if (this->bUseCorral) {
            if(this->bUseDihedralAngles) {
                double angles[2];
                if ( **xi_projected < 0 )
                    angles[0] = **xi_projected+2*PI;
                else
                    angles[0] = **xi_projected;
                if ( *(*xi_projected+1) < 0 )
                    angles[1] = *(*xi_projected+1)+2*PI;
                else
                    angles[1] = *(*xi_projected+1);
                cout<<"corralAngles = [ "<<angles[0]<<" , "<<angles[1]<<" ]"<<endl;
                this->corral->ComputeCorral( angles );
            } else {
                this->corral->ComputeCorral( *xi_projected );
            }
            *dC_dxi = this->corral->GetCorralGradient();
            double C_xi = this->corral->GetCorralValue();
            smIO::WriteMatrix( 1, 1, &C_xi, "corral_val.txt" );
        }
        else
            *dC_dxi = NULL;
    }
    return this->ierr;
}




int smAcceleratedMD::GetForces( double *x_new, double **x_projected, double **xi_projected, double **dxi_dx)
{
    //____________________________________________________________________________________
    // [1] Find the closest point projection of the new configuration
    this->GetProjection( x_new );
    *x_projected  = this->x_projected ;
    *xi_projected = this->xi_projected;

    if ( this->ierr == 0 )
    {
        //____________________________________________________________________________________
        // [2] Apply Curvature Correction
        this->ComputeJacobian( this->CorrectedJacobian );

        //____________________________________________________________________________________
        // [3] Apply Dealignment
        if (this->bUseContactMap)
        {
            this->ContactMap->DeContactMap( this->tDim, x_new, this->x_aligned, this->CorrectedJacobian, this->dxi_dx_dealigned );
        }
        else if (this->bUseProcrustesAnalysis)
        {
            this->ProcrustesAnalysis->Dealignment( this->tDim, x_new, this->CorrectedJacobian, this->dxi_dx_dealigned );
        }
        else
        {
            printf("ERROR::smAcceleratedMD::GetForces :\n");
            printf("No dealignment method has been set.\n");
            return 1;
        }

        *dxi_dx = this->dxi_dx_dealigned;
    }
    return this->ierr;
}




int smAcceleratedMD::GetProjection( double *x_new, double **x_projected, double **xi_projected )
{
    GetProjection( x_new );
    *x_projected  = this->x_projected ;
    *xi_projected = this->xi_projected;
    return this->ierr;
}


void smAcceleratedMD::GetProjection( double *x_new )
{
    //____________________________________________________________________________________
    // [1] Find the closest point projection of the new configuration
    if (this->bUseContactMap)
    {
        this->ContactMap->ComputeContactMap( x_new, this->x_aligned );
    }
    else if (bUseProcrustesAnalysis)
    {
        if (this->bUseAlignedCPP)
        {
            if (this->nActiveAtoms>0)
                this->ProcrustesAnalysis->ComputeRigidAlignment3D( x_new, this->x_aligned );
            else
                this->ProcrustesAnalysis->ComputeRigidAlignmentAll3D( x_new, this->x_aligned );
        }
        else
        {
            if (this->nActiveAtoms>0)
                this->ProcrustesAnalysis->ComputeRigidAlignment3DandDerivatives( x_new, this->x_aligned );
            else
                this->ProcrustesAnalysis->ComputeRigidAlignmentAll3DandDerivatives( x_new, this->x_aligned );
        }
//         smIO::WriteMatrix(this->nAtoms, 3, x_new, "x_new.txt");
//         smIO::WriteMatrix(this->nAtoms, 3, x_aligned, "x_aligned.txt");
    }
    else
    {
        printf("ERROR::smAcceleratedMD::GetForces :\n");
        printf("No alignment method has been set.\n");
        exit(1);
    }

    this->cpp->ComputeClosestPointProjection( this->x_aligned );
    this->x_projected  = (double *)this->cpp->GetProjectedSample();
    this->xi_projected = (double *)this->cpp->GetSamplesEmbedding();

    if (this->bUseAlignedCPP)
    {
        if (this->nActiveAtoms>0)
        {
            this->ProcrustesAnalysis->SetReferencePoints( this->nActiveAtoms, this->activeAtomsIndex, (double *)this->cpp->GetProjectedSample() );
            this->ProcrustesAnalysis->ComputeRigidAlignment3DandDerivatives( x_new, this->x_aligned );
        }
        else
        {
            this->ProcrustesAnalysis->SetReferencePoints( (double *)this->cpp->GetProjectedSample() );
            this->ProcrustesAnalysis->Update();
            this->ProcrustesAnalysis->ComputeRigidAlignmentAll3DandDerivatives( x_new, this->x_aligned );
        }
    }
    this->ierr = this->cpp->GetError();
}



int smAcceleratedMD::ComputeJacobian( double *CorrectedJacobian )
{
    int    i; //Counter  over nDim
    int    j; //Counters over tDim
    //____________________________________________________________________________________
    // [1] Obtain closest point projection information
    this->x_projected   = (double *)this->cpp->GetProjectedSample();
    this->dx_projected  = (double *)this->cpp->GetSampleGradients();
    this->ddx_projected = (double *)this->cpp->GetSampleHessians();

    //____________________________________________________________________________________
    // [2] Compute G inverse
    smMath::Tensor(this->tDim, this->nDim, this->dx_projected, this->Gtemp);

    if ( this->tDim == 1 )
        this->Ginv[0] = 1.0/this->Gtemp[0];
    else if ( this->tDim == 2 )
        smMath::Invert2x2( this->Gtemp, this->Ginv );
    else if ( this->tDim == 3 )
        smMath::Invert3x3( this->Gtemp, this->Ginv );

    //____________________________________________________________________________________
    // [3] Compute Jacobian (Psuedo Inverse Mapping)
    smMath::MatrixMultiplication( this->tDim, this->tDim, this->nDim,
                                  this->Ginv, this->dx_projected, this->Jacobian );

    //____________________________________________________________________________________
    // [4] Compute Curvature Correction
    //_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
    // [4.1] Compute distance vector between the aligned configuration and its projection
    for ( i = 0; i < this->nDim; i++ )
        this->u[i] = this->x_aligned[i] - this->x_projected[i];
    this->distance = smMath::Norm(this->nDim, this->u);
    if (this->distance < 1e-10 || !this->bUseCurvatureCorrection)
    {
        smArray::Copy ( this->tDim*this->nDim, this->Jacobian, CorrectedJacobian );
        return 0;
    }
    for ( i = 0; i < this->nDim; i++ )
        this->u[i] /= this->distance;
    //_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
    // [4.2] Compute some du/dxi
    smMath::MatrixMultiplication( this->tDim*this->tDim, this->nDim, 1,
                                  this->ddx_projected, this->u, this->Gtemp );
    smMath::MatrixMultiplication( this->tDim, this->tDim, this->tDim,
                                  this->Ginv, this->Gtemp, this->ddu );
    smMath::MatrixTransposition (this->tDim, this->nDim, this->dx_projected );
    smMath::MatrixMultiplication( this->nDim, this->tDim, this->tDim,
                                  this->dx_projected, this->ddu, this->du_xi );
    //_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
    // [4.3] Compute the correction value matrix
    smMath::MatrixMultiplication( this->tDim, this->nDim, this->tDim,
                                  this->Jacobian, this->du_xi, this->correction );
    for ( i = 0; i < this->tDim*this->tDim; i++ )
        this->correction[i] *= -this->distance;

    for ( j = 0; j < tDim; j++ ) {
        this->correction[ j*tDim + j ] += 1.0;
    }
    if ( this->tDim == 1 )
        this->correctionInv[0] = 1.0/this->correction[0];
    else if ( this->tDim == 2 )
        smMath::Invert2x2( this->correction, this->correctionInv );
    else if ( this->tDim == 3 )
        smMath::Invert3x3( this->correction, this->correctionInv );
    //_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
    // [4.4] Applying the computed correction
    smMath::MatrixMultiplication( this->tDim, this->tDim, this->nDim,
                                  this->correctionInv, this->Jacobian, this->CorrectedJacobian );

    return 0;
}


