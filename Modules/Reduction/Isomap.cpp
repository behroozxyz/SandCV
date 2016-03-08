/*=========================================================================

Module    : Solid Mechanics
File      : Isomap.cpp
Copyright : (C)opyright 2011++
            See COPYRIGHT statement in top level directory.
Authors   : D. Millan
Modified  :
Purpose   : 
Date      :
Version   :
Changes   :

    This software is distributed WITHOUT ANY WARRANTY; without even 
    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
    PURPOSE.  See the above copyright notices for more information.
=========================================================================*/
#include "Isomap.h"

Isomap::Isomap(int knn, int lowDim, int highDim, const vector<double> &xNodes)
    : nPts(xNodes.size()/highDim),
      tDim(lowDim),
      nDim(highDim),
      knn(knn),
      _searcher(knn+1, highDim, xNodes)
{
    for (int i=0; i<xNodes.size(); i++ )
    {
        _idxNN  = _searcher->FindNearestNeighbors( xNodes.at(i*nDim) );
        _ddxNN  = _searcher->GetNearestDistances();
    }
}

Isomap::Isomap()
{
    // Inputs
    this->nodeProj  = true;
    this->sampProj  = false;
    this->ierr      = 0;
    this->nDim      = 0;
    this->tDim      = 0;
    this->sPts      = 0;
    this->nPts      = 0;
    this->kNN       = 0;
    
    this->jn_knn    = NULL;
    this->in_knn    = NULL;
    this->js_knn    = NULL;
    this->is_knn    = NULL;

    this->in_nodes  = NULL;
    this->jn_nodes  = NULL;
    this->dn_nodes  = NULL;
    this->x_nodes   = NULL;
    this->x_samples = NULL;
    
    this->nears     = NULL;//false;
    this->nodes     = NULL;//false;
    this->samples   = NULL;//false;
    
    this->nears     = NULL;//false;
    this->geoDist   = NULL;//false;

    this->eigVal    = NULL;//array with the eigenvalues
    this->eigVec    = NULL;//array with the eigenvectors, these are stored consecutively [v1 v2 ... vd]
    
    // Outputs
    this->xi_nodes  = NULL;
    this->xi_samples= NULL;
    this->ierr      = 0;
}

Isomap::~Isomap()
{
    this->Clear();
}

//cleaning preallocated member variables
void Isomap::Clear()
{
    if ( this->eigVal     != NULL ) delete [] this->eigVal;
    if ( this->eigVec     != NULL ) delete [] this->eigVec;
    this->eigVal    = NULL;
    this->eigVec    = NULL;
    
    if ( this->xi_nodes   != NULL ) delete [] this->xi_nodes;
    if ( this->xi_samples != NULL ) delete [] this->xi_samples;
    this->xi_nodes  = NULL;
    this->xi_samples= NULL;
    this->ierr      = 0;
}

void Isomap::Update( )
{

    if ( this->CheckSettings() )
    {
        printf("ERROR:ISOMAP settings are wrong\n");
        this->ierr = 1;
        return;
    }

    printf("Step 1: Computing shortest paths\n");
    //Step 1: Computing shortest paths for the connected points
    if ( this->ierr == 0 )
        this->GeodesicDistance();

    printf("Step 2: Computing double centring matrix\n");
    // Step 2: Double Centering
    if ( this->ierr == 0 )
        this->DoubleCentering();
    
    printf("Step 3: Eigenvalue Problem\n");
    //Step 3: Eigenvalue Problem
    if ( this->ierr == 0 )
        this->EigenvalueProblem();
    
    printf("Step 4: Compute Embedding\n");
    //Step 4: Compute Embedding
    if ( this->ierr == 0 )
        this->Embedding();

    return;
}

//check the setting before to compute
int Isomap::CheckSettings()
{
    if ( this->nears.IsNotNull() )
    {
        this->nDim      = this->nears->GetDimension();
        this->nPts      = this->nears->GetNumberOfNodePoints();
        this->sPts      = this->nears->GetNumberOfSamplePoints();
        this->x_nodes   = this->nears->GetNodePoints();
        this->x_samples = this->nears->GetSamplePoints();
        this->dn_nodes  = this->nears->GetDnNodesAdjacency();
        this->jn_nodes  = this->nears->GetJnNodesAdjacency();
        this->in_nodes  = this->nears->GetInNodesAdjacency();
    }
    if ( this->nodes.IsNotNull() )
    {
        this->nDim      = this->nodes->GetDimension();
        this->nPts      = this->nodes->GetNumberOfPoints();
        this->x_nodes   = this->nodes->GetPoints();
    }
    if ( this->samples.IsNotNull() )
    {
        this->nDim      = this->samples->GetDimension();
        this->sPts      = this->samples->GetNumberOfPoints();
        this->x_samples = this->samples->GetPoints();
    }
    if ( this->nDim < 1 )
    {
        printf("ERROR::Spatial dimension is %d < 1\n", this->nDim);
        return 1;
    }
    if ( this->tDim < 1 )
    {
        printf("ERROR::Spatial dimension for the tangent projection is %d < 1\n", this->tDim);
        return 1;
    }
    if ( this->tDim > this->nDim )
    {
        printf("ERROR::Spatial dimension for the tangent projection is bigger than the spatial one [tDim=%d > %d=nDim]\n", this->tDim, this->nDim);
        return 1;
    }
    if ( this->nPts < 1 )
    {
        printf("ERROR::Number of nodes either is not defined or is less than 1 (nPts=%d)\n",
            this->nPts);
        return 1;
    }
    if ( this->x_nodes == NULL )
    {
        printf("ERROR::Node points coordinates are not set\n");
        return 1;
    }
        if (this->dn_nodes == NULL)
    {
        printf("ERROR::Distances between each node and its nearest neighbors has not been defined\n");
        return 1;
    }
    if (this->jn_nodes == NULL)
    {
        printf("ERROR::List indices with the nodes closest to each node point is not set\n");
        return 1;
    }
    if (this->in_nodes == NULL)
    {
        printf("ERROR::Index list for counting the number of nodes at each node point is not set\n");
        return 1;
    }
    for ( int i=0; i<this->nPts; i++ )
    {
        for ( int j=this->in_nodes[i]; j<this->in_nodes[i+1]; j++ )
        {
            if ( this->jn_nodes[j] == i )
            {
                printf("ERROR::The index list is non well defined for the ith-node ('i' is present)\n");
                return 1;
            }
        }
    }
    if ( this->sPts < 1 && this->sampProj == true )
    {
        printf("ERROR::Number of sample points is not defined or is less than 1 (sPts=%d)\n", this->sPts);
        return 1;
    }

    if ( this->x_samples == NULL && this->sampProj == true )
    {
        printf("ERROR::Sample points coordinates are not specified\n");
        return 1;
    }
    
    // added in order to avoid memory leak, i.e. for recalling the same object with
    // different parameters as the dimension (see SimplexRegion.cpp example)
    this->Clear();

    return 0;
}

//step 2
void Isomap::GeodesicDistance()
{
    this->geoDist = smGeodesicDistance::New();
    this->geoDist->SetAdjacencyStructure(this->nears);
    this->geoDist->SetNumberOfNearestNeighbours(this->nPts-1);
    this->geoDist->SetSelfIdentifierON();
    this->geoDist->SetUpperStorageON();
    this->geoDist->Update();

    if ( this->geoDist->GetError() == 0 )
        this->M = this->geoDist->GetGeodesicDistanceMatrix();
    else
        this->ierr = 2;

    return;
}

// Double center of a distance matrix: Classical MDS
// Center the matrix by substracting the mean of each column, row,
// adding the overall mean, and multiplying by -0.5
//DD = -0.5*[D.^2 - sum(D.^2)'*ones(1,N)/N - ones(N,1)*sum(D.^2)/N + sum(sum(D.^2))/(N^2)]
void Isomap::DoubleCentering()
{
    int    i, j, jj;
    int    N, NNZ;
    int    iInit, iEnd, iNN;
    double fullmean;
    double *row_means;

    N   = this->M->n;
    NNZ = this->M->nnz;

    //[1] Square of the distance matrix
    //    D -> D.^2
    for ( i=0; i<NNZ; i++ )
        this->M->an[i] *= this->M->an[i];

    //[2] average the rows
    //    D -> -sum(D.^2,2)*ones(1,N)/N
    row_means = new double [N];
    smArray::Fill(N, 0.0, row_means);
    for ( i=0; i<N; i++ )
    {
        iInit = this->M->ia[i];
        iNN   = this->M->ia[i+1] - iInit;
        row_means[i] += smMath::Sum(iNN, &this->M->an[iInit]);
        for ( j=0; j<i; j++ )
        {
            jj = this->M->ia[j] + i - j;
            row_means[i] += this->M->an[jj];
        }
    }
    for ( i=0; i<N; i++ )
        row_means[i] /= ((double)N);

    //[3] average the columns
    //    D -> -ones(N,1)*sum(D.^2,1)/N
    // NOTE its are equal to the row mean because of the GeodesicDistance matrix is symmetric

    //[4] full average
    //    D -> sum(sum(D.^2))/(N^2)
    // NOTE because the matrix is symmetric the mean is computed in this way, to avoid
    // a double counting of the diagonal terms. The diagonal terms are zero.
    fullmean  = 2.0*smMath::Sum(NNZ, this->M->an);
    fullmean /= (double)(2*NNZ - N);

    //[5] D -> *(-0.5)
    for ( i=0; i<N; i++ )
    {
        iInit = this->M->ia[i];
        iEnd  = this->M->ia[i+1];
        for ( jj=iInit; jj<iEnd; jj++ )
        {
            j = this->M->ja[jj];
            this->M->an[jj] += fullmean - row_means[i] - row_means[j];
            this->M->an[jj] *= (-0.5);
        }
    }
    
    delete [] row_means;
    return;
}


//Step 3: Eigenvalue Problem
// we need to find the smallest eigenvalues of the system M = (I-W)' (I-W)
// and then map the data to the low dimensional embedding defined by its eigenvectors
void Isomap::EigenvalueProblem()
{     
    // Solve eigenproblem: dspevx
    // Computes selected eigenvalues and eigenvectors of a real symmetric matrix in packed storage.
    
    //INPUT PARAMETERS
    int    mat_layout;
    char   jobz;
    char   range;
    char   uplo;
    int    n;           //The order of the matrix A (n ≥ 0)
    double *ap;         //ap stores the packed upper triangular part of A
    double vl, vu;      //lower and upper bounds of the interval to be searched for eigenvalues
    int    il, iu;      //index in ascending order of the smallest/largest eigenvalues to be returned
    double abstol;      //the absolute error tolerance to which each eigenvalue is required
    int    ldz;         //the leading dimension of the output array z, max(1,n)

    //OUTPUT PARAMETERS
    //ap: this array is overwritten by the values of the reduction to tridiagonal form.
    int     m;          //total number of eigenvalues found, if range = 'I', m = iu-il+1.
    double *w;          //selected eigenvalues of the matrix A in ascending order. Size max(1, n)
    double *z;          //the first m columns of z contain the orthonormal eigenvectors of 
                        //the matrix A corresponding to the selected eigenvalues, with the
                        //i-th column of z holding the eigenvector associated with w(i).
                        //Size of z is max(1, m).
    int    *ifail;      //contains the index of the eigenvectors that failed to converge    
    int    info;        //error identifier:
                        //  If info = 0, the execution is successful.
                        //  If info = -i, the i-th parameter had an illegal value.
                        //  If info = i, then i eigenvectors failed to converge; their indices are stored
                        //  in the array ifail.
    
    this->eigVal  = new double [this->tDim];                //eigenvalues
    this->eigVec  = new double [this->nPts*(this->tDim+1)]; //egenvectors, stored as [v1 v2 ... vd]
    
    //Initialization
    mat_layout = LAPACK_COL_MAJOR; //row-major ordering (C default in language)
    jobz       = 'V';              //eigenvalues and eigenvectors are computed
    range      = 'I';              //the routine computes eigenvalues with index il to iu
    uplo       = 'L';              //ap stores the packed upper triangular part of A
    n          = this->M->n;       
    ap         = this->M->an;      //packed upper triangle of the symmetric matrix A
    vl         = 1.E0;             //lower bound of the interval to be searched for eigenvalues (vl<vu)
    vu         = 1.E8;             //upper bound of the interval to be searched for eigenvalues
    il         = n-this->tDim;     //index for the smallest eigenvalue to be returned
    iu         = n;                //index for the larger eigenvalue to be returned
    abstol     = 1.E-12;           //absolute error tolerance to which each eigenvalue is required
    w          = new double[n];    //eigenvalues
    z          = this->eigVec;     //eigenvectors in column wise
    ldz        = n;                //the leading dimension of the output array z
    ifail      = new int [n];
    for ( int i=0; i<n; i++ )
        ifail[i] = -1;
    
    info = LAPACKE_dspevx( mat_layout, jobz, range, uplo, n, ap, vl, vu, il, iu, abstol,
                           &m, w, z, ldz, ifail );
    
    // Check for convergence
    if( info != 0 )
    {
        printf( "The algorithm failed to solve the eigenvalue problem.\n" );
        if ( info < 0 )
            printf( "The %d-th parameter of LAPACKE_dspevx had an illegal value\n", ABS(info));
        else
            printf( "There are %d eigenvectors that failed to converge. Their indices are:\n", info);
        for ( int i=0; i<n; i++ )
        {
            if (ifail[i] > 0 )
                printf("eigenvector %4d - ifail=[%4d]\n", i, ifail[i]);
        }
        this->ierr = 3;
    }
    else
    {
        //printf("The total number of eigenvalues found is m=%d\n", m);
        smArray::Copy(this->tDim, &w[1], this->eigVal);
        //printf("eval : ");smIO::PrintVector(this->tDim, this->eigVal);
    }
    
    delete [] w;
    delete [] ifail;
    return;
}

// The function returns the low-dimensional representation of the data x_nodes
// in matrix form
void Isomap::Embedding()
{
    register int i, j, k, iD;

    double  val;       //temporal auxiliary value

    //The mapped points are extracted from the eigen-analysis
    this->xi_nodes = new double [this->nPts*this->tDim];

    for ( j=0; j<this->tDim; j++ )
    {
        if ( std::isnan(this->eigVal[j]) == 1 || std::isinf(this->eigVal[j]) == 1 )
        {
            printf("ERROR::Eigenvalues from the ISOMAP problem are NaN or Inf val[%d]=%f\n",
                   j, this->eigVal[j]);
            this->ierr = 11;
            j = this->tDim;
        }
        else if ( this->eigVal[j] < 0.0 )
        {
            printf("ERROR::Eigenvalue from the ISOMAP problem is negative: val[%d]=%f\n",
                   j, this->eigVal[j]);
            printf("eval : ");smIO::PrintVector(this->tDim, this->eigVal);
            this->ierr = 12;
            j = this->tDim;
        }
    }

    if ( this->tDim > 0 && this->ierr == 0 )
    {
        //scaling embedded dimensional coordinates
        for ( i=0; i<this->nPts; i++ )
        {
            iD = i*this->tDim;
            for ( j=0; j<this->tDim; j++ )
            {
                //each dimension is scaled
                //val = this->eigVec[j*this->nPts+i]*sqrt(this->eigVal[j]);
                val = this->eigVec[(j+1)*this->nPts+i]*sqrt(this->eigVal[j]);
                if ( std::isnan(val) == 1 || std::isinf(val) == 1 )
                {
                    printf("ERROR::Values from the eigenvectors are NaN or Inf val[%d,%d]=%f\n",
                        i, j, val);
                    this->ierr = 13;
                    j = this->tDim;
                    i = this->nPts;
                }
                else
                {
                    k   = this->tDim-1-j;
                    this->xi_nodes[iD+k] = val;
                }
            }
        }
    }
    return;
}
