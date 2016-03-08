/*=========================================================================

  Module    : Solid Mechanics
  File      :
  Copyright : (C)opyright 2011++
              See COPYRIGHT statement in top level directory.
  Authors   : D. Millan
  Modified  :
  Purpose   : Dijkstra algorithm to compute the graph geodesic.
  Date      :
  Version   :
  Changes   :

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.
=========================================================================*/

#include "GeodesicDistance.h"

GeodesicDistance::GeodesicDistance()
{
    this->upper      = true;
    this->selfId     = true;
    this->ierr       = 0;
    this->nPts       = 0;
    this->in_nodes   = NULL;
    this->jn_nodes   = NULL;
    this->dn_nodes   = NULL;
    this->kNN_       = 0;
    this->kNN        = 0;
    this->nears      = NULL;

    // sparse Geodesic distance matrix
    this->Gmat.nnz   = 0;
    this->Gmat.m     = 0;
    this->Gmat.n     = 0;
    //this->Gmat.ia_diag = NULL;
    this->Gmat.ia    = NULL;
    this->Gmat.ja    = NULL;
    this->Gmat.an    = NULL;
}

GeodesicDistance::~GeodesicDistance()
{
    this->Clear();
}

void GeodesicDistance::Clear()
{
    if ( this->ierr != 2 )
    {
        if ( this->Gmat.ia    != NULL ) delete [] this->Gmat.ia;
        if ( this->Gmat.ja    != NULL ) delete [] this->Gmat.ja;
        if ( this->Gmat.an    != NULL ) delete [] this->Gmat.an;
        this->Gmat.ia    = NULL;
        this->Gmat.ja    = NULL;
        this->Gmat.an    = NULL;
    }
}

void GeodesicDistance::Update()
{
    if ( this->CheckSettings() )
    {
        printf("ERROR:DijkstraGraph settings are wrong\n");
        this->ierr = 1;
        return;
    }


    //in the case that kNN=(nPts-1) and selfId=true the full Geodesic Distance Matrix is computed
    if ( this->selfId == true && this->kNN == (this->nPts-1))
    {
        this->ShortestFullPath();
        this->kNN_ = this->kNN+1;
    }
    else
    {
        this->ShortestPartialPath();
        this->kNN_ = this->kNN;
    }
    return;
}


int GeodesicDistance::CheckSettings()
{
    //basic checking for the input settings
    if ( this->nears.IsNotNull() )
    {
        this->nPts      = this->nears->GetNumberOfNodePoints();
        this->dn_nodes  = this->nears->GetDnNodesAdjacency();
        this->jn_nodes  = this->nears->GetJnNodesAdjacency();
        this->in_nodes  = this->nears->GetInNodesAdjacency();
    }
    if ( this->nPts < 1 )
    {
        printf("ERROR::Number of nodes either is not defined or is less than 1 (nPts=%d)\n",
               this->nPts);
        return 1;
    }
    if ( (this->kNN < 1) || (this->kNN >= this->nPts) )
    {
        printf("ERROR::Number of k-nearest neighbors is bad defined (kNN=%d)\n", this->kNN);
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
            if ( this->dn_nodes[j] <= EPSILON )
            {
                printf("ERROR::The distance between nearest neighbors is lower than EPSILON\n");
                return 1;
            }
        }
    }

    // the following lines are added in order to avoid memory leak, i.e. for recalling the
    // same object with different parameters
    this->Clear();

    //All is OK
    return 0;
}

void GeodesicDistance::PartialDijkstra( int s, int *rangeID, int *mark_neigh, int *Jn, double *Dn)
{
    int     finished;
    int     i, startInd, endInd;
    int     neighID, closest, minID, maxID;
    int     nDone, counter;
    double  closestD, arcLength, oldDist;
    smHeapNode *A, *hnMin, hnTmp;
    smFibonacciHeap *heap;

    // setup heap
    if ((heap = new smFibonacciHeap) == NULL || (A = new smHeapNode[this->nPts+1]) == NULL )
    {
        printf( "Memory allocation failed-- ABORTING.\n" );
        this->ierr = 2;
        return;
    }
    heap->ClearHeapOwnership();

    // initialize
    for ( i=0; i<this->nPts; i++ )
    {
        if ( i!=s )
            A [i] = (double) INFINITY;
        else
            A [i] = (double) EPSILON;
        if ( i!=s )
            Dn[i] = (double) INFINITY;
        else
            Dn[i] = (double) EPSILON;
        Jn[i] = -1;
        heap->Insert( &A[i] );
        A[i].SetIndexValue( i );
    }

    // Insert 0 then extract it, which will balance heap
    heap->Insert(&hnTmp);
    heap->ExtractMin();

    smArray::Fill(this->nPts, 0, mark_neigh);

    // loop over nonreached nodes
    finished = 0;
    nDone    = 0;
    counter  = 0;
    minID    = s;
    maxID    = s;
    while ( (finished==0) && (nDone < this->nPts) )
    {
        hnMin    = (smHeapNode *) heap->ExtractMin();
        closest  = hnMin->GetIndexValue();
        closestD = hnMin->GetKeyValue();
        if ( (closest<0) || (closest>=this->nPts) )
        {
            printf( "ERROR::Minimum Index out of bound...\n");
            this->ierr = 3;
            return;
        }
        Dn[closest] = closestD;
        if ( closestD == INFINITY )
            finished = 1;
        else
        {
            // relax all nodes adjacent to closest
            nDone++;
            startInd = this->in_nodes[closest  ];
            endInd   = this->in_nodes[closest+1];
            for( i=startInd; i<endInd; i++ )
            {
                neighID   = this->jn_nodes[i];
                arcLength = this->dn_nodes[i];
                oldDist   = Dn[neighID];
                if ( oldDist > ( closestD + arcLength ))
                {
                    //printf("closest=%2d  j=%2d  dist=%f -> %f\n", closest+1, neighID+1, oldDist, closestD+arcLength);
                    Dn[neighID] = closestD + arcLength;
                    Jn[neighID] = closest + 1;
                    hnTmp = A[ neighID ];
                    hnTmp.SetKeyValue( closestD + arcLength );
                    heap->DecreaseKey( &A[ neighID ], hnTmp );

                    if ( mark_neigh[neighID] == 0 )
                    {
                        mark_neigh[neighID] = 1;
                        minID = MIN(minID, neighID);
                        maxID = MAX(maxID, neighID);
                        counter++;

                        if ( counter > 2*this->kNN ) //it must be kNN and NOT kNN_
                            finished = 1;
                    }
                }
            }//endfor
        }//endif
    }

    rangeID[0] = minID;
    rangeID[1] = maxID;

    // cleanup
    delete heap;
    delete [] A;
    return;
}


void GeodesicDistance::ShortestPartialPath()
{

    int i, j;
    int start;

    //in the case that kNN=(nPts-1) and selfId=true the full Geodesic Distance Matrix is computed
    if ( this->selfId == true )
    {
        start = 0;
        this->kNN_ = this->kNN+1;
    }
    else
    {
        start = 1;
        this->kNN_ = this->kNN;
    }

    this->Gmat.n   = this->nPts;
    this->Gmat.m   = this->nPts;
    this->Gmat.nnz = this->nPts*this->kNN_;
    
    if ( this->Gmat.nnz < 1 )
    {
        printf("ERROR::The number of non-zeros for the distance matrix is bad defined NNZ=%d (<1)\n",
               this->Gmat.nnz);
        this->ierr = 2;
        return;
    }
    
    this->Gmat.ia  = new    int [this->nPts+1];
    this->Gmat.ja  = new    int [this->Gmat.nnz];
    this->Gmat.an  = new double [this->Gmat.nnz];

    this->Gmat.ia[0] = 0;
    for( i=0; i<this->nPts; i++ )
        this->Gmat.ia[i+1] = this->Gmat.ia[i] + this->kNN_;

    smArray::Fill(this->Gmat.nnz, 0.0, this->Gmat.an);

#ifdef _OPENMP
    int          nPts = this->nPts;
    int          kNN_ = this->kNN_;
    SparseMatrix Gmat = this->Gmat;
#pragma omp parallel               \
    default(none)                  \
    private(i, j)                  \
    shared (start, nPts, kNN_, Gmat)
#endif
    {
        int    minID, maxID;
        // allocate memory for single source results (automatically recycled)
        int    rangeID[2];
        int    *mark_neigh = new    int [nPts];
        int    *Jn_neigh   = new    int [nPts];
        double *Dn_neigh   = new double [nPts];
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
        // loop over sources
        for( i=0; i<nPts; i++ )
        {
            // run the Dijkstra code for single source (0 indexed)
            PartialDijkstra( i, rangeID, mark_neigh, Jn_neigh, Dn_neigh );

            minID = rangeID[0];
            maxID = rangeID[1];

            //sort values
            //  wa2[]   - the work array whose elements are to be sorted
            //  iwa2[]  - the associated index work array to the wa2 array
            //  N       - index of the last element to be sorted
            for( j=minID; j<=maxID; j++ )
                Jn_neigh[j] = j;

            smHeapSort::HeapSort( maxID-minID+1, &Dn_neigh[minID], &Jn_neigh[minID] );
            smHeapSort::HeapSort( kNN_, &Jn_neigh[minID+start], &Dn_neigh[minID+start] );

            // store results in the full distance matrix
            smArray::Copy(kNN_, &Jn_neigh[minID+start], &Gmat.ja[i*kNN_]);
            smArray::Copy(kNN_, &Dn_neigh[minID+start], &Gmat.an[i*kNN_]);
        }
        delete [] mark_neigh;
        delete [] Jn_neigh;
        delete [] Dn_neigh;
    }
    return;
}



void GeodesicDistance::FullDijkstra( int s, int *Jn, double *Dn)
{
    int     finished;
    int     i, startInd, endInd;
    int     neighID, closest;
    int     nDone;//, counter;
    double  closestD, arcLength, oldDist;
    smHeapNode *A, *hnMin, hnTmp;
    smFibonacciHeap *heap;

    // setup heap
    if ((heap = new smFibonacciHeap) == NULL || (A = new smHeapNode[this->nPts+1]) == NULL )
    {
        printf( "Memory allocation failed-- ABORTING.\n" );
        this->ierr = 2;
        return;
    }
    heap->ClearHeapOwnership();

    // initialize
    for ( i=0; i<this->nPts; i++ )
    {
        if ( i!=s )
            A [i] = (double) INFINITY;
        else
            A [i] = (double) EPSILON;
        if ( i!=s )
            Dn[i] = (double) INFINITY;
        else
            Dn[i] = (double) EPSILON;
        Jn[i] = -1;
        heap->Insert( &A[i] );
        A[i].SetIndexValue( i );
    }

    // Insert 0 then extract it, which will balance heap
    heap->Insert(&hnTmp);
    heap->ExtractMin();

    // loop over nonreached nodes
    finished = 0;
    nDone    = 0;
    //counter  = 0;
    while ( (finished==0) && (nDone < this->nPts) )
    {
        hnMin    = (smHeapNode *) heap->ExtractMin();
        closest  = hnMin->GetIndexValue();
        closestD = hnMin->GetKeyValue();
        if ( (closest<0) || (closest>=this->nPts) )
        {
            printf( "ERROR::Minimum Index out of bound...\n");
            this->ierr = 3;
            return;
        }
        Dn[closest] = closestD;
        if ( closestD == INFINITY )
            finished = 1;
        else
        {
            // relax all nodes adjacent to closest
            nDone++;
            startInd = this->in_nodes[closest  ];
            endInd   = this->in_nodes[closest+1];
            for( i=startInd; i<endInd; i++ )
            {
                neighID   = this->jn_nodes[i];
                arcLength = this->dn_nodes[i];
                oldDist   = Dn[neighID];
                if ( oldDist > ( closestD + arcLength ))
                {
                    //printf("closest=%2d  j=%2d  dist=%f -> %f\n", closest+1, neighID+1, oldDist, closestD+arcLength);
                    Dn[neighID] = closestD + arcLength;
                    Jn[neighID] = closest + 1;
                    hnTmp = A[ neighID ];
                    hnTmp.SetKeyValue( closestD + arcLength );
                    heap->DecreaseKey( &A[ neighID ], hnTmp );
                }
            }//endfor
        }//endif
    }

    // cleanup
    delete heap;
    delete [] A;
    return;
}


void GeodesicDistance::ShortestFullPath()
{

    int i, j, jj;

    //in the case that kNN=(nPts-1) and selfId=true the full Geodesic Distance Matrix is computed
    this->kNN_ = this->nPts;

    this->Gmat.n   = this->nPts;
    this->Gmat.m   = this->nPts;
    if (this->upper == true)
        this->Gmat.nnz = (int)this->nPts*((this->nPts+1)*0.5);
    else
        this->Gmat.nnz = this->nPts*this->nPts;

    if ( this->Gmat.nnz < 1 )
    {
        printf("ERROR::The number of non-zeros for the distance matrix is bad defined NNZ=%d (<1)\n",
               this->Gmat.nnz);
        this->ierr = 2;
        return;
    }
    this->Gmat.ia  = new    int [this->nPts+1];
    this->Gmat.ja  = new    int [this->Gmat.nnz];
    this->Gmat.an  = new double [this->Gmat.nnz];

    this->Gmat.ia[0] = 0;
    for( i=0; i<this->nPts; i++ )
    {
        this->Gmat.ia[i+1] = this->Gmat.ia[i] + this->nPts;
        if ( this->upper==true )
        {
            this->Gmat.ia[i+1] -= i;
            for ( j=i, jj=0; j<nPts; j++, jj++ )
                this->Gmat.ja[this->Gmat.ia[i]+jj] = j;
        }
    }
    smArray::Fill(this->Gmat.nnz, 0.0, this->Gmat.an);

#ifdef _OPENMP
    bool         upper = this->upper;
    int          nPts  = this->nPts;
    SparseMatrix Gmat  = this->Gmat;
#pragma omp parallel               \
    default(none)                  \
    private(i, j, jj)          \
    shared (upper, nPts, Gmat)
#endif
    {
        // allocate memory for single source results (automatically recycled)
        int    *Jn = new    int [nPts];
        double *Dn = new double [nPts];
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
        // loop over sources
        for( i=0; i<nPts; i++ )
        {
            // run the Dijkstra code for single source (0 indexed)
            FullDijkstra( i, Jn, Dn );

            //sort values
            //  wa2[]   - the work array whose elements are to be sorted
            //  iwa2[]  - the associated index work array to the wa2 array
            //  N       - index of the last element to be sorted
            for( j=0; j<nPts; j++ )
                Jn[j] = j;

            smHeapSort::HeapSort( nPts, Dn, Jn );
            smHeapSort::HeapSort( nPts, Jn, Dn );

            // store results in the full distance matrix
            if ( upper == true )
            {
#ifdef _OPENMP
#pragma omp critical //section of code that must be executed by a single thread at a time
#endif
                {
                    //triangular lower part
                    for( j=0; j<i; j++ )
                    {
                        jj = Gmat.ia[j] + i - j;
                        Gmat.an[jj] += Dn[j];
                    }

                    //diagonal and triangular upper part
                    for( j=i, jj=Gmat.ia[i]; j<nPts; j++, jj++ )
                    {
                        Gmat.an[jj] += Dn[j];
                    }

                }
            }
            else
            {
                smArray::Copy(nPts, Jn, &Gmat.ja[i*nPts]);
                smArray::Copy(nPts, Dn, &Gmat.an[i*nPts]);
            }
        }
        delete [] Jn;
        delete [] Dn;
    }

    //averaging the computed distances
    if ( this->upper == true )
    {
        for( i=0; i<this->Gmat.nnz; i++ )
            this->Gmat.an[i] *= 0.5;
    }

    return;
}
