%module Interface
%include "carrays.i"
%include "std_vector.i"
%{
#include "ZoneLock.h"
#include "PuckerABF.h"
#include "SandABF.h"
#include "SandSMD.h"
#include "SandCMapSMD.h"
#include "MPSandABF.h"
#include "MPSandRest.h"
%}

// Instantiate templates used by example
namespace std {
   %template(IntVector) vector<int>;
   %template(DoubleVector) vector<double>;
}

%array_functions(int   , IntArray   );
%array_functions(double, DoubleArray);

class ZoneLock
{
public:
    ZoneLock(int nAtoms, int rdim, int ldim,
             char *xRefFilename, int nInvolvedAtoms, int *involvedAtoms,
             char *zoneCenterFilename, double kBias, double bBias);

    void CalculateForces(int stepNum, double *masses, double *positions, double *forces);
};



class PuckerABF
{
public:
    PuckerABF(int bApplyABF,
             // simulation information
             int nAtoms, int rdim, double timeStep,
             int logFrequency, int histFrequency,
             char *coreName,
             // corral potential
             double corralHeight, char *corralFilename,
             // histogram
             int bLoadHistogram , double *xiMin, double *xiMax, int *nBins
             , int fullSample, int rampSample);

    PuckerABF(int bApplyABF,
             // simulation information
             int nAtoms, int rdim, double timeStep,
             int logFrequency, int histFrequency,
             char *coreName,
             // semi-harmonic potential
             double kWall, double *lowerWall, double *upperWall,
             // histogram
             int bLoadHistogram , double *xiMin, double *xiMax, int *nBins
             , int fullSample, int rampSample);

    void CalculateForces(int stepNum, double *masses, double *positions, double *forces);
};


class SandABF
{
public:
    SandABF( int bApplyABF,
             // simulation information
             int nAtoms, int rdim, double timeStep,
             int logFrequency, int histFrequency,
             char *coreName,
             // alignment
             char *xRefFilename, int nInvolvedAtoms, int *involvedAtoms,
             // parametrization
             int ldim, int bDistance,
             char *xiNodesFilename, char *xNodesFilename,
             char *xiSeedsFilename,
             // corral potential
             double corralHeight, char *corralFilename,
             // histogram
             int bLoadHistogram , double *xiMin, double *xiMax, int *nBins
             );

    SandABF( int bApplyABF,
             // simulation information
             int nAtoms, int rdim, double timeStep,
             int logFrequency, int histFrequency,
             char *coreName,
             // alignment
             char *xRefFilename, int nInvolvedAtoms, int *involvedAtoms,
             // parametrization
             int ldim, int bDistance,
             char *xiNodesFilename, char *xNodesFilename,
             char *xiSeedsFilename,
             // semi-harmonic potential
             double kWall, double *lowerWall, double *upperWall,
             // histogram
             int bLoadHistogram , double *xiMin, double *xiMax, int *nBins
             );

    void CalculateForces(int stepNum, double *masses, double *positions, double *forces);
};


class SandSMD
{
public:
    SandSMD(int bApplySMD,
                 // simulation information
                 int nAtoms, int rdim, double timeStep,
                 int logFrequency,
                 char *coreName,
                 // alignment
                 char *xRefFilename, int nInvolvedAtoms, int *involvedAtoms,
                 // parametrization
                 int ldim, int bDistance,
                 char *xiNodesFilename, char *xNodesFilename,
                 char *xiSeedsFilename,
                 // Pulling SMD
                 int bWithReturn,
                 double pullingRate, double pullingConstant, double pullingMaxLength,
                 double kWall, double *lowerWall, double *upperWall,
                 double xiMin, double xiMax
                 );

    void CalculateForces(double *positions, double *forces);
};



class SandCMapSMD
{
public:
    SandCMapSMD(int bApplySMD,
                 // simulation information
                 int nAtoms, int sdim, double timeStep,
                 int logFrequency,
                 char *coreName,
                 // contact map
                 double d0,
                 // parametrization
                 int ldim, int bDistance,
                 char *xiNodesFilename, char *xNodesFilename,
                 char *xiSeedsFilename,
                 // Pulling SMD
                 int bWithReturn,
                 double pullingRate, double pullingConstant, double pullingMaxLength,
                 double kWall, double *lowerWall, double *upperWall,
                 double xiMin, double xiMax
                 );

    void CalculateForces(double *positions, double *forces);
};


class MPSandABF
{
public:
    MPSandABF(int nPartitions, int bApplyABF,
              // simulation information
              int nAtoms, int rdim, double timeStep,
              int logFrequency, int histFrequency,
              char *coreName,
              // alignment
              char *xRefFilename, int nInvolvedAtoms, int *involvedAtoms,
              // parametrization
              int ldim, int bDistance,
              char *xiNodesFilename, char *xNodesFilename,
              char *xiSeedsFilename,
              // Markers
              char *xMarkersFilename,
              // histogram
              int bLoadHistogram , double *xiMin, double *xiMax, int *nBins,
              double gammaPU, double spacingPU, double tol0PU);

    ~MPSandABF();

    void CalculateForces(int stepNum, double *masses, double *positions, double *forces);

    void CalculateUniqueForces(int stepNum, double *masses, double *positions, double *forces);
};


class MPSandRest
{
public:
    MPSandRest(int nPartitions,
               // simulation information
               int nAtoms, int rdim,
               int logFrequency,
               // alignment
               char *xRefFilename, int nInvolvedAtoms, int *involvedAtoms,
               // parametrization
               int ldim, int bDistance,
               char *xiNodesFilename, char *xNodesFilename,
               char *xiSeedsFilename,
               // Markers
               char *xMarkersFilename,
               // Bias force
               char *xCenterFilename, double kBias, double rBias,
               double *xiMin, double *xiMax, double gammaPU, double spacingPU, double tol0PU);

    void CalculateForces(int stepNum, double phiAngle, double psiAngle, double *masses, double *positions, double *forces);
};
