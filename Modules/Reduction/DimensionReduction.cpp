/*=========================================================================

  Module    : Solid Mechanics
  File      : smDimensionReduction.cpp
  Copyright : (C)opyright 2009++
              See COPYRIGHT statement in top level directory.
  Authors   : D. Millan, A. Rosolen
  Modified  :
  Purpose   :
  Date      :
  Version   :
  Changes   :

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.
=========================================================================*/
#include "smDimensionReduction.h"

// Constructor
smDimensionReduction::smDimensionReduction()
{
	// Inputs
	this->nodeProj  = true;
	this->sampProj  = false;
	this->ierr		= 0;
	this->nDim		= 0;
	this->tDim		= 0;
	this->sPts		= 0;
	this->nPts		= 0;
	this->kNN		= 0;

	this->jn_knn	= NULL;
	this->in_knn	= NULL;
	this->js_knn	= NULL;
	this->is_knn	= NULL;

	this->x_nodes   = NULL;
	this->x_samples = NULL;
	
	this->nears     = NULL;//false;
	this->nodes     = NULL;//false;
	this->samples   = NULL;//false;

	// Outputs
	this->xi_nodes   = NULL;
	this->xi_samples = NULL;
	this->ierr      = 0;
}

//Destructor
smDimensionReduction::~smDimensionReduction()
{
}
