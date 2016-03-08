/*=========================================================================

Module    : Solid Mechanics
File      :
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
#ifndef FIBONACCIHEAP_H
#define FIBONACCIHEAP_H

/** \class smFibonacciHeap
 * \brief The classes in this package are designed to allow the package user
* to quickly and easily develop applications that require a heap data
* structure.  Using amortized analysis, the asymptotically fastest heap
* data structure is the Fibonacci heap.  The constants are a little
* high so the real speed gain will not be seen until larger data sets
* are required, but in most cases, if the data set is small, then the
* run-time will be neglible anyway.
*
* To use this heap class you need do only two things.  First, subclass
* the smFibonacciHeapNode class to create the class of objects that you'd
* like to store in a heap.  Second, create an instance of the smFibonacciHeap
* class, which can then be used to Insert(), ExtractMin(), etc.,
* instances of your smFibonacciHeapNode subclass.  Notice that you don't need
* to create a subclass of the smFibonacciHeap class.
*
* The application-specific data object that you'd like to store in a heap
* will have a key value.  In the class that you derive from smFibonacciHeapNode,
* you will need to define the key structure then provide assignment (=),
* equality (==) and less-than operators and a destructor.  These functions
* are declared virtual so that the code in the smFibonacciHeap class can compare,
* assign and destroy your objects by calling on your code.
*
* The overloaded operators in your defined class MUST call functions in
* the Fibonacci heap node class first.  For assignment, the function
* FHN_Assign() should be called before code that deals with the copy of
* the key value.  For comparison operators, the function FHN_Cmp() should
* appear first.  If it returns 0, then keys can be compared as normal.
* The following indicates what the three most common operators must do
* based on the return value of FHN_Cmp()
*
* For ==, if zero returned, then compare keys \n
*     if non-zero X returned, then return 0 \n
* For <,  if zero returned, then compare keys \n
*         if non-zero X returned, then return X<0?1:0 \n
* For >,  if zero returned, then compare keys \n
*         if non-zero X returned, then return X>0?1:0
*
* \verbatim The Fibonacci heap implementation contained in FIBHEAP.H and FIBHEAP.CPP
* is Copyright (c) 1996 by John Boyer
*
* Once this Fibonacci heap implementation (the software) has been published
* by Dr. Dobb's Journal, permission to use and distribute the software is
* granted provided that this copyright notice remains in the source and
* and the author (John Boyer) is acknowledged in works that use this program.
*
* Every effort has been made to ensure that this implementation is free of
* errors.  Nonetheless, the author (John Boyer) assumes no liability regarding
* your use of this software.
*
* The author would also be very glad to hear from anyone who uses the
* software or has any feedback about it.
* \e Email: \b jboyer@gulf.csc.uvic.ca
* \endverbatim
*/

#define OK      0
#define NOTOK   -1

#include <cstdio>
#include <cstdlib>
#include <iostream>

//======================================================
// Fibonacci Heap Node Class
//======================================================

class smFibonacciHeap;

class smFibonacciHeapNode
{
	friend class smFibonacciHeap;
	smFibonacciHeapNode *Left, *Right, *Parent, *Child;
	short Degree, Mark, NegInfinityFlag;

protected:
    /** To be used as the first step of ALL comparators in a derived class.
    *
    * Compares the relative state of the two neg. infinity
    * flags.  Note that 'this' is the left hand side.  If
    * LHS neg. infinity is set, then it will be less than (-1)
    * the RHS unless RHS neg. infinity flag is also set.
    * Only if function returns 0 should the key comparison
    * defined in the derived class be performed, e.g.
    *
    * For ==, if zero returned, then compare keys \n
    *     if non-zero X returned, then return 0 \n
    * For <,  if zero returned, then compare keys \n
    *         if non-zero X returned, then return X<0?1:0 \n
    * For >,  if zero returned, then compare keys \n
    *         if non-zero X returned, then return X>0?1:0 */
	inline int  FHN_Cmp(smFibonacciHeapNode& RHS)
    {
		if (NegInfinityFlag)
			return RHS.NegInfinityFlag ? 0 : -1;
		return RHS.NegInfinityFlag ? 1 : 0;
	}

	/** To be used as first step of an assignment operator in a
    * derived class.  The derived class will handle assignment
    * of key value, and this function handles copy of the
    * NegInfinityFlag (which overrides the key value if it is set).*/
	inline void FHN_Assign(smFibonacciHeapNode& RHS)
    {
		NegInfinityFlag = RHS.NegInfinityFlag;
	}

public:
	inline smFibonacciHeapNode()
    {
		Left = Right = Parent = Child = NULL;
		Degree = Mark = NegInfinityFlag = 0;
	}
	virtual ~smFibonacciHeapNode();
    
    /** We do, on occasion, compare and assign objects of type smFibonacciHeapNode, but
    * only when the NegInfinityFlag is set.  See for example smFibonacciHeap::Delete().
    *
    * Also, these functions exemplify what a derived class should do.*/
	virtual void operator =(smFibonacciHeapNode& RHS);
	virtual int  operator ==(smFibonacciHeapNode& RHS);
	virtual int  operator <(smFibonacciHeapNode& RHS);

	virtual void Print();
};

//========================================================================
// Fibonacci Heap Class
//========================================================================

class smFibonacciHeap
{
	smFibonacciHeapNode *MinRoot;
	long NumNodes, NumTrees, NumMarkedNodes;
	int  HeapOwnershipFlag;

public:
	smFibonacciHeap();
	virtual ~smFibonacciHeap();

	/** @name The Standard Heap Operations */
    //@{ 

    /** O(1) actual; O(2) amortized
    *
    * I am using O(2) here to indicate that although Insert() is
    * constant time, its amortized rating is more costly because some
    * of the work of inserting is done by other operations such as
    * ExtractMin(), which is where tree-balancing occurs.
    *
    * The child pointer is deliberately not set to NULL because Insert()
    * is also used internally to help put whole trees onto the root list.*/
	void Insert(smFibonacciHeapNode *NewNode);

    /** O(1) actual; O(1) amortized */
	void Union(smFibonacciHeap *OtherHeap);

    /** O(1) actual; O(1) amortized */
	inline smFibonacciHeapNode *Minimum()
    {   return MinRoot; }

	/** O(n) worst-case actual; O(lg n) amortized */
	smFibonacciHeapNode *ExtractMin();

    /** O(lg n) actual; O(1) amortized
    *
    * The O(lg n) actual cost stems from the fact that the depth, and
    * therefore the number of ancestor parents, is bounded by O(lg n).*/
	int DecreaseKey(smFibonacciHeapNode *theNode, smFibonacciHeapNode& NewKey);

    /** O(lg n) amortized; ExtractMin() dominates
    *
    * Notice that if we don't own the heap nodes, then we clear the
    * NegInfinityFlag on the deleted node.  Presumably, the programmer
    * will be reusing the node.*/
	int Delete(smFibonacciHeapNode *theNode);
    //@}


	/** @name Extra utility functions */
    //@{
	int  GetHeapOwnership() { return HeapOwnershipFlag; };
	void SetHeapOwnership() { HeapOwnershipFlag = 1; };
	void ClearHeapOwnership() { HeapOwnershipFlag = 0; };
	long GetNumNodes() { return NumNodes; };
	long GetNumTrees() { return NumTrees; };
	long GetNumMarkedNodes() { return NumMarkedNodes; };

    /** Used internally for debugging purposes.  The function prints the key
    * value for each node along the root list, then it calls itself on each
    * child list.*/
	void Print(smFibonacciHeapNode *Tree = NULL, smFibonacciHeapNode *theParent=NULL);
    //@}
    
private:
	// Internal functions that help to implement the Standard Operations
	inline void _Exchange(smFibonacciHeapNode* &N1, smFibonacciHeapNode* &N2)
    {
		smFibonacciHeapNode *Temp; Temp = N1; N1 = N2; N2 = Temp;
	}
	
	/** Internal function that reorganizes heap as part of an ExtractMin().
    * We must find new minimum in heap, which could be anywhere along the
    * root list.  The search could be O(n) since all nodes could be on
    * the root list.  So, we reorganize the tree into more of a binomial forest
    * structure, and then find the new minimum on the consolidated O(lg n) sized
    * root list, and in the process set each Parent pointer to NULL, and count
    * the number of resulting subtrees.
    *
    * Note that after a list of n inserts, there will be n nodes on the root
    * list, so the first ExtractMin() will be O(n) regardless of whether or
    * not we consolidate.  However, the consolidation causes subsequent
    * ExtractMin() operations to be O(lg n).  Furthermore, the extra cost of
    * the first ExtractMin() is covered by the higher amortized cost of
    * Insert(), which is the real governing factor in how costly the first
    * ExtractMin() will be.*/
	void _Consolidate();
    
	void _Link(smFibonacciHeapNode *, smFibonacciHeapNode *);
    
	void _AddToRootList(smFibonacciHeapNode *);

    /** Remove node x from the child list of its parent node y */
	void _Cut(smFibonacciHeapNode *, smFibonacciHeapNode *);

    /** Cuts each node in parent list, putting successive ancestor nodes on the
    * root list until we either arrive at the root list or until we find an
    * ancestor that is unmarked.  When a mark is set (which only happens during
    * a cascading cut), it means that one child subtree has been lost; if a
    * second subtree is lost later during another cascading cut, then we move
    * the node to the root list so that it can be re-balanced on the next
    * consolidate.*/
	void _CascadingCut(smFibonacciHeapNode *);
};

#endif // FIBONACCIHEAP_H
