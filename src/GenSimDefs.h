
#ifndef ARRAYTYPEDEFS_H_
#define ARRAYTYPEDEFS_H_

/** Parallelization options
 *
 * 		OMPHYPER:		parallelize all the double-for loops of the Hyperbolic Solver
 * 		OMPELLIPTIC:	parallelize all the double-for loops of the Elliptic Solver
 */
/*
	#define OMPHYPER
	#define OMPELLIPTIC
	#define OMPLEVELSET
	#define OMPIMPEC
	#define OMPUTIL
	#define OMPPROJECTOR
*/
#define SPLIT2ND


	// Enable thread-safety
	#define BZ_THREADSAFE
	#include "blitz/array.h"

	typedef blitz::TinyVector<double, 6> CurvVec;
	typedef blitz::Array<CurvVec,1>  CurvArray;
	typedef blitz::Array<double, 2> Array2d;		    
	typedef blitz::Array<double, 3> Array3d;		    
	typedef blitz::Array<bool, 2> Mask2d;
	typedef blitz::TinyVector<int,2> Idx2d;
	typedef blitz::TinyVector<int,3> Idx3d;
#endif

