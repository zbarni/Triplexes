// ==========================================================================
//                                triplexator
// ==========================================================================
// Copyright (c) 2011,2012, Fabian Buske, UQ
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Fabian Buske or the University of Queensland nor 
//       the names of its contributors may be used to endorse or promote products 
//       derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Fabian Buske <fbuske@uq.edu.au>
// ==========================================================================

#ifndef FBUSKE_APPS_TRIPLEXATOR_HEADER_GARDENER_H
#define FBUSKE_APPS_TRIPLEXATOR_HEADER_GARDENER_H

#include <limits>
#include "find_index_qgrams.h"
#include "triplex_alphabet.h"
#include "helper.h"
#include "edlib.h"
#include <seqan/seeds2.h>  // Include module under test.
#include <seqan/sequence/adapt_std_list.h>
#include <seqan/misc/priority_type_base.h>
#include <seqan/misc/priority_type_heap.h>
#include <seqan/misc/misc_dequeue.h>

#if SEQAN_ENABLE_PARALLELISM
#include <seqan/parallel.h>
#endif  // #if SEQAN_ENABLE_PARALLELISM

#ifndef SEQAN_PRAGMA_IF_PARALLEL
#if SEQAN_ENABLE_PARALLELISM
#define STRINGIFY(a) #a
#define SEQAN_PRAGMA_IF_PARALLEL(code) \
_Pragma(STRINGIFY(code))
#else // SEQAN_ENABLE_PARALLELISM
#define SEQAN_PRAGMA_IF_PARALLEL(code)
#endif // SEQAN_ENABLE_PARALLELISM
#endif // SEQAN_PRAGMA_IF_PARALLEL
#define DEBUG
using namespace seqan;

namespace SEQAN_NAMESPACE_MAIN
{    
	// ============================================================================
	// Tags, Classes, Enums
	// ============================================================================
	double timeCollectSeeds 	= 0;
	double timeCollectSeedsLoop = 0;
	double timeGardenerFind 	= 0;
	double timePutSeedsInMap 	= 0;
	double timeCSFreeSpace		= 0;
	long long cntCSFind 		= 0;
	long long cntCSNewNeedle	= 0;
	long long cntCSExistingNeedle	= 0;


	struct _MULTIPLE_WORKER;
	typedef Tag<_MULTIPLE_WORKER> MULTIPLE_WORKER; // tag for parallel execution
	
	struct _SINGLE_WORKER;
	typedef Tag<_SINGLE_WORKER> SINGLE_WORKER; // tag for serial execution
	
	/**
	 * Class that holds a hit between a needle and a haystack
	 *
	 **/
	template <typename TSpec, typename TPos, typename TId>
	class GardenerHit_
	{
	public:
		TId			hstId;			// haystack sequence id 
		TId			ndlSeqNo;		// needle sequence number
		TPos		hstkPos;		// begin in haystack 
		TPos		ndlPos;			// begin position of hit in needle
		TPos		diag;			// the diagonal
		TPos		score;			// the score
		TPos		hitLength;		// length of the hit
		
		bool operator==(const GardenerHit_<TSpec, TPos, TId>& b) const;
		bool operator!=(const GardenerHit_<TSpec, TPos, TId>& b) const;
		bool operator<(const GardenerHit_<TSpec, TPos, TId>& b) const;
		bool operator>(const GardenerHit_<TSpec, TPos, TId>& b) const;
		
		GardenerHit_(){}
		
		GardenerHit_(TId _hstId, TId _ndlSeqNo, TPos _hstkPos, TPos _ndlPos, TPos _diag, TPos _score, TPos _hitLength):
		hstId(_hstId),
		ndlSeqNo(_ndlSeqNo),
		hstkPos(_hstkPos),
		ndlPos(_ndlPos),
		diag(_diag),
		score(_score),
		hitLength(_hitLength)
		{}
		
		GardenerHit_(GardenerHit_ const &orig): 
		hstId(orig.hstId),
		ndlSeqNo(orig.ndlSeqNo),
		hstkPos(orig.hstkPos),
		ndlPos(orig.ndlPos),
		diag(orig.diag),
		score(orig.score),
		hitLength(orig.hitLength)
		{};
		
		GardenerHit_ & operator = (GardenerHit_ const &orig){
			hstId = orig.hstId;
			ndlSeqNo = orig.ndlSeqNo;
			hstkPos = orig.hstkPos;
			ndlPos = orig.ndlPos;
			diag = orig.diag;
			score = orig.score;
			hitLength = orig.hitLength;
			return *this;
		}
		
		~GardenerHit_(){}
		
		inline TId getHstId(){
			return hstId;
		}
		
		inline TId getNdlSeqNo(){
			return ndlSeqNo;
		}
		
		inline TPos getHstkPos(){
			return hstkPos;
		}
		
		inline TPos getNdlPos(){
			return ndlPos;
		}
		
		inline TPos getDiag(){
			return diag;
		}
		
		inline TPos getScore(){
			return score;
		}
		
		inline TPos getHitLength(){
			return hitLength;
		}
	};
	
	//____________________________________________________________________________
	
	/**
	 * specializations
	 **/
	struct GardenerUngapped_;
	typedef Tag<GardenerUngapped_> GardenerUngapped;
	
	struct GardenerUngappedSegmented_;
	typedef Tag<GardenerUngappedSegmented_> GardenerUngappedSegmented;
	
	//____________________________________________________________________________
	
	template <typename TId, typename TGardenerSpec>
	class Gardener
	{
	public:
		typedef __int64											TPos;
		typedef GardenerHit_<TGardenerSpec, TPos, TId>			TGardenerHit;
		typedef Map<TGardenerHit, Skiplist<> >					THitSet;
		typedef Pair<TId, THitSet* >							THitMapPair;
		typedef Map<THitMapPair, Skiplist< > >					THitMap;	
		typedef typename Iterator<THitSet, Standard>::Type		THitIterator;
		
		double timeQgramFind;
		double timeQgramMultiSeed;
		double timeCollectSeeds;
		double timeGardenerFind;
		double timePutSeedsInMap;

		THitMap hits; // containing for each  duplex sequence (TId) the list of detected
		
		Gardener<TId, TGardenerSpec>() :
				timeQgramFind(0),
				timeCollectSeeds(0),
				timeGardenerFind(0),
				timeQgramMultiSeed(0),
				timePutSeedsInMap(0) {}
		
		Gardener<TId, TGardenerSpec>(Gardener<TId, TGardenerSpec> const &orig): hits(orig.hits),
				timeQgramFind(0),
				timeCollectSeeds(0),
				timeGardenerFind(0),
				timeQgramMultiSeed(0),
				timePutSeedsInMap(0) {};
		
		Gardener<TId, TGardenerSpec> & operator = (Gardener<TId, TGardenerSpec> const &orig){
			hits = orig.hits;
			timeQgramFind 		= orig.timeQgramFind;
			timeGardenerFind 	= orig.timeGardenerFind;
			timeCollectSeeds	= orig.timeCollectSeeds;
			timeQgramMultiSeed	= orig.timeQgramMultiSeed;
			timePutSeedsInMap	= orig.timePutSeedsInMap;

			return *this;
		}
		
		~Gardener<TId, TGardenerSpec>()
		{
		}

	};

	// seed comparator
	// ATTENTION assumes same diagonal and size of seeds
	template <typename TSeed>
	struct LessRSeed : public ::std::binary_function < TSeed, TSeed, bool >
	{
		inline bool operator() (TSeed const &a, TSeed const &b) const 
		{
			if (getBeginDim0(a) > getBeginDim0(b)) return true;
			return false;
		}
	};
	
	// indicates if two seeds overlap
	// ATTENTION assumes same diagonal and size of seeds
	template <typename TSeed>
	bool isOverlapping(TSeed const & a, TSeed const & b){
		if (getBeginDim0(a) < getBeginDim0(b) and getBeginDim0(b) < getEndDim0(a)) return true;
		if (getBeginDim0(b) < getBeginDim0(a) and getBeginDim0(a) < getEndDim0(b)) return true;
		return false;
	}
	
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//											Meta Functions		                                                  //
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///.Metafunction.Spec.param.T.type:Class.Gardener
	
	template <typename TId, typename TGardenerSpec>
	struct Spec<Gardener<TId, TGardenerSpec> >
	{
		typedef TGardenerSpec Type;
	};
	
	//____________________________________________________________________________
	
	template <typename TId, typename TGardenerSpec>
	struct Iterator<Gardener<TId, TGardenerSpec> >
	{
		typedef typename Gardener<TId, TGardenerSpec>::THitIterator Type;
	};
	
	//____________________________________________________________________________
	
	template <typename TId, typename TGardenerSpec>
	struct Cargo<Gardener<TId, TGardenerSpec> >
	{
		typedef typename Gardener<TId, TGardenerSpec>::TGardenerHit Type;
	};
	
	
	//____________________________________________________________________________
	
	template <typename TId, typename TGardenerSpec>
	struct Id<Gardener<TId, TGardenerSpec> >
	{
		typedef typename Gardener<TId, TGardenerSpec>::TId Type;
	};
	
	//____________________________________________________________________________
	
	template <typename TId, typename TGardenerSpec>
	struct Value<Gardener<TId, TGardenerSpec> >
	{
		typedef typename Gardener<TId, TGardenerSpec>::THitMap Type;
	};
	
	//____________________________________________________________________________
	
	template <typename TId, typename TGardenerSpec>
	struct Position<Gardener<TId, TGardenerSpec> >
	{
		typedef typename Gardener<TId, TGardenerSpec>::TPos Type;
	};
	
	//____________________________________________________________________________
	
	template <typename TId, typename TGardenerSpec>
	inline void reset(Gardener<TId, TGardenerSpec> & me)
	{
		eraseAll(me);
	}
	
	//____________________________________________________________________________
	
	template <typename TId, typename TGardenerSpec>
	inline typename Value<Gardener<TId, TGardenerSpec> >::Type & hits(Gardener<TId, TGardenerSpec> & me)
	{
		return me.hits;
	}
	
	//____________________________________________________________________________
	
	template <typename TSpec, typename TPos, typename TId>
	bool GardenerHit_<TSpec, TPos, TId>::operator==(const GardenerHit_<TSpec, TPos, TId>& b) const {
		if (hstId != b.hstId) return false;
		if (ndlSeqNo != b.ndlSeqNo) return false;
		if (hitLength != b.hitLength) return false;		
		if (ndlPos != b.ndlPos) return false;
		if (hstkPos != b.hstkPos) return false;
		return true;
	}
	
	//____________________________________________________________________________
	
	template <typename TSpec, typename TPos, typename TId>
	bool GardenerHit_<TSpec, TPos, TId>::operator!=(const GardenerHit_<TSpec, TPos, TId>& b) const {
		return !(*this == b);
	}
	
	//____________________________________________________________________________
	
	template <typename TSpec, typename TPos, typename TId>
	bool GardenerHit_<TSpec, TPos, TId>::operator<(const GardenerHit_<TSpec, TPos, TId>& b) const {
		if (hstId < b.hstId) return true;
		if (hstId > b.hstId) return false;
		if (ndlSeqNo < b.ndlSeqNo) return true;
		if (ndlSeqNo > b.ndlSeqNo) return false;
		if (hstkPos < b.hstkPos) return true;
		if (hstkPos > b.hstkPos) return false;
		if (ndlPos < b.ndlPos) return true;
		if (ndlPos > b.ndlPos) return false;
		if (hitLength < b.hitLength) return true;
		return false;
	}
	
	//____________________________________________________________________________
	
	template <typename TSpec, typename TPos, typename TId>
	bool GardenerHit_<TSpec, TPos, TId>::operator>(const GardenerHit_<TSpec, TPos, TId>& b) const {
		return b<*this;
	}


	//____________________________________________________________________________
	
	/** 
	 * returns an iterator to the begin of the results for a given query id
	 */
	template<typename TSpec, typename TId >
	inline bool
	hasAnyHit(Gardener<TId, TSpec> &gardener, 
			  TId &queryid
			  ){ ENTER
		if (hasKey(gardener.hits, queryid))
			return true;
		else
			return false;
	}
	
	/** 
	 * returns an iterator to the begin of the results for a given query id
	 */
	template<typename TSpec, typename TId >
	inline typename Iterator<typename Gardener<TId, TSpec>::THitSet>::Type 
	harvestBegin(Gardener<TId, TSpec> &gardener, 
				 TId &queryid
				 ){ ENTER
		if (hasKey(gardener.hits, queryid))
			return begin(*cargo(gardener.hits, queryid));
		else
			return NULL;
	}
	
	//____________________________________________________________________________
	
	/** 
	 * returns an iterator to the end of the results for a given query id
	 */
	template<typename TSpec, typename TId >
	inline typename Iterator<typename Gardener<TId, TSpec>::THitSet>::Type 
	harvestEnd(Gardener<TId, TSpec> &gardener, 
			   TId &queryid
			   ){ ENTER
		if (hasKey(gardener.hits, queryid))
			return end(*cargo(gardener.hits, queryid));
		else
			return NULL;
	}
	
	//____________________________________________________________________________
	
	/** 
	 * returns the size of the result container
	 */
	template<typename TSpec, typename TId >
	inline typename Size<typename Gardener<TId, TSpec>::THitSet>::Type 
	length(Gardener<TId, TSpec> &gardener, 
		   TId &queryid
		   ){ ENTER
		if (hasKey(gardener.hits, queryid))
			return length(*cargo(gardener.hits, queryid));
		else
			return 0;
	}
	
	//____________________________________________________________________________
	
	/** 
	 * indicates if a specific query id is registered in the result array 
	 */
	template<typename TSpec, typename TId >
	inline bool hasHit(Gardener<TId, TSpec> &gardener, 
					   TId &queryid, 
					   typename Gardener<TId, TSpec>::TGardenerHit hit
					   ){ ENTER
		if (hasKey(gardener.hits, queryid) && hasKey(*cargo(gardener.hits,queryid), hit))
			return true;
		else {
			return false;
		}
	}
	

	
	//____________________________________________________________________________
	
	/** 
	 * erase all hits
	 */
	template<typename TSpec, typename TId >
	inline void eraseAll(Gardener<TId, TSpec>  &gardener){ ENTER
		typedef typename Value<Gardener<TId, TSpec> >::Type	THitMap;
		typedef typename Iterator<THitMap>::Type				THitMapIter;
		typedef typename Value<THitMap>::Type					THitMapPair;
		
		THitMap hitmap = hits(gardener);
		for (THitMapIter it = begin(hitmap);it != end(hitmap);++it){
			THitMapPair hm = *it;
			delete(hm.i2);
		}
		clear(gardener.hits);
	}	

	/**
	 * create runs of mismatches
	 * Triple contains <seed offset mismatch begin, seed offset mismatch end, total errors in sequence>
	 */
	template<typename THaystack, typename TNeedle, typename TSeed, typename TTripleSet>
	inline void	_fillRuns(THaystack const	&haystack,
						  TNeedle const		&needle,
						  TSeed				&seed,
						  TTripleSet		&runs
						  )
	{ ENTER
		SEQAN_CHECKPOINT
		typedef typename Value<TTripleSet>::Type	TRun;
		typedef typename Value<TRun, 1>::Type		TPos;
		TPos totalErrors = 0;
		TPos mismatchBegin = 0;
		TPos i = mismatchBegin;
		TPos seedLength = getEndDim0(seed)-getBeginDim0(seed);
		
		// append mismatch run starting at beginPosition
		while(i < seedLength && !isMatch(haystack[getBeginDim0(seed)+i], needle[getBeginDim1(seed)+i])) {
			++i;
			++totalErrors;
		}
		appendValue(runs, TRun(mismatchBegin, i, totalErrors));		
		
		// iterate over alignment and append mismatch runs
		while (i < seedLength) {
			// skip matches
			while(i < seedLength && isMatch(haystack[getBeginDim0(seed)+i], needle[getBeginDim1(seed)+i])) {
				++i;
			}
			mismatchBegin = i;
			// skip and count mismatches
			while(i < seedLength && !isMatch(haystack[getBeginDim0(seed)+i], needle[getBeginDim1(seed)+i])) {
				++i;
				++totalErrors;
			}
			appendValue(runs, TRun(mismatchBegin, i, totalErrors));			
		}
	#ifdef TRIPLEX_DEBUG	
		std::cout << infix(needle, getBeginDim1(seed), getEndDim1(seed)) << "\n";
		std::cout << infix(haystack, getBeginDim0(seed), getEndDim0(seed)) << "\n";
		for(unsigned l = 0; l < length(runs); ++l) {
			std::cout << getValue(runs,l).i1 << "  " << getValue(runs,l).i2 << "  " << getValue(runs,l).i3 << std::endl;
		}
	#endif		
	}


	//____________________________________________________________________________


	// checks the error rate of the fragment between end of left and start of right
	template<typename TPos, typename TError>
	inline bool
	_isEpsMatch(Triple<TPos, TPos, TPos> const &left,
				Triple<TPos, TPos, TPos> const &right,
				TError &errorRate
				) { ENTER
		SEQAN_CHECKPOINT
		// compute mismatches and length
		TPos errors = right.i3 - left.i3 - (right.i2 - right.i1);
		TPos length = right.i1 - left.i2;
		
		// check error rate
		return errors/(TError)(length) <= errorRate;
	}

	//____________________________________________________________________________


	// counts the number of matches
	template<typename TPos>
	inline TPos
	_countMatches(Triple<TPos, TPos, TPos> const &left,
				  Triple<TPos, TPos, TPos> const &right
				  ) { ENTER
		SEQAN_CHECKPOINT
		// compute mismatches and length
		TPos errors = right.i3 - left.i3 - (right.i2 - right.i1);
		TPos length = right.i1 - left.i2;
		
		// check error rate
		return length-errors;
	}


	//____________________________________________________________________________


	/**
	 * Identifies the longest epsilon match and adjusts the seed accordingly
	 */
	template<typename THaystack, typename TNeedle, typename TSeed, typename TPos, typename TError>
	inline TPos _longestEpsMatch(THaystack const	&haystack,
								 TNeedle const		&needle,
								 TSeed				&seed,
								 TPos const			&matchMinLength,
								 TError	const		&errorRate
								 ){ ENTER
		SEQAN_CHECKPOINT
		// Preprocessing: compute and store mismatch and lengths
		// A run is a triple of mismatch begin position, mismatch end position, 
		// and total number of errors in sequence from begin to end position of this mismatch.
		typedef std::vector<Triple<TPos, TPos, TPos> >	TRuns;
		typedef typename Iterator<TRuns, Rooted>::Type	TRunsIter;
		
		TRuns runs;
		_fillRuns(haystack, needle, seed, runs);
		
		// Identify longest eps match by iterating over combinations of left and right positions
		TRunsIter leftIt = begin(runs, Rooted());
		TRunsIter rightIt = end(runs, Rooted());
		--rightIt;
		
		TPos beginOffset = 0;
		TPos endOffset = 0;
		TPos minLength = matchMinLength - 1;
		TPos matches = 0;
		
		// border up to which the begin of the seed is able to fulfil the minimum length constraint
		TPos seedBeginBorder = max(TPos(0),(*rightIt).i1 - matchMinLength);
		TPos seedEndBorder = (*leftIt).i2 + matchMinLength;
		
		while ((*leftIt).i2 + minLength < (*rightIt).i1 && (*leftIt).i2 <= seedBeginBorder) {
			while ((*leftIt).i2 + minLength < (*rightIt).i1 && (*rightIt).i2 >= seedEndBorder) {
				if(_isEpsMatch(*leftIt, *rightIt, errorRate)) {
					beginOffset = (*leftIt).i2;
					endOffset = (*rightIt).i1; 
					minLength = endOffset - beginOffset;
					matches = _countMatches(*leftIt, *rightIt);
					break;
				}
				--rightIt;
			}
			rightIt = end(runs);
			--rightIt;
			++leftIt;
		}
		
		// adjust the seed borders
		setEndDim0(seed, getBeginDim0(seed) + endOffset);					
		setEndDim1(seed, getBeginDim1(seed) + endOffset);		
		setBeginDim0(seed, getBeginDim0(seed) + beginOffset);
		setBeginDim1(seed, getBeginDim1(seed) + beginOffset);
		return matches;
	}

	//____________________________________________________________________________
	
	/**	
	 * Clips the ends of seeds and add to hitlist
	 */
	template <
	typename THitSet,
	typename TId,
	typename TSeqNo,
	typename TDiag,
	typename THaystack,
	typename TNeedle,
	typename TSeedSet,
	typename TError,
	typename TSize
	>
	inline void _clipSeedEndsAndAdd(THitSet		&hitSet,
									TId			const &queryid,
									TSeqNo		const &seqno,
									TDiag		const &diag,
									THaystack	const &haystack,
									TNeedle		const &needle,
									TSeedSet	&seedset,
									TError		&errorRate,
									TSize		&minLength,
									SINGLE_WORKER const & 
									){ ENTER
		typedef typename Iterator<TSeedSet, Rooted>::Type	TIter;
		typedef typename Value<TSeedSet>::Type				TSeed;
		typedef typename Position<TSeed>::Type				TPos;
		typedef typename Value<THitSet>::Type				THit;
		
		// cycle through each seedset at a time
		TIter it = begin(seedset, Rooted());
		TIter itEnd = end(seedset, Rooted());
		
		while(it != itEnd){		
			// cut ends to obtain longest epsilon-match that contains the seed alignment
			TPos matches = _longestEpsMatch(haystack, needle, *it, minLength, errorRate);
			
			// size constaint fulfilled ?
			TPos seedLength = getEndDim0(*it)-getBeginDim0(*it);	
			if (seedLength  >= (TPos)minLength){
				
				// create a new hit and append it to the gardeners hit list
				THit hit(queryid,
						 seqno,					// needle seq. number            
						 getBeginDim0(*it),      // begin in haystack      
						 getBeginDim1(*it),		// needle position
						 diag,					// the diagonal
						 matches,
						 seedLength
						 );
				
				// append the hit to the finders hit list if not already contained
				if (!hasKey(hitSet, hit)){
					add(hitSet, hit);
#ifdef TRIPLEX_DEBUG
					::std::cout << "new hit :" << ::std::endl;
					_printHit(hit);
#endif			
				}
			}
			++it;
		}	
	}
	

	//____________________________________________________________________________

	/**	
	 * Clips the ends of seeds and add to hitlist
	 */
	template <
	typename THitSet,
	typename TId,
	typename TSeqNo,
	typename TDiag,
	typename THaystack,
	typename TNeedle,
	typename TSeedSet,
	typename TError,
	typename TSize
	>
	inline void _clipSeedEndsAndAdd(THitSet		&hitSet,
									TId			const &queryid,
									TSeqNo		const &seqno,
									TDiag		const &diag,
									THaystack	const &haystack,
									TNeedle		const &needle,
									TSeedSet	&seedset,
									TError		&errorRate,
									TSize		&minLength,
									MULTIPLE_WORKER const & 
									){
		typedef typename Iterator<TSeedSet, Rooted>::Type	TIter;
		typedef typename Value<TSeedSet>::Type				TSeed;
		typedef typename Position<TSeed>::Type				TPos;
		typedef typename Value<THitSet>::Type				THit;
		
		// cycle through each seed at a time
		TIter it = begin(seedset, Rooted());
		TIter itEnd = end(seedset, Rooted());
		
		while(it != itEnd){		
			// cut ends to obtain longest epsilon-match that contains the seed alignment
			TPos matches = _longestEpsMatch(haystack, needle, *it, minLength, errorRate);
			
			// size constaint fulfilled ?
			TPos seedLength = getEndDim0(*it)-getBeginDim0(*it);	
			if (seedLength  >= (TPos)minLength){
				
				// create a new hit and append it to the gardeners hit list
				THit hit(queryid,
						 seqno,					// needle seq. number            
						 getBeginDim0(*it),      // begin in haystack      
						 getBeginDim1(*it),		// needle position
						 diag,					// the diagonal
						 matches,
						 seedLength
						 );
				
				// append the hit 
				appendValue(hitSet, hit);
			}
			++it;
		}	
	}
	

	//____________________________________________________________________________
	
	/**	
	 * Extends all seeds according to the scoring schema and a X-dropoff 
	 */
	template <
	typename THitSet,
	typename THaystack,
	typename TIndex,
	typename TSpec,
	typename TMap,
	typename TSize,
	typename TScore,
	typename TId
	>
	inline void _extendSeedlings(THitSet					&hitSet,
								 Finder<THaystack, TSpec>	&finder,
								 Pattern<TIndex, TSpec> const &pattern,
								 TMap						&seqmap,
								 Score<TScore, Simple> const &scoreMatrix, 
								 TSize const				&minLength,
								 TScore	const				&scoreDropOff, 
								 TId						&queryid
								 ){
		typedef typename Iterator<TMap>::Type				TMapIter;
		typedef Finder<THaystack, TSpec>					TFinder;
		typedef typename Value<TMap>::Type					TMapPair;
		typedef typename Value<TMapPair,1>::Type			TSeqNo;
		typedef typename Value<TMapPair,2>::Type			TDiagMapPointer;
		typedef typename Value<TDiagMapPointer>::Type		TDiagMap;
		typedef typename Iterator<TDiagMap>::Type			TDiagMapIter;
		typedef typename Value<TDiagMap>::Type				TDiagMapPair;		
		typedef typename Value<TDiagMapPair,1>::Type		TDiag;
		typedef typename Value<TDiagMapPair,2>::Type		TSeedSetPointer;
		typedef typename Value<TSeedSetPointer>::Type		TSeedSet;
		typedef typename Value<TSeedSet>::Type				TSeed;
		typedef typename Position<TSeedSet>::Type			TPos;
		typedef Dequeue<TSeed>								TSeedList;// @TODO known memory leak in seqan string http://trac.mi.fu-berlin.de/seqan/ticket/364
		typedef typename Host<TFinder>::Type				THost;
		typedef typename Value<THitSet>::Type				THit;
		
		// process all needles with entries in the map
		TMapIter itsme = end(seqmap);
		for (TMapIter itsmb = begin(seqmap); itsmb != itsme; ++itsmb){
			TSeqNo &seqno = key(*itsmb);
			// check that seqno is valid
			if (seqno < (TId)countSequences(needle(pattern))){
				TDiagMapPointer diagmapP = cargo(*itsmb);
				TDiagMap diagmap = *diagmapP;
				// process each diagonal at a time
				TDiagMapIter itdme = end(diagmap) ;
				for (TDiagMapIter itdmb = begin(diagmap); itdmb != itdme; ++itdmb){
					TDiag &diag = key(*itdmb);
					// get the seed set
					TSeedSetPointer seedsetP = cargo(*itdmb);
					TSeedSet seedset = *seedsetP;
					THost tmp = host(finder);
					// extend all seeds by first overlapping with succeeding seeds
					TSeedList newset;
					while (!empty(seedset)){
						TSeed seed = front(seedset);
						popFront(seedset);
#ifdef TRIPLEX_DEBUG
						::std::cout << "seed_1:" << getBeginDim0(seed) << ":" << getEndDim0(seed) << " " << getBeginDim1(seed) << ":" << getEndDim1(seed) << ::std::endl;
						cout << seed << endl;
#endif
						
						while (!empty(seedset) and isOverlapping(seed, front(seedset))){
							TSeed seed2 = front(seedset);
							setBeginDim0(seed, min(getBeginDim0(seed),getBeginDim0(seed2)));
							setBeginDim1(seed, min(getBeginDim1(seed),getBeginDim1(seed2)));
							setEndDim0(seed, max(getEndDim0(seed),getEndDim0(seed2)));
							setEndDim1(seed, max(getEndDim1(seed),getEndDim1(seed2)));
							popFront(seedset);
						}
						
#ifdef TRIPLEX_DEBUG
						::std::cout << "seed_2:" << getBeginDim0(seed) << ":" << getEndDim0(seed) << " " << getBeginDim1(seed) << ":" << getEndDim1(seed) << ::std::endl;
						cout << seed << endl;
#endif
						// extend seed to both sides as far as possible
						extendSeed(seed, tmp, getSequenceByNo(seqno,needle(pattern)), EXTEND_BOTH, scoreMatrix, scoreDropOff, UnGappedXDrop());

#ifdef TRIPLEX_DEBUG
						::std::cout << "seed_3:" << getBeginDim0(seed) << ":" << getEndDim0(seed) << " " << getBeginDim1(seed) << ":" << getEndDim1(seed) << ::std::endl;
						cout << seed << endl;
#endif
						
						// merge overlapping windows
						if (getEndDim0(seed)-getBeginDim0(seed) >=  (TPos)minLength){
							if (!empty(newset) && isOverlapping(seed, back(newset))){
								setBeginDim0(seed, min(getBeginDim0(seed),getBeginDim0(back(newset))));
								setBeginDim1(seed, min(getBeginDim1(seed),getBeginDim1(back(newset))));
								setEndDim0(seed, max(getEndDim0(seed),getEndDim0(back(newset))));
								setEndDim1(seed, max(getEndDim1(seed),getEndDim1(back(newset))));
							} else {
								pushBack(newset,seed);
							}
						}
					}
					
					// add all now non-overlapping windows to hitlist
					while(!empty(newset)){
//						std::cout << "Diagonal: " << diag << std::endl;
						// create a new hit and append it to the gardeners hit list
						THit hit(queryid,				// hstId
								 seqno,					// needle seq. number            
								 getBeginDim0(front(newset)),      // begin in haystack (hstkPos)
								 getBeginDim1(front(newset)),		// needle position 	(ndlPos)
								 diag,					// the diagonal
								 0,
								 getEndDim0(front(newset))-getBeginDim0(front(newset))
								 );

#ifdef TRIPLEX_DEBUG
						::std::cout << "newset:" << length(newset) << " hitSet:" << length(hitSet) << ::std::endl;
#endif
						
						// append the hit 
						add(hitSet, hit);
						popFront(newset);
#ifdef TRIPLEX_DEBUG
						::std::cout << "newset:" << length(newset) << " hitSet:" << length(hitSet) << ::std::endl;
#endif
					}
					
//					// find the longest match conform with error rate
//					_clipSeedEndsAndAdd(hitSet, queryid, seqno, diag, host(finder), getSequenceByNo(seqno,needle(pattern)), newset, errorRate, minLength, SINGLE_WORKER() );
					
				} //diagmap
			} else {
				::std::cerr << "Sequence no " << seqno << " exceeds index " << ::std::endl;
			}
		}
	}
	
	//____________________________________________________________________________

	/**
	 * adjust the seed position with respect to the host sequence
	 */
	template <
	typename TSeedSet,
	typename TSequence,
	typename TSpec,
	typename TSize
	>
	inline void _adjustSeeds(TSeedSet					&seedset, 
							 Segment<TSequence, TSpec>	&sequence, 
							 TSize						dim
							 ){ ENTER
		typedef typename Iterator<TSeedSet>::Type					TSeedSetIter;
		typedef typename Position<Segment<TSequence, TSpec> >::Type	TPos;
		
		TSeedSetIter itsE = end(seedset);
		TPos offset = beginPosition(sequence); 
		for (TSeedSetIter it = begin(seedset); it != itsE; ++it){
			if (dim){
				setBeginDim1(*it, getBeginDim1(*it)+offset);
				setEndDim1(*it, getEndDim1(*it)+offset);
			} else {
				setBeginDim1(*it, getBeginDim0(*it)+offset);
				setEndDim0(*it, getEndDim0(*it)+offset);
			}
		}
	}

	//____________________________________________________________________________
	/**
	 * copy seed passing the qgram lemma
	 */
	template <
	typename TMap,
	typename TId,
	typename TDiag,
	typename TSet
	>
	inline void _putSeedsInMap(TMap			&seqmap,
							   TId	const	&seqNo,
							   TDiag const	&diag,
							   TSet		&posSet
							   ){
        SEQAN_PROTIMESTART(time_putseeds);

		typedef typename Value<TMap>::Type				TMapPair;
		typedef typename Cargo<TMapPair>::Type			TDiagMapPointer;
		typedef typename Value<TDiagMapPointer>::Type 	TDiagMap;
		typedef typename Value<TDiagMap>::Type			TDiagMapPair;
		typedef typename Cargo<TDiagMapPair>::Type		TSeedSetPointer;
		typedef typename Value<TSeedSetPointer>::Type	TSeedSet;
		typedef typename Value<TSeedSet>::Type			TSeed;
		typedef typename Iterator<TSet>::Type			TPosIter;
		typedef typename Value<TSet>::Type				TPos;
		
#ifdef TRIPLEX_DEBUG			
		::std::cout << "add new window:" << length(posSet) << "-" << front(posSet) << " " << back(posSet) << ::std::endl;
#endif	
		// new needle sequence that has no entries yet -- add new needle, and seedset corresponding to diagonal
		if ( ! hasKey(seqmap, seqNo)){		
			// create new needle map
			{
				TDiagMapPointer diagMapPointer = new TDiagMap;
				insert(seqmap, seqNo, diagMapPointer);
			}
			// create new seedset
			{
				TDiagMapPointer tmp_diagMapPointer = cargo(seqmap, seqNo);
				TSeedSetPointer seedSetPointer = new TSeedSet;
				TSeed seed(diag+front(posSet), front(posSet), back(posSet)-front(posSet)+1);
				pushBack(*seedSetPointer,seed);
				insert(*tmp_diagMapPointer, diag, seedSetPointer);
			}
		}
		// needle sequence is known
		else { 
			TDiagMap* diagMapPointer = cargo(seqmap, seqNo);
			// but diagonal index is new
			if (!hasKey(*diagMapPointer, diag)){
				// create new seedset for this diagonal
				TSeedSet* seedSetPointer = new TSeedSet;
				TSeed seed2(diag+front(posSet), front(posSet), back(posSet)-front(posSet)+1);
				pushBack(*seedSetPointer,seed2);
				insert(*diagMapPointer, diag, seedSetPointer);
			} 
			else { // diagonal index is known -- push seed on heap
				TSeedSet* seedSetPointer = cargo(*diagMapPointer, diag);
				TSeed seed3(diag+front(posSet), front(posSet), back(posSet)-front(posSet)+1);
				// check if window simply extends previous one
#ifdef TRIPLEX_DEBUG			
				::std::cout << "extend? " << getBeginDim0(back(*seedSetPointer)) << "=" << getBeginDim0(seed3) << " " << getEndDim0(back(*seedSetPointer)) << " " << getEndDim0(seed3) << ::std::endl;
#endif	
				if (getBeginDim0(back(*seedSetPointer)) <= getBeginDim0(seed3) && getBeginDim0(seed3) <= getEndDim0(back(*seedSetPointer))){
#ifdef TRIPLEX_DEBUG			
					::std::cout << "was: " << getBeginDim0(back(*seedSetPointer)) << "-" << getEndDim0(back(*seedSetPointer)) << " " << getBeginDim1(back(*seedSetPointer)) << "-" << getEndDim1(back(*seedSetPointer)) <<  ::std::endl;
#endif	
					setEndDim0(back(*seedSetPointer), getEndDim0(seed3));
					setEndDim1(back(*seedSetPointer), getEndDim1(seed3));					
#ifdef TRIPLEX_DEBUG			
					::std::cout << "now: " << getBeginDim0(back(*seedSetPointer)) << "-" << getEndDim0(back(*seedSetPointer)) << " " << getBeginDim1(back(*seedSetPointer)) << "-" << getEndDim1(back(*seedSetPointer)) <<  ::std::endl;
#endif	
					
				} else { // add
					pushBack(*seedSetPointer,seed3);
				}
			}
		}	
        timePutSeedsInMap += SEQAN_PROTIMEDIFF(time_putseeds);
	}
	
	//____________________________________________________________________________
	/**
	 * get all the hits between needles and haystack
	 * TMap seqmap will be filled using the new operator 
	 * freeing the memory is up to the caller
	 */
	template <
	typename THaystack,
	typename TSpec,
	typename TIndex,
	typename TPos,
	typename TMap
	>
	inline void _collectSeeds(Finder<THaystack, QGramsLookup<TSpec> >		&finder,
							  Pattern<TIndex,  QGramsLookup<TSpec> > const	&pattern,
							  TPos const									&seedsThreshold,
							  TPos const									&minLength,
							  TMap											&seqmap
							  ){
        SEQAN_PROTIMESTART(time_collectseeds);
		typedef typename Value<TMap>::Type				TMapPair;
		typedef typename Key<TMapPair>::Type			TId;
		typedef typename Cargo<TMapPair>::Type			TDMPointer;
		typedef typename Value<TDMPointer>::Type		TDM;
		typedef typename Value<TDM>::Type				TDMPair;
		typedef typename Key<TDMPair>::Type				TDiag;
		
		// tmpSeqMap structure constains list of needle position rather than seeds
		typedef TPos									TSeed;
		typedef Dequeue<TSeed>							TSeedSet;
		typedef Pair<TDiag, TSeedSet*>					TDiagMapPair;
		typedef Map<TDiagMapPair, Skiplist< > >			TDiagMap;
		typedef Pair<TId, TDiagMap*>					TSeqMapPair;
		typedef Map<TSeqMapPair, Skiplist< > >			TSeqMap;
		typedef typename Iterator<TDiagMap>::Type		TDiagMapIter;
		typedef typename Iterator<TSeqMap>::Type		TIterM;
		typedef typename Iterator<TDiagMap>::Type		TIterD;
		typedef typename Iterator<TSeedSet>::Type		TIterS;		
		
		TSeqMap tmpSeqmap;
		while (find(finder, pattern)) {
			cntCSFind++;
	        SEQAN_PROTIMESTART(time_collectseeds_loop);
#ifdef TRIPLEX_DEBUG
//	        ::std::cout << "\n_collectSeeds(): new seed found: ******************************\n";
//			::std::cout << "Q (inv = tfo):" << infix(finder) << ::std::endl; // The Infix of the match in the haystack.
//			::std::cout << "T (inv = tts):" << infix(pattern, *finder.curHit) << ::std::endl;
//			::std::cout << "H:" << (*finder.curHit).hstkPos << "-N" << (*finder.curHit).ndlSeqNo << ":P" << (*finder.curHit).ndlPos << ":D" << (*finder.curHit).diag << ::std::endl;
            //TODO@barni remove
            ::std::cout << "Q:" << infix(finder) << ::std::endl;
            ::std::cout << "T:" << infix(pattern, *finder.curHit) << ::std::endl;
            ::std::cout << "H:" << (*finder.curHit).hstkPos << "-N" << (*finder.curHit).ndlSeqNo << ":P" << (*finder.curHit).ndlPos << ":D" << (*finder.curHit).diag << ::std::endl;
            ::std::cout << "tmpseqmanSize:" << length(tmpSeqmap) << " ,KeyKnown:" << hasKey(tmpSeqmap, (*finder.curHit).ndlSeqNo) << ::std::endl;
            ::std::cout << "seqmanSize:" << length(seqmap) << " ,KeyKnown:" << hasKey(seqmap, (*finder.curHit).ndlSeqNo) << ::std::endl;
#endif			
			// new needle sequence that has no entries yet -- add new needle, and seedset corresponding to diagonal
			if ( ! hasKey(tmpSeqmap, (*finder.curHit).ndlSeqNo)){
				cntCSNewNeedle ++;
				// create new needle map
				{
					TDiagMap* diagMapPointer = new TDiagMap;
					insert(tmpSeqmap,(*finder.curHit).ndlSeqNo, diagMapPointer);
				}
				// create new seedset
				{
					TDiagMap* tmp_diagMapPointer = cargo(tmpSeqmap, (*finder.curHit).ndlSeqNo);
					TSeedSet* seedSetPointer = new TSeedSet;
					pushBack(*seedSetPointer,(*finder.curHit).ndlPos);
					if (length(*seedSetPointer) >= seedsThreshold){
						_putSeedsInMap(seqmap, (*finder.curHit).ndlSeqNo, (*finder.curHit).diag, *seedSetPointer);
					}
					
					insert(*tmp_diagMapPointer, (*finder.curHit).diag, seedSetPointer);
				}
			}
			// needle sequence is known
			else {
				TDiagMap* diagMapPointer = cargo(tmpSeqmap, (*finder.curHit).ndlSeqNo);
				// but diagonal index is new
				if (!hasKey(*diagMapPointer, (*finder.curHit).diag)){
					// create new seedset for this diagonal
					TSeedSet* seedSetPointer = new TSeedSet;
					pushBack(*seedSetPointer,(*finder.curHit).ndlPos);
					if (length(*seedSetPointer) >= seedsThreshold){
						_putSeedsInMap(seqmap, (*finder.curHit).ndlSeqNo, (*finder.curHit).diag, *seedSetPointer);
					}					
					insert(*diagMapPointer, (*finder.curHit).diag, seedSetPointer);
					// diagonal index is known -- push seed on heap
				} 
			 	else { 
					TSeedSet* seedSetPointer = cargo(*diagMapPointer, (*finder.curHit).diag);
					// remove positions outside the window
					while( !empty(*seedSetPointer) && front(*seedSetPointer)+minLength < (*finder.curHit).ndlPos+weight(pattern.shape)){
#ifdef TRIPLEX_DEBUG
						::std::cout << "Poping " << (*finder.curHit).ndlSeqNo << " 1st:" << front(*seedSetPointer) << " cur:" << (*finder.curHit).ndlPos << " diag:" << (*finder.curHit).diag << " " << (*finder.curHit).ndlPos+weight(pattern.shape)-minLength << " l" << length(*seedSetPointer) <<::std::endl;
#endif
						popFront(*seedSetPointer);
					}
					pushBack(*seedSetPointer,(*finder.curHit).ndlPos);
					if (length(*seedSetPointer) >= seedsThreshold){
						_putSeedsInMap(seqmap, (*finder.curHit).ndlSeqNo, (*finder.curHit).diag, *seedSetPointer);
					}
				}
			}
			timeCollectSeedsLoop += SEQAN_PROTIMEDIFF(time_collectseeds_loop);
		}

		// housekeeping
		// empty tmp seqmap properly
		SEQAN_PROTIMESTART(time_delete);
		for (TIterM sit=begin(tmpSeqmap); sit != end(tmpSeqmap); ++sit){
			TDiagMap* diagmapPointer = (*sit).i2;
			for (TIterD dit=begin(*diagmapPointer); dit != end(*diagmapPointer); ++dit){
				delete (*dit).i2;
			}
			delete diagmapPointer;
		}
		timeCSFreeSpace 	+= SEQAN_PROTIMEDIFF(time_delete);
		timeCollectSeeds 	+= SEQAN_PROTIMEDIFF(time_collectseeds);
	}

		
	//____________________________________________________________________________
	
	/**	
	 * Performs a q-gram search of the needle set and the haystack and merges the resulting seeds 
	 * subject to a scoring scheme, and an error rate.
	 * Reports hits subject to a minimum length constraint. 
	 * 
	 * Important note 1: The gardener is an ungapped search algorithm
	 * Important note 2: The algorithm is designed so it can be called using multiple haystacks 
	 * but the same set of needles (needles shared in threads)
	 * Important note 3: Gardener should be used with relatively short needles and haystacks only, since its
	 * memory consumption is O(|haystack| * sum(|needle|)) in worst case
	 */
	template <
	typename THitSet,
	typename THaystack,
	typename TSpec,
	typename TIndex,
	typename TError,
	typename TPos,
	typename TDrop,
	typename TId
	>
	inline bool _find(THitSet							&hitSet,
					  Finder<THaystack,  TSpec >		&finder,
					  Pattern<TIndex,  TSpec > const	&pattern,
					  TError const						&errorRate,
					  TPos const						&minLength,
					  TPos const						&seedsThreshold, // q-gram lemma
					  TDrop const						&xDrop,
					  TId								&queriyid
					  ){
        SEQAN_PROTIMESTART(time_find);
		// used datastructure: Map < NeedleSeqNo, < Map < diagonal, SeedSet > > 
		// for each needle all q-gram hits are stored according to the diagonal they reside in
		typedef int											TScore;
		typedef Seed<Simple, DefaultSeedConfig>				TSeed;
		//		typedef Seed<TPos, SimpleSeed>						TSeed;		
		//		typedef LessRSeed<TSeed>							TLess;
		//		typedef PriorityType<TSeed, TLess, PriorityHeap>	TSeedSet;
		typedef Dequeue<TSeed>								TSeedSet;
		typedef typename MakeSigned_<TPos>::Type			TDiag;
		typedef Pair<TDiag, TSeedSet*>						TDiagMapPair;
		typedef Map<TDiagMapPair, Skiplist< > >				TDiagMap;
		typedef Pair<TId, TDiagMap*>						TSeqMapPair;
		typedef Map<TSeqMapPair, Skiplist< > >				TSeqMap;
		typedef typename Iterator<TSeqMap>::Type			TSeqMapIter;
		typedef typename Iterator<TDiagMap>::Type			TDiagMapIter;
		
		// run gardener on first call
		if (empty(finder) ){
			TSeqMap seqmap;
			// get all maxed seeds for any needle in the haystack (flanked by mismatches)
			_collectSeeds(finder, pattern, seedsThreshold, minLength, seqmap);
			
			// define a scoring scheme
			TScore match = 1;
			TScore mismatch = (TScore)_max((TScore) (-1.0/(errorRate+0.00000001)) + 1, -(TScore)length(host(finder)));
			Score<TScore> scoreMatrix(match, mismatch, std::numeric_limits<int>::max());
			TScore scoreDropOff = (TScore) _max((TScore) xDrop * (-mismatch), minValue<TScore>()+1);
			// extend seeds
			_extendSeedlings(hitSet, finder, pattern, seqmap, scoreMatrix, minLength, scoreDropOff, queriyid);
			
			// housekeeping
			// free memory from gram-hitlist
			clear(finder.hits);
			// empty seqmap properly
			for (TSeqMapIter sit=begin(seqmap); sit != end(seqmap); ++sit){
				TDiagMap* diagmapPointer = (*sit).i2;
				for (TDiagMapIter dit=begin(*diagmapPointer); dit != end(*diagmapPointer); ++dit){
					delete (*dit).i2;
				}
				delete diagmapPointer;
			}
			timeGardenerFind += SEQAN_PROTIMEDIFF(time_find);
			return true;
		} else
			timeGardenerFind += SEQAN_PROTIMEDIFF(time_find);
			return false;
	}
	
	//____________________________________________________________________________

	/** 
	 * prints all hits for a query
	 */
	template< 
	typename TSpec, 
	typename TIndex,		// index 
	typename TPatternSpec,	// pattern spec
	typename TQuerySet,		// query set (needle)
	typename TId 
	>	
	inline void _printHits(Gardener<TId, TSpec>					&gardener, 
						   Pattern<TIndex, TPatternSpec> const	&pattern,
						   TQuerySet							&queries,
						   TId									&queryid
						   ){ ENTER
		typedef Gardener<TId, TSpec>						TGardener;
		typedef typename Iterator<TGardener >::Type			TIter;
		typedef typename Value<TGardener>::Type				THitMap;
		typedef typename Cargo<TGardener>::Type				THit;
		
		if (hasAnyHit(gardener, queryid)){
			for (TIter it = harvestBegin(gardener,queryid); it != harvestEnd(gardener, queryid); ++it){
				THit hit = *it;
				::std::cout << hit.ndlSeqNo << "("<< hit.ndlPos <<":" << (hit.ndlPos+hit.hitLength) << ")\t";
				::std::cout << hit.hstId << "("<< hit.hstkPos <<":" << (hit.hstkPos+hit.hitLength) << ")\t";
				::std::cout << "d:" << hit.diag << "\t";
				std::cout << infix(getSequenceByNo(hit.ndlSeqNo, needle(pattern)), hit.ndlPos, hit.ndlPos+hit.hitLength) << "\t";
				::std::cout << hit.ndlSeqNo << "("<< beginPosition(getSequenceByNo(hit.ndlSeqNo, needle(pattern)))+hit.ndlPos;
				::std::cout <<":" << beginPosition(getSequenceByNo(hit.ndlSeqNo, needle(pattern)))+(hit.ndlPos+hit.hitLength) << ")\t";
				std::cout << infix(queries[hit.hstId], hit.hstkPos, hit.hstkPos+hit.hitLength) << "\t";
				::std::cout << hit.hstId << "("<< beginPosition(queries[hit.hstId])+hit.hstkPos;
				::std::cout <<":" << beginPosition(queries[hit.hstId])+(hit.hstkPos+hit.hitLength) << ")\t";
				::std::cout << ::std::endl;
			}
		} else{
			::std::cout << "No Hits found" << ::std::endl;
		}
	}
	
	template< 
	typename THit 
	>	
	inline void _printHit(THit	&hit){ ENTER
		::std::cout << hit.ndlSeqNo << "("<< hit.ndlPos <<":" << (hit.ndlPos+hit.hitLength) << ")\t";
		::std::cout << hit.hstId << "("<< hit.hstkPos <<":" << (hit.hstkPos+hit.hitLength) << ")\t";
		::std::cout << "d:" << hit.diag << "\t";
		::std::cout << ::std::endl;	
	}
	//////////////////////////////////////////////////////////////////////////////


	//____________________________________________________________________________
	// entry functions

	
	/** 
	 * start gardening by planting from scratch preparing patterns (needles) first
	 */
	template< 
	typename TMotif,		// target set (pattern)
	typename TQuerySet,		// query set (needle)
	typename TError,		// error rate
	typename TSize,			// minimum hit size
	typename TDrop,			// xdrop
	typename TShape,		// shape
	typename TSpec,			// specialization
	typename TId,			// sequence id
	typename TWORKER		// Worker Tag
	>
	void plant(Gardener<TId, TSpec>		&gardener,
			   StringSet<TMotif> const	&targets,
			   TQuerySet				&queries,
			   TError const				&errorRate,
			   TSize const				&minLength,
			   TDrop const				&xDrop,
			   TShape const				&shape,
			   TWORKER
			   ){ ENTER
		typedef StringSet<TMotif>													TTargetSet;
		typedef Index<TTargetSet, IndexQGram<TShape, OpenAddressing> >				TQGramIndex;
		// configure Pattern
		typedef typename Value<TQuerySet>::Type										TSequence;
		typedef Finder<TSequence, QGramsLookup< TShape, Standard_QGramsLookup> >	TFinder;
		typedef typename Infix<typename GetSequenceByNo<TQGramIndex const>::Type >::Type	TInfix;
		typedef typename Iterator<TQuerySet>::Type									TQueryIter;
		
		typedef typename Fibre<TQGramIndex, QGramCounts>::Type						TSA;
		typedef typename Iterator<TSA, Standard>::Type								TSAIter;
		typedef typename Position<TSequence>::Type									TPosition;
		typedef typename Gardener<TId, TSpec>::THitSet								THits;
		
		// create index
		TQGramIndex index_qgram(targets);
		resize(indexShape(index_qgram), weight(shape));
		
		// create pattern
		Pattern<TQGramIndex, QGramsLookup< TShape, Standard_QGramsLookup> > pattern(index_qgram,shape);
		plant(gardener, pattern, queries, errorRate, minLength, xDrop, TWORKER() );
	}

	/**
	 * Binary search in sorted haystack segment map / vector.
	 * Returns the index (sequence number) of a duplex fiber in the haystack, given its starting
	 * position (bPos) in the merged haystack.
	 */
	template< typename TSegment	>
	int getHaystackFiberSeqNo(int bPos, TSegment segmentMap) {
		int l = 0;
		int r = segmentMap.size() - 1;
		int m;

		while (l <= r) {
			m = (l + r) / 2;
			if (bPos == segmentMap[m]) {
				break;
			}
			if (bPos < segmentMap[m]) {
				r = m - 1;
			}
			else {
				l = m + 1;
			}
		}

		if (segmentMap[m] > bPos) {
			m--;
		}

#ifdef DEBUG
		cout << "getHaystackFiberSeqNo: " << m << "(from " << segmentMap.size() << ") for bPos: " << bPos << endl;
#endif
		return m;
	}

	template<
	typename TSeedHashMap,
	typename TSeed,
	typename TSeedList,
	typename THaystackFiber,
	typename TQuery,
	typename TError>
	void extendSimpleSeedRevised(
			TSeedHashMap		  &addedSeedHashes,
			TSeed 			const &seed,
			TSeedList			  &extendedSeeds,
			THaystackFiber 	const &fiber,
			TQuery 			const &query,
			TError  		const &errorRate,
			int				const &k,
			int				const &minLength,
			int				const &consErrors) {
		extendedSeeds.clear();
		int posFiber;
		int posQuery;
		int mismatches;
		int consecutiveMismatches;
		int bDim0;
		int bDim1;
		int eDim0;
		int eDim1;

		int localK = 1000; // should be inf but this is enough

		bool isLeftZeroIncluded = false;
		bool isRightEndIncluded = false;

		// keep list of where mismatches occur, left and right from seed
		std::vector<int> lMmOffsets;
		std::vector<int> rMmOffsets;

		// go left
		mismatches = 0;
		consecutiveMismatches = 0;
		posFiber = getBeginDim0(seed);
		posQuery = getBeginDim1(seed);
		while (posFiber >= 0 && posQuery >= 0 && mismatches <= localK + 1
				&& consecutiveMismatches <= consErrors) {
			if (fiber[posFiber] != query[posQuery]) {
				lMmOffsets.push_back(getBeginDim0(seed) - posFiber);
				++mismatches;
				if (posFiber == 0 || posQuery == 0) {
					isLeftZeroIncluded = true;
				}
				++consecutiveMismatches;
			}
			else {
				consecutiveMismatches = 0;
			}
			localK = std::max((int)floor((getEndDiagonal(seed) - posFiber + 1) * errorRate), k);
			--posFiber;
			--posQuery;
		}

		// add corresponding end points if required
		if ((posFiber == -1 || posQuery == -1) && !isLeftZeroIncluded) {
			lMmOffsets.push_back(std::min((int)(getBeginDim0(seed)),(int)(getBeginDim1(seed))));
		}

		// go right
		mismatches = 0;
		consecutiveMismatches = 0;
		posFiber = getEndDim0(seed);
		posQuery = getEndDim1(seed);
		int endFiber = length(fiber);
		int endQuery = length(query);
		// only go until the second last element, the last one is added later if needed (end* - 1)
		while (posFiber < endFiber && posQuery < endQuery && mismatches <= localK + 1
				&& consecutiveMismatches <= consErrors) {
			if (fiber[posFiber] != query[posQuery]) {
				rMmOffsets.push_back(posFiber - getEndDim0(seed));
				mismatches++;
				if (posFiber == endFiber - 1 || posQuery == endQuery - 1) {
					isRightEndIncluded = true;
				}
				++consecutiveMismatches;
			}
			else {
				consecutiveMismatches = 0;
			}
			localK = std::max((int)floor((posFiber - getBeginDim0(seed) + 1) * errorRate), k);
			++posFiber;
			++posQuery;
		}
		// add corresponding end points if required
		if ((posFiber == endFiber || posQuery == endQuery) && !isRightEndIncluded) {
			rMmOffsets.push_back(std::min((int)(length(query) - getEndDim1(seed) - 1), (int)(length(fiber) - getEndDim0(seed) - 1)));
		}

#ifdef DEBUG
		// ------------------------------------------------------------
		// debug print
		cout << "lMmOffsets: ";
		for (int i = 0; i < lMmOffsets.size(); ++i) {
			cout << lMmOffsets[i] << " ";
		}
		cout << endl << "rMmOffsets: ";
		for (int i = 0; i < rMmOffsets.size(); ++i) {
			cout << rMmOffsets[i] << " ";
		}
		cout << endl;
#endif
		// ============================================================================
		// START REAL SHIT
		// init left and right mismatch numbers
		int lMmSize = lMmOffsets.size();
		int rMmSize = rMmOffsets.size();
		int previousR = -1;

		for (int l = lMmSize - 1; l >= 0; --l) {
			for (int r = rMmSize - 1; r >= 0; --r) {
				bDim0 = getBeginDim0(seed) - lMmOffsets[l];
				bDim1 = getBeginDim1(seed) - lMmOffsets[l];
				eDim0 = getEndDim0(seed) + rMmOffsets[r];
				eDim1 = getEndDim1(seed) + rMmOffsets[r];
        		int cnt = 0;

//        		cout << "bDim0, bDim1:" << bDim0 << ", " << bDim1 << endl << std::flush;
//        		cout << "eDim0, eDim1:" << eDim0 << ", " << eDim1 << endl << std::flush;
        		// skip beginning positions until match found
        		while (fiber[bDim0] != query[bDim1]) {
        			if (cnt) {
        				--l;
//        				cout << "skip 1 bDim* with consequences, lower l: " << l << endl << std::flush;
        			}
        			else {
//        				cout << "skip 1 bDim* for FREE, they differ" << endl << std::flush;
        			}
        			++bDim0;
        			++bDim1;
        			++cnt;
        		}

        		cnt = 0;
        		while (fiber[eDim0] != query[eDim1]) {
        			if (cnt) {
        				--r;
//        				cout << "skip 1 eDim* with consequences" << endl << std::flush;
        			}
        			else {
//        				cout << "skip 1 eDim* for FREE, they differ" << endl << std::flush;
        			}
        			--eDim0;
        			--eDim1;
        			++cnt;
        		}

        		// bDim* and eDim* include the respective elements as match
        		int diag 	= eDim0 - bDim0 + 1; //
        		localK		= std::max((int)floor(diag * errorRate), k);
        		mismatches 	= l + r;

#ifdef DEBUG
        		cout << "eDim0, eDim1 after skipping:" << eDim0 << ", " << eDim1 << endl << std::flush;
        		cout << "bDim0, bDim1 after skipping:" << bDim0 << ", " << bDim1 << endl << std::flush;
        		cout << "diagonal (+1 incl.): " << diag << endl;
        		cout << "#mismatches " << mismatches << ", localK: " << localK << endl;
        		cout << "previousR: " << previousR << endl;
#endif

        		// make sure #mismatches is valid
        		if (localK < mismatches) {
//        			cout << "localK < mismatches" << endl;
        			continue;
        		}

        		if (r <= previousR) {
//        			cout << "r <= previousR: " << r << " <= " << previousR << endl;
        			continue;
        		}

        		// eDim* still include the last match, to report we add +1
        		if (eDim0 + 1 - bDim0 >= minLength) {
        			if (bDim0 >= 65536 || bDim1 >= 65536 || eDim0 >= 65536 || eDim0 >= 65536) {
        				cout << "YOU SUCK! find a new hash function." << endl << std::flush;
        				exit(-1);
        			}
        			long long hash = bDim0;
        			hash = (((((hash << 16) + eDim0) << 16) + bDim1) << 16) + eDim1;

        			// if new match
        			if (!addedSeedHashes.count(hash)) {
        				// add +1 to endDim*, it's how seeds / matches are reported
        				extendedSeeds.push_back(TSeed(bDim0, bDim1, eDim0 + 1, eDim1 + 1));
        				addedSeedHashes.insert(hash);
        				previousR = r;
//        				cout << "l: " << l << ", r: " << r << endl;
//        				cout << "valid & new seed extension: " << extendedSeeds.back() << endl;
        			}
        			else {
//        				cout << "Hmm, this is wrong... same hash already existing." << endl;
        			}
        		}
        		else {
//        			cout << "length of match is: " << eDim0 + 1 - bDim0 << " < " << minLength << endl;
        		}
			}
		}
	}



	/**
	 * If seed is new, add it to seedMap and return True, False otherwise.
	 * param: addedSeedHashes is the set of hashes for a given tts-tfo pair
	 */
	template<
	typename TId,
	typename TSeed,
	typename TSeedHashesSet
	>
	bool addIfNewSeed(
			TId 	const 	&haystackFiberSeqNo,
			TId 	const 	&tfoSeqNo,
			TSeed  	const	&seed,
			TSeedHashesSet	&addedSeedHashes) {

		unsigned long long hash = getBeginDim0(seed);
		hash = (((((hash << 16) + getEndDim0(seed)) << 16) + getBeginDim1(seed)) << 16) + getEndDim1(seed);

		if (addedSeedHashes.count(hash)) {
#ifdef DEBUG
			std::cout << "== XXX >> Seed is already in map: " << seed << std::endl
					<< "\thash: " << hash << std::endl
					<< "\ttfoSeqNo: " << tfoSeqNo << std::endl
					<< "\thaystackFiberSeqNo: " << haystackFiberSeqNo << std::endl;
#endif
			return false;
		}
		addedSeedHashes.insert(hash);
#ifdef DEBUG
		std::cout << "===>>> New seed is added to map: " << seed << std::endl
				<< "\thash: " << hash << std::endl
				<< "\ttfoSeqNo: " << tfoSeqNo << std::endl
				<< "\thaystackFiberSeqNo: " << haystackFiberSeqNo << std::endl;
#endif
		return true;
	}

	/**
	 *
	 */
	template<
	typename TSeedMap,
	typename THitList,
//	typename THitSetPointerMap,
	typename THaystack,
	typename TMergedHaystack,
	typename TSegment,
	typename TMotifSet,
	typename TSuffix,
	typename TSAIter,
	typename TError,
	typename TSize,
	typename TConsError,
	typename THit,
	typename THitListKey
	>
	void verifyMatches(
			TSeedMap	   			&addedSeedHashMap,
			THitList				&hitList,
//			THitSetPointerMap 		&hitSetMap,
			THaystack 		const 	&haystack,
			TMergedHaystack const 	&mergedHaystack,
			TSegment 		const 	&segmentMap,
			TMotifSet 		const	&needleSet,
			TSuffix 		const	&suffixQGram,
			TSAIter 		const	&itStartBucket,
			TSAIter 		const	&itEndBucket,
			TError 			const   &errorRate,
			int 		   	const	&numLocations,
	        int 			const	endLocations[],
	        TSize 			const	&minLength,
			TConsError		const	&consErrors,
			THit,
			THitListKey)
	{
		typedef int								TScore;
		typedef Seed<Simple, DefaultSeedConfig>	TSeed;
		typedef std::vector< Seed<Simple, DefaultSeedConfig> > 			 TSeedList;
		typedef std::vector< Seed<Simple, DefaultSeedConfig> >::iterator TSeedListIterator;

	    int ePos;
	    int bPos;
	    int qPos;
	    int misM;
	    int haystackFiberSeqNo;
	    // max. seed positions
	    int maxSeedLength;
	    int maxSeedMergedEnd;	// global position in merged haystack
	    int maxSeedQGramEnd;	// local position in
	    int maxSeedFiberEnd; 	// local begin position in fiber (elem in haystack)
	    int k = ceil(errorRate * minLength); // #mismatches allowed
	    TSAIter itSB, itEB;
	    TSeedList extendedSeeds;

	    // iterate over all putative matches (end locations), find max seed and then extend
	    for (int match = 0; match < numLocations; match++) {
#ifdef DEBUG
//	    	cout << endl << "Iteration #" << match << " =====================================" << endl;
#endif
	        ePos = endLocations[match]; 	// end of current putative match in duplex (mergedHaystack)
	        bPos = ePos - minLength + 1; 	// beginning of --||--
	        if (bPos < 0) {
	        	continue;
	        }

	        // reset position variables
	        misM = 0;
	        qPos = minLength - 1;
	        maxSeedLength 		= 0;
	        maxSeedMergedEnd	= 0;
	        maxSeedQGramEnd 	= 0;
	        haystackFiberSeqNo 	= getHaystackFiberSeqNo(bPos, segmentMap);
	        assert(haystackFiberSeqNo >= 0);

	        // make sure the found match doesn't span 2 different fibers in the original haystack
	        if ((bPos - segmentMap[haystackFiberSeqNo]) + minLength > length(haystack[haystackFiberSeqNo])) {
	        	continue;
	        }

	        // reset SA iterators
	        itSB = itStartBucket;
	        itEB = itEndBucket;

	        /////////////////////////////
	        // validate TODO check for # of consecutive mismatches
	        int consecutiveMismatches = 0;
	        for (int pos = ePos; pos > bPos && misM <= k; --pos) {
	        	if (mergedHaystack[pos] != suffixQGram[qPos--]) {
	        		++misM;
	        		++consecutiveMismatches;
	        		if (consecutiveMismatches > consErrors) {
	        			continue;
	        		}
	        	}
	        	else {
	        		consecutiveMismatches = 0;
	        	}
	        }
	        // invalid alignment
	        if (misM > k) {
	        	continue;
	        }

#ifdef DEBUG
	        cout 	<< endl << "=====>>>>> Found valid match @ " << endl;
//					<< "\tmismatches: " << misM << endl
//	        		<< "haystackFiberSeqNo (in original haystack): " << haystackFiberSeqNo << endl
	        cout
	        		<< "Seed match is (Myers):" << endl << "\t";
	        for (int pos = bPos; pos <= ePos; ++pos) {
	        	std::cout << mergedHaystack[pos];
	        }
	        std::cout << "\ttts" << endl << "\t";
	        for (int pos = 0; pos < minLength; ++pos) {
	        	std::cout << suffixQGram[pos];
	        }
	        std::cout << "\ttfo" << endl;

#endif
	        // find the largest seed within the q-gram
	        qPos = minLength - 1;
	        for (int pos = ePos; pos > bPos && qPos >= 0; --pos, --qPos) {
	        	int matchLength = 0;

	        	while (qPos >= 0 && mergedHaystack[pos] == suffixQGram[qPos]) {
	        		matchLength++;
	        		--pos;
	        		--qPos;
	        	}

	        	if (matchLength > maxSeedLength) {
	        		maxSeedLength 		= matchLength;
	        		maxSeedMergedEnd 	= pos + matchLength;
	        		maxSeedQGramEnd		= qPos + matchLength;
	        		maxSeedFiberEnd		= maxSeedMergedEnd - segmentMap[haystackFiberSeqNo];
	        	}
	        }

#ifdef DEBUG
	        cout << endl << "Max. seed: " << endl << "\t";
	        for (int pos = maxSeedMergedEnd - maxSeedLength + 1; pos <= maxSeedMergedEnd; ++pos) {
	        	cout << mergedHaystack[pos];
	        	assert(mergedHaystack[pos] == haystack[haystackFiberSeqNo][pos - segmentMap[haystackFiberSeqNo]]);
	        }
	        cout 	<< endl << "MaxSeed Length: " << maxSeedLength << std::endl
	        		<< "MaxSeed Position (mergedHaystack): " <<  maxSeedMergedEnd - maxSeedLength + 1
					<< " - " << maxSeedMergedEnd << std::endl;
#endif

	        // found a valid match, extend for each tfo which had this q-gram suffix
	        for (; itSB != itEB; ++itSB) {
	            int ndlSeqNo      	= itSB->i1;
	            int qGramOffset  	= itSB->i2; // start position / offset of q-gram in tfo
	            int qGramSeedBegin  = maxSeedQGramEnd - maxSeedLength + 1 + qGramOffset;
	            int qGramSeedEnd 	= maxSeedQGramEnd + qGramOffset;

	            THitListKey seqNoKey(haystackFiberSeqNo, ndlSeqNo);
	            TSeed seed(maxSeedFiberEnd - maxSeedLength + 1, qGramSeedBegin, maxSeedFiberEnd, qGramSeedEnd);
#ifdef DEBUG
	            std::cout << endl << "Original seed: " << seed << endl
	            		<< "Actual right ends (including this pos.): " << getEndDim0(seed) << ", " << getEndDim1(seed) << endl
	            		<< "seedFiber: " << infix(haystack[haystackFiberSeqNo], getBeginDim0(seed), getEndDim0(seed) + 1) << "\n"
						<< "seedQuery: " << infix(needleSet[ndlSeqNo], getBeginDim1(seed), getEndDim1(seed) + 1) << "\n";
#endif
	            extendSimpleSeedRevised(addedSeedHashMap[seqNoKey], seed, extendedSeeds, haystack[haystackFiberSeqNo],
	            		needleSet[ndlSeqNo], errorRate, k, minLength, consErrors);

	            for (TSeedListIterator seed = extendedSeeds.begin(); seed != extendedSeeds.end(); ++seed) {
#ifdef DEBUG
	            	cout << "++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
					cout << "one seed after extension: " << *seed << endl;
					cout << "Match: " << endl;
	            	cout << "\tseedFiber: " << infix(haystack[haystackFiberSeqNo], getBeginDim0(*seed), getEndDim0(*seed)) << endl;
					cout << "\tseedQuery: " << infix(needleSet[ndlSeqNo], getBeginDim1(*seed), getEndDim1(*seed)) << endl;
	            	cout << "needle (tfo - isparallel): " << isParallel(needleSet[ndlSeqNo]) << "\t\t\t" << needleSet[ndlSeqNo] << endl;
	            	cout << "match length: " << getEndDim0(*seed) - getBeginDim0(*seed) << endl;
	            	cout << "tfo sequence id (index): \t\t" << ndlSeqNo << endl;
	            	cout << "+++++ ------- +++++++ ------ ++++++ will try to add new? seed" << endl;
#endif
	            	// if new seed, add to seedMap and to hitSet
	            	if (addIfNewSeed(haystackFiberSeqNo, ndlSeqNo, *seed, addedSeedHashMap[seqNoKey])) {
	            		THit *hit = new THit(
	            				haystackFiberSeqNo,
								ndlSeqNo,
								getBeginDim0(*seed),
								getBeginDim1(*seed),
								-1, //diag,					// the diagonal
								0,
								getEndDim0(*seed)-getBeginDim0(*seed)
	            		);

	            		THitListKey matchKey(getBeginDim0(*seed), getEndDim0(*seed));

	            		cout << "+++ adding new hit for matchKey<int,int> : " << matchKey.first << ", " << matchKey.second << endl;
	            		((hitList[seqNoKey])[matchKey]).push_back(hit);
	            	}
	            }
	        }
	    }
	}

	/**
	 *
	 */
	template<
	typename THaystack,		// haystack spec - double stranded sequences (tts)
	typename TMergedHaystack,
	typename TSegmentMap
	>
	void mergeHaystack(
			THaystack const 	&haystack,
			TMergedHaystack 	&mergedHaystack,
			TSegmentMap  		&segmentMap,
			unsigned char  	   *&bitTarget,
			unsigned int 		&totalLength
			) {
		assert(!length(mergedHaystack));

		// calculate total length of segments in haystack
		for (int i = 0; i < length(haystack); ++i) {
			totalLength += length(haystack[i]);
		}

#ifdef DEBUG
		assert(totalLength == length(mergedHaystack));
		std::cout 	<< "Merged haystack has length: " << totalLength << std::endl;
		std::cout 	<< "Fibers in merged haystack: " << std::endl;
#endif

		// resize merged haystack to totalLength
		resize(mergedHaystack, totalLength, Exact());
		bitTarget = new unsigned char[totalLength];
		totalLength = 0;

		for (int i = 0; i < length(haystack); ++i) {
			// add starting position in mergedHaystack for current fiber
			segmentMap.push_back(totalLength);

#ifdef DEBUG
			std::cout 	<< "seq no: " << i << " (len = " << length(haystack[i]) << "): ";
#endif
			// fill mergedHaystack and Myers vector
			for (int j = 0; j < length(haystack[i]); ++j) {
				mergedHaystack[totalLength + j] = haystack[i][j];
				// rewrite double stranded seq (tts) in a proper form for Myers (char -> index)
				bitTarget[totalLength + j] = charToBit(haystack[i][j]);//static_cast<unsigned int>(haystack[i][j]);
			}
			totalLength += length(haystack[i]);
#ifdef DEBUG
			std::cout 	<< haystack[i] << std::endl;
#endif
		}

#ifdef DEBUG
		std::cout << "Segment map / vector:\t";
		for (unsigned i = 0; i < segmentMap.size(); ++i) {
			cout << segmentMap[i] << " x ";
		}
		cout << endl;
#endif
	}

	template<typename TChar>
	unsigned int charToBit(TChar const &c) {
		switch(char(c)) {
		case 'A': return 0;
		case 'C': return 1;
		case 'G': return 2;
		case 'T': return 3;
		case 'Y': return 4;
		case 'N': return 5;
		default: assert(false && "Invalid character.");
		}
	}

	template<
		typename THitSetPointerMap,
		typename THitList,
		typename THit
	>
	void mergeOverlappingHits(
			THitList 			&hitList,
			THitSetPointerMap 	&hitSetMap,
			THit)
	{
		typedef std::pair<int,int>										THitListKey;
		typedef typename 	THitList::iterator 							THitListIterator;
		typedef typename 	std::vector<THit*>::iterator 				THitIterator;
		typedef typename 	std::map<THitListKey, std::vector<THit*> >::iterator TDim0Iterator;

		int haystackFiberSeqNo;
		TDim0Iterator dim0It;
		TDim0Iterator nextDim0It;

		cout << endl << endl << "Merging overlaps" << endl;

		for (THitListIterator it = hitList.begin(); it != hitList.end(); ++it) {
			haystackFiberSeqNo = (it->first).first;

			cout << "Starting new tts-tfo pair" << endl;

			for (dim0It = (it->second).begin(); dim0It != (it->second).end(); ++dim0It) {
				int currBegDim0 = dim0It->first.first;
				int currEndDim0 = dim0It->first.second;

				nextDim0It = dim0It;
				bool noMoreOverlap = false;
				for (++nextDim0It; nextDim0It != (it->second).end() && !noMoreOverlap; ++nextDim0It) {

					int nextBegDim0 = nextDim0It->first.first;
					int nextEndDim0 = nextDim0It->first.second;

					cout << endl << endl << "current matchKey<int,int> : " << currBegDim0 << ", " << currEndDim0 << endl << std::flush;
					cout << "next matchKey<int,int> : " << nextBegDim0 << ", " << nextEndDim0 << endl << std::flush;

					// check if dim0 indices overlap
					if (nextBegDim0 >= currBegDim0 && nextBegDim0 < currEndDim0) {
						// iterate over all hit matches
						for (THitIterator currHitIt = dim0It->second.begin(); currHitIt != dim0It->second.end();) {
							bool moveIt = true;
							for (THitIterator nextHitIt = nextDim0It->second.begin(); nextHitIt != nextDim0It->second.end();) {
								THit *currHit = *currHitIt;
								THit *nextHit = *nextHitIt;
								cout 	<< "c Hit: ndlPos = " << currHit->ndlPos << endl
										<< "n Hit: ndlPos = " << nextHit->ndlPos << endl
										<< "c Hit: hstkPos = " << currHit->hstkPos << endl
										<< "n Hit: hstkPos = " << nextHit->hstkPos << endl
										<< "c Hit: hstkPos - ndlPos = " << currHit->hstkPos - currHit->ndlPos << endl
										<< "n Hit: hstkPos - ndlPos = " << nextHit->hstkPos - nextHit->ndlPos << endl;


								// case I: currHit contains nextHit
								if (currHit->ndlPos <= nextHit->ndlPos
										&& currHit->ndlPos + currHit->hitLength >= nextHit->ndlPos + nextHit->hitLength
										&& currHit->hstkPos - currHit->ndlPos == nextHit->hstkPos - nextHit->ndlPos 	// alignment must also match
								) {
									cout << "current hit contains next hit, delete next" << std::flush;
									delete nextHit;
									nextHitIt = nextDim0It->second.erase(nextHitIt);
								}

								// case II: nextHit contains currHit
								else if (currHit->ndlPos >= nextHit->ndlPos
										&& currHit->ndlPos + currHit->hitLength <= nextHit->ndlPos + nextHit->hitLength
										&& currHit->hstkPos - currHit->ndlPos == nextHit->hstkPos - nextHit->ndlPos		// alignment must also match
								) {
									cout << "next hit contains current hit, just delete current and skip to next current" << std::flush;
									delete currHit;
									currHitIt = dim0It->second.erase(currHitIt);
									moveIt = false;
									break;
								}
								else {
									++nextHitIt;
								}
							}
							if (moveIt) {
								++currHitIt;
							}
						}
					}
					// no overlap --> add all hits and stop verifying current dim0It
					else {
						cout << "no overlap! stop after this" << endl;
						noMoreOverlap = true;
					}
				}

				cout << "will add all remaining current hits, nextDim0It loop exited (or didn't enter)" << endl;
				// add all remaining hits and stop verifying current dim0It
				for (THitIterator hit = dim0It->second.begin(); hit != dim0It->second.end(); ++hit) {
					add((*hitSetMap[haystackFiberSeqNo]), **hit);
				}
			} // end dim0It for loop
		} // end tts-tfo match for loop
	}

	/**
	 * One parameter is the q-gram index of the TFO set. This function iterates over each double
	 * stranded sequence in the haystack, and for each of these sequences iterates over each substring
	 * s of length minLength (== len(q-gram)). For each s we do Myers' approximate matching with maximal
	 * k (# allowed errors), and then verify / extend these matches (seeds).
	 * TODO assume we have only one duplex sequence
	 */
	template<
	typename THaystack,		 // haystack spec - double stranded sequences (tts)
	typename TQGramIndex,	 // q-gram index - (tfo)
	typename TQuerySet,		 // query set (needle) - single stranded sequences (tfo)
	typename TError,		 // error rate
	typename TSize,			 // minimum hit size
	typename TSpec,			 // specialization
	typename TId,			 // sequence id
	typename TWorker
	>
	void plantMyers(
			Gardener<TId, TSpec>	&gardener,
			THaystack const			&haystack, 			// TTS set
			TQGramIndex				&index,				// q-gram index of tfoSet (for suffixes)
			TQuerySet				&needles,			// TFO set
			TError const			&errorRate,
			TSize const				&minLength,
			int const				&consErrors,
			TWorker
	){
	    double t = sysTime();
		typedef typename Value<Gardener<TId, TSpec> >::Type		THitMap;
		typedef typename Value<THitMap>::Type					THitMapEntry;
		typedef typename Value<THitMapEntry,2>::Type			THitSetPointer;
		typedef typename Value<THitSetPointer>::Type			THitSet;
		typedef typename Value<THitSet>::Type					THit;

	    typedef typename Value<THaystack>::Type               	THaystackValue;
	    typedef typename Fibre<TQGramIndex, QGramSA>::Type      TSA;
	    typedef typename Fibre<TQGramIndex, QGramSADir>::Type   TSADir;
	    typedef typename Value<TSADir>::Type                    TSADirValue;
	    typedef typename Iterator<TSA, Standard>::Type          TSAIter;
	    typedef typename Iterator<TSADir, Standard>::Type       TSADirIter;
		typedef Seed<Simple, DefaultSeedConfig>					TSeed;


		typedef 		 std::vector<int>								TSegments;
		typedef  		 String<char> 									TMergedHaystack;
		typedef 		 std::map<int, THitSetPointer> 					THitSetPointerMap;
		typedef			 std::pair<int, int> 							THitListKey;
		typedef 		 std::map<THitListKey, std::set<unsigned long long> > 				TSeedHashesMap;
		typedef typename Suffix<typename GetSequenceByNo<TQGramIndex const>::Type >::Type	TSuffix;
		typedef			 std::map<THitListKey, std::map<THitListKey, std::vector<THit*> > >	THitList;

	    int score;
	    int numLocations;
	    int *endLocations;
	    int *startLocations;
	    int alignmentLength;
	    unsigned char *alignment;
		unsigned char *bitTarget;
		unsigned char bitQGram[minLength];

	    int k 			 = ceil(errorRate * minLength);
	    int cnt 		 = 0;
		unsigned targetLength 	= 0;
		unsigned alphabetSize	= 6;	// 4 + Y and N

		THitList		hitList;		// for each tts-tfo pair: keeps a
										//
		TSeedHashesMap 	seedHashMap;
	    TSADirIter  	itBucketDir;
	    TSADirIter  	itBucketDirEnd;
	    TSAIter     	saBegin;
	    TSADirValue 	bucketBegin;	// ptr to each unique qgram in SA of tfoSet (each bucket
	    								// may have several entries
	    TMergedHaystack mergedHaystack;	// current item in haystack (tts)
	    TSuffix 		suffixQGram;	// q-gram of appropriate (minimum) length within suffix
	    TSegments		segmentMap; 	// for duplex each segment, the map stores its starting
	    								// position in the merged haystack; sorted --> binary search

		THitSetPointerMap hitSetPointerMap; // for each query in queries (tfoSet), stores a
											// hitSetPointer containing all corresponding hits

		// add a new hitset object for each fiber in haystack
		for (int i = 0; i < length(haystack); ++i) {
			hitSetPointerMap[i] = new THitSet;
		}

	    // merge duplex segments
	    mergeHaystack(haystack, mergedHaystack, segmentMap, bitTarget, targetLength);

	    // reset SA iterators to beginning of index
	    itBucketDir         = begin(indexDir(index), Standard());
	    itBucketDirEnd      = end(indexDir(index), Standard());
	    saBegin             = begin(indexSA(index), Standard());
	    bucketBegin 		= *itBucketDir;

	    // iterate over each substring of length 'minLength' == same as iterating over each
	    // bucket (unique q-gram) in the q-gram index
	    for (++itBucketDir; itBucketDir != itBucketDirEnd; ++itBucketDir)
	    {
	    	const TSADirValue bucketEnd = *itBucketDir; // end of this bucket == beginning of next one
	    	if (bucketBegin != bucketEnd)
	    	{
	    		const TSAIter itBucketItem = saBegin + bucketBegin;
	    		const TSAIter itEndBucket  = saBegin + bucketEnd;
	    		cnt ++;
	    		suffixQGram = suffix(indexText(index), *itBucketItem);

	    		// transform query into index for Myers and update array
				for (int i = 0; i < minLength; ++i) {
					bitQGram[i] = charToBit(suffixQGram[i]);//static_cast<unsigned int>(suffixQGram[i]);
				}

	    		//calculate Myers distance - this yields putative matches in endLocations
	    		edlibCalcEditDistance(bitQGram, minLength, bitTarget, targetLength, alphabetSize,
	    				k, EDLIB_MODE_HW, true, true, &score,
						&endLocations, &startLocations, &numLocations,
						&alignment, &alignmentLength);

#ifdef DEBUG
//	    		cout 	<< endl << endl
//	    				<< "----------------------------------------" << endl
//	    				<< "----------- NEW QGRAM-SUFFIX -----------" << endl
//	    				<< "----------------------------------------" << endl;
//	    		std::cout << "current suffix (length == " << minLength << "): ";
//	    		for (int i = 0; i < minLength; ++i) {
//	    			printf("%c",(char)suffixQGram[i]);//(unsigned int)(char)(suf[i]));
//	    		}
//	    		cout << endl << "current target (length == " << targetLength << "):";
//	    		for (int i = 0; i < targetLength; ++i) {
//	    			cout << (int)((char)bitTarget[i]);
//	    			assert((int)((char)bitTarget[i]) <=5 && (int)((char)bitTarget[i]) >= 0);
//	    		}
	    		cout << endl << "Number of hits (numLocations): " << numLocations << endl;
#endif

	    		verifyMatches(seedHashMap, hitList, haystack, mergedHaystack, segmentMap, needles,
	    				suffixQGram, itBucketItem, itEndBucket, errorRate, numLocations,
						endLocations, minLength, consErrors, THit(), THitListKey());
	    	}
	    	bucketBegin = bucketEnd;
	    }
	    // keep only maximum segments, remove those included in larger hits
	    mergeOverlappingHits(hitList, hitSetPointerMap, THit());
	    // add hits to gardener
	    for (int i = 0; i < length(haystack); ++i) {
	    	insert(gardener.hits, i, hitSetPointerMap[i]);
	    }
	    delete bitTarget;
	}

	/** 
	 * start gardening by planting
	 * prepeare pattern beforehand for reuse
	 * all putative TTSs per duplex will be processed in parallel mode
	 * (for short duplex sequences with presumably few TTSs)
	 */
	template< 
	typename TIndex,		// index 
	typename TPatternSpec,	// pattern spec
	typename TShape,		// shape
	typename TQuerySet,		// query set (needle)
	typename TError,		// error rate
	typename TSize,			// minimum hit size
	typename TDrop,			// xdrop
	typename TSpec,			// specialization
	typename TId,			// sequence id
	typename TWorker
	>
	void plant(Gardener<TId, TSpec>	&gardener,
			   Pattern<TIndex, QGramsLookup< TShape, TPatternSpec> > const	&pattern, // index of tfoSet patterns (originally), now tts
			   TQuerySet			&queries,   // tts
			   TError const			&errorRate,
			   TSize const			&minLength,
			   TDrop const			&xDrop,
			   TWorker
			   ){
		// reset local timers
		timeCollectSeeds = 0.0;
		timeGardenerFind = 0.0;
		timePutSeedsInMap = 0.0;

		typedef typename Iterator<TQuerySet>::Type									TQueryIter;
		typedef typename Value<Gardener<TId, TSpec> >::Type							THitMap;
		typedef typename Value<THitMap>::Type										THitMapEntry;
		typedef typename Value<THitMapEntry,2>::Type								THitSetPointer;
		typedef typename Value<THitSetPointer>::Type								THitSet;
		typedef typename Value<TQuerySet>::Type										TSequence;
		typedef Finder<TSequence, QGramsLookup< TShape, Standard_QGramsLookup > >	TFinder;
		typedef typename Position<TFinder>::Type									TPos;

		typedef Pattern<TIndex, QGramsLookup< TShape, TPatternSpec> >			 	TPattern;
		typedef typename Iterator< TPattern >::Type									TPatternIter;
		
		// q-gram lemma
		// w+1-(k+1)q | w=minimum length, k=errors, q=weight(q-grams)
		TPos minSeedsThreshold = static_cast<TPos>(minLength+1-(ceil(errorRate*minLength)+1)*weight(pattern.shape));
#ifdef TRIPLEX_DEBUG
        ::std::cout << "\n\nPlant called ------------------------------ ==============================================\n";
        std::cout << "min length: " << minLength << std::endl
        		<< "errors (k): " << ceil(errorRate*minLength)+1 << std::endl
        		<< "weight (q-grams): " << weight(pattern.shape) << std::endl;
#endif						
		
		// serial processing
		TId querylen = (TId)length(queries);
		for (TId queryid=0; queryid<querylen; ++queryid){
#ifdef TRIPLEX_DEBUG
			::std::cout << "TTS (originally, now TFO) " << queryid << " : " << queries[queryid] << ::std::endl;
			::std::cout << "TTS (originally, now TFO) as ttsString: " << queryid << " : " << ttsString(queries[queryid]) << ::std::endl;
			_printHits(gardener, pattern, queries, queryid);
#endif

			THitSetPointer hitsPointer = new THitSet;
			TFinder finder(queries[queryid]); 
			_find(*hitsPointer, finder, pattern, errorRate, (TPos) minLength, minSeedsThreshold, xDrop, queryid );	
			insert(gardener.hits, queryid, hitsPointer);
			gardener.timeQgramFind += finder.timeFind;
		}		
		gardener.timeCollectSeeds += timeCollectSeeds;
		gardener.timeGardenerFind += timeGardenerFind;
		gardener.timePutSeedsInMap+= timePutSeedsInMap;
	}
		
	/** 
	 * start gardening by planting
	 * prepeare pattern beforehand for reuse
	 * all putative TTSs per duplex will be processed in parallel mode
	 * (for short duplex sequences with presumably few TTSs)
	 */
	template< 
	typename TIndex,		// index 
	typename TPatternSpec,	// pattern spec
	typename TShape,		// shape
	typename TQuerySet,		// query set (needle)
	typename TError,		// error rate
	typename TSize,			// minimum hit size
	typename TDrop,			// xdrop
	typename TSpec,			// specialization
	typename TId,			// sequence id
	typename TRepeat,		// repeat minimum length
	typename TWorker
	>
	void plant(Gardener<TId, TSpec>	&gardener,
			   Pattern<TIndex, QGramsLookup< TShape, TPatternSpec> > const	&pattern,
			   TQuerySet			&queries,
			   TError const			&errorRate,
			   TSize const			&minLength,
			   TDrop const			&xDrop,
			   TRepeat const		&minRepeatLength,
			   TRepeat const		&maxRepeatPeriod,
			   TWorker
			   ){
		typedef typename Iterator<TQuerySet>::Type									TQueryIter;
		typedef typename Value<Gardener<TId, TSpec> >::Type							THitMap;
		typedef typename Value<THitMap>::Type										THitMapEntry;
		typedef typename Value<THitMapEntry,2>::Type								THitSetPointer;
		typedef typename Value<THitSetPointer>::Type								THitSet;
		typedef typename Value<TQuerySet>::Type										TSequence;
		typedef Finder<TSequence, QGramsLookup< TShape, Standard_QGramsLookup > >	TFinder;
		typedef typename Position<TFinder>::Type									TPos;
		
		// q-gram lemma
		// w+1−(k+1)q | w=minimum length, k=errors, q=weight(q-grams)
		TPos minSeedsThreshold = static_cast<TPos>(minLength+1-(ceil(errorRate*minLength)+1)*weight(pattern.shape));
#ifdef TRIPLEX_DEBUG
		::std::cout << "minLength:" << minLength << " errorRate:" << errorRate << " qgram:" << weight(pattern.shape) << ::std::endl;
		::std::cout << (ceil(errorRate*minLength)+1) << " " << ((ceil(errorRate*minLength)+1)*weight(pattern.shape)) << " " << minLength+1-(ceil(errorRate*minLength)+1)*weight(pattern.shape) << ::std::endl;
#endif								
		// serial processing
		TId querylen = (TId)length(queries);
		for (TId queryid=0; queryid<querylen; ++queryid){
			THitSetPointer hitsPointer = new THitSet;
			TFinder finder(queries[queryid], minRepeatLength, maxRepeatPeriod); 
			_find(*hitsPointer, finder, pattern, errorRate, (TPos) minLength, minSeedsThreshold, xDrop, queryid );	
			insert(gardener.hits, queryid, hitsPointer);
#ifdef TRIPLEX_DEBUG
			::std::cout << "TTS " << queryid << " : " << queries[queryid] << ::std::endl;
			_printHits(gardener, pattern, queries, queryid);
#endif
		}		
	}
		
#if SEQAN_ENABLE_PARALLELISM	
	/** 
	 * start gardening by planting
	 * prepeare pattern beforehand for reuse
	 * all putative TTSs per duplex will be processed in parallel mode
	 * (for long duplex sequences with presumably many TTSs)
	 */
	template< 
	typename TIndex,		// index 
	typename TPatternSpec,	// pattern spec
	typename TShape,		// shape
	typename TQuerySet,		// query set (needle)
	typename TError,		// error rate
	typename TSize,			// minimum hit size
	typename TDrop,			// xdrop
	typename TSpec,			// specialization
	typename TId			// sequence id
	>
	void plant(Gardener<TId, TSpec>	&gardener,
			   Pattern<TIndex, QGramsLookup< TShape, TPatternSpec> > const &pattern,
			   TQuerySet			&queries,
			   TError const			&errorRate,
			   TSize const			&minLength,
			   TDrop const			&xDrop,
			   MULTIPLE_WORKER
			   ){ ENTER
		typedef typename Iterator<TQuerySet>::Type									TQueryIter;
		typedef typename Value<Gardener<TId, TSpec> >::Type							THitMap;
		typedef typename Value<THitMap>::Type										THitMapEntry;
		typedef typename Value<THitMapEntry,2>::Type								THitSetPointer;
		typedef typename Value<THitSetPointer>::Type								THitSet;
		typedef typename Value<TQuerySet>::Type										TSequence;
		typedef Finder<TSequence, QGramsLookup< TShape, Standard_QGramsLookup > >	TFinder;
		typedef typename Position<TFinder>::Type									TPos;
		
		TId querylen = (TId)length(queries);

		// q-gram lemma
		// w+1−(k+1)q | w=minimum length, k=errors, q=weight(q-grams)
		TPos minSeedsThreshold = static_cast<TPos>(minLength+1-(ceil(errorRate*minLength)+1)*weight(pattern.shape));
#ifdef TRIPLEX_DEBUG
        std::cout << "bla\n";
		::std::cout << "minLength:" << minLength << " errorRate:" << errorRate << " qgram:" << weight(pattern.shape) << ::std::endl;
		::std::cout << (ceil(errorRate*minLength)+1) << " " << ((ceil(errorRate*minLength)+1)*weight(pattern.shape)) << " " << minLength+1-(ceil(errorRate*minLength)+1)*weight(pattern.shape) << ::std::endl;
#endif		
		
		SEQAN_PRAGMA_IF_PARALLEL(omp parallel) 
		{
			SEQAN_PRAGMA_IF_PARALLEL(omp for schedule(dynamic) )
			for (TId queryid=0; queryid<querylen; ++queryid){
				THitSetPointer hitsPointer = new THitSet;
				TFinder finder(queries[queryid]);
				_find(*hitsPointer, finder, pattern, errorRate, (TPos) minLength, minSeedsThreshold, xDrop, queryid );	
				
				SEQAN_PRAGMA_IF_PARALLEL(omp critical(addhitmap)  )
				insert(gardener.hits, queryid, hitsPointer);
			}
		}
	}
	
	/** 
	 * start gardening by planting
	 * prepeare pattern beforehand for reuse
	 * all putative TTSs per duplex will be processed in parallel mode
	 * (for long duplex sequences with presumably many TTSs)
	 */
	template< 
	typename TIndex,		// index 
	typename TPatternSpec,	// pattern spec
	typename TShape,		// shape
	typename TQuerySet,		// query set (needle)
	typename TError,		// error rate
	typename TSize,			// minimum hit size
	typename TDrop,			// xdrop
	typename TSpec,			// specialization
	typename TId,			// sequence id
	typename TRepeat		// repeat minimum length
	>
	void plant(Gardener<TId, TSpec>	&gardener,
			   Pattern<TIndex, QGramsLookup< TShape, TPatternSpec> > const &pattern,
			   TQuerySet			&queries,
			   TError const			&errorRate,
			   TSize const			&minLength,
			   TDrop const			&xDrop,
			   TRepeat const		&minRepeatLength,
			   TRepeat const		&maxRepeatPeriod,
			   MULTIPLE_WORKER
			   ){
		typedef typename Iterator<TQuerySet>::Type									TQueryIter;
		typedef typename Value<Gardener<TId, TSpec> >::Type							THitMap;
		typedef typename Value<THitMap>::Type										THitMapEntry;
		typedef typename Value<THitMapEntry,2>::Type								THitSetPointer;
		typedef typename Value<THitSetPointer>::Type								THitSet;
		typedef typename Value<TQuerySet>::Type										TSequence;
		typedef Finder<TSequence, QGramsLookup< TShape, Standard_QGramsLookup > >	TFinder;
		typedef typename Position<TFinder>::Type									TPos;
		
		TId querylen = (TId)length(queries);
		
		// q-gram lemma
		// w+1−(k+1)q | w=minimum length, k=errors, q=weight(q-grams)
		TPos minSeedsThreshold = static_cast<TPos>(minLength+1-(ceil(errorRate*minLength)+1)*weight(pattern.shape));
#ifdef TRIPLEX_DEBUG
        std::cout << "qwe\n";
		::std::cout << "minLength:" << minLength << " errorRate:" << errorRate << " qgram:" << weight(pattern.shape) << ::std::endl;
		::std::cout << (ceil(errorRate*minLength)+1) << " " << ((ceil(errorRate*minLength)+1)*weight(pattern.shape)) << " " << minLength+1-(ceil(errorRate*minLength)+1)*weight(pattern.shape) << ::std::endl;
#endif	
		
		// create hits pointer for all query ids
		::std::vector<THitSetPointer> tmpPointerList;
		
		for (TId queryid=0; queryid<querylen; ++queryid){
			THitSetPointer hitsPointer = new THitSet;
			appendValue(tmpPointerList, hitsPointer);
		}
		
		SEQAN_PRAGMA_IF_PARALLEL(omp parallel)
		{
			SEQAN_PRAGMA_IF_PARALLEL(omp for schedule(dynamic))
			for (TId queryid=0; queryid<querylen; ++queryid){
				TFinder finder(queries[queryid], minRepeatLength, maxRepeatPeriod); 
				_find(*tmpPointerList[queryid], finder, pattern, errorRate, (TPos) minLength, minSeedsThreshold, xDrop, queryid );	
			}
		}
		
		// copy all hits to the gardener
		for (TId queryid=0; queryid<querylen; ++queryid){
			insert(gardener.hits, queryid, tmpPointerList[queryid]);
		}
			
	}
		
#endif  // SEQAN_ENABLE_PARALLELISM
	
} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef FBUSKE_APPS_TRIPLEXATOR_HEADER_GARDENER_H


//	/**
//	 * Given a maximum seed match between fiber and query, this function returns all maximally
//	 * extended seeds which don't overlap entirely, w.r.t. to k and maximum lengths of the segments.
//	 * The ends of the seed to be extended are included in the match.
//	 */
//	template<
//	typename TSeed,
//	typename TSeedList,
//	typename THaystackFiber,
//	typename TQuery>
//	void extendSimpleSeed(
//			TSeed 			const &seed,
//			TSeedList			  &extendedSeeds,
//			THaystackFiber 	const &fiber,
//			TQuery 			const &query,
//			int 			const &k,
//			int				const &minLength) {
//		extendedSeeds.clear();
//
//		int posFiber;
//		int posQuery;
//		int mismatches;
//		int bDim0;
//		int bDim1;
//		int eDim0;
//		int eDim1;
//
//		bool isLeftZeroIncluded = false;
//		bool isRightEndIncluded = false;
//
//		// keep list of where mismatches occur, left and right from seed
//		std::vector<int> lMmOffsets;
//		std::vector<int> rMmOffsets;
//
//        // go left
//        mismatches = 0;
//        posFiber = getBeginDim0(seed);
//        posQuery = getBeginDim1(seed);
//        while (posFiber >= 0 && posQuery >= 0 && mismatches <= k) {
//            if (fiber[posFiber] != query[posQuery]) {
//                lMmOffsets.push_back(getBeginDim0(seed) - posFiber);
//                mismatches++;
//                if (posFiber == 0 || posQuery == 0) {
//                	isLeftZeroIncluded = true;
//                }
//            }
//            posFiber--;
//            posQuery--;
//        }
//
//        // go right
//        mismatches = 0;
//        posFiber = getEndDim0(seed);
//        posQuery = getEndDim1(seed);
//        int endFiber = length(fiber);
//        int endQuery = length(query);
//        // only go until the second last element, the last one is added later if needed (end* - 1)
//        while (posFiber < endFiber && posQuery < endQuery && mismatches <= k) {
//            if (fiber[posFiber] != query[posQuery]) {
//                rMmOffsets.push_back(posFiber - getEndDim0(seed));
//                mismatches++;
//                if (posFiber == endFiber - 1 || posQuery == endQuery - 1) {
//                	isRightEndIncluded = true;
//                }
//            }
//            posFiber++;
//            posQuery++;
//        }
//
//        // init left and right mismatch numbers
//        int lMmSize = lMmOffsets.size();
//        int rMmSize = rMmOffsets.size();
//
//        // if #mismatches less than k, whole segments match
//        if (lMmSize + rMmSize <= k) {
//#ifdef DEBUG
//            cout << "Whole segments match in seed extension, limits are included." << endl;
//#endif
//        	bDim0 = std::max(0, (int)(getBeginDim0(seed) - getBeginDim1(seed)));
//        	bDim1 = std::max(0, (int)(getBeginDim1(seed) - getBeginDim0(seed)));
//        	eDim0 = std::min((int)(length(fiber)), (int)(getEndDim0(seed) + length(query) - getEndDim1(seed)));
//        	eDim1 = std::min((int)(length(query)), (int)(getEndDim1(seed) + length(fiber) - getEndDim0(seed)));
//
//        	while (fiber[bDim0] != query[bDim1]) {
//        		++bDim0;
//        		++bDim1;
//        	}
//
//        	while (fiber[eDim0 - 1] != query[eDim1 - 1]) {
//        		--eDim0;
//        		--eDim1;
//        	}
//
//            if (eDim0 - bDim0 >= minLength) {
//            	extendedSeeds.push_back(TSeed(bDim0, bDim1, eDim0, eDim1));
//            }
//#ifdef DEBUG
//            else {
//            	cout << "length of match is: " << eDim0 - bDim0 << " < " << minLength << endl;
//            }
//#endif
//            return;
//        }
//
//        // add corresponding end points if required
//        // if leftmost offset (beg. of seed) was not included
//        if (lMmSize <= k && !isLeftZeroIncluded) {
//        	lMmOffsets.push_back(std::min((int)(getBeginDim0(seed)),(int)(getBeginDim1(seed))));
//        	++lMmSize;
//        }
//        // if rightmost offset (beg. of seed) was not included
//        if (rMmSize <= k && !isRightEndIncluded) {
//        	rMmOffsets.push_back(std::min((int)(length(query) - getEndDim1(seed) - 1), (int)(length(fiber) - getEndDim0(seed) - 1)));
//        	++rMmSize;
//        }
//
//#ifdef DEBUG
//        cout << "leftMismatches: " << lMmSize << endl;
//        cout << "rightMismatches: " << rMmSize << endl;
//        cout << "isLeftZeroIncluded: " << isLeftZeroIncluded << endl;
//        cout << "isRightEndIncluded: " << isRightEndIncluded << endl;
//#endif
//
//        std::set<int> addedSeedHashes;
//        std::vector<int> offsetPairs(lMmSize + 1, -1);
//
//        for (int l = lMmSize - 1; l >= 0; --l) {
//        	for (int r = rMmSize - 1; r >= 0 && l + r >= k; --r) {
//        		// lower r until match is valid, i.e., <= k
//        		if (l + r > k) {
//        			continue;
//        		}
//#ifdef DEBUG
//        		cout << "iter start l: " << l << endl;
//        		cout << "iter start r: " << r << endl;
//        		cout << "l + r: " << l + r << endl;
//#endif
//        		bDim0 = getBeginDim0(seed) - lMmOffsets[l];
//        		bDim1 = getBeginDim1(seed) - lMmOffsets[l];
//        		eDim0 = getEndDim0(seed) + rMmOffsets[r];
//        		eDim1 = getEndDim1(seed) + rMmOffsets[r];
//        		int cnt = 0;
//
//        		cout << "bDim0, bDim1:" << bDim0 << ", " << bDim1 << endl << std::flush;
//        		cout << "eDim0, eDim1:" << eDim0 << ", " << eDim1 << endl << std::flush;
//        		// skip beginning positions until match found
//        		while (fiber[bDim0] != query[bDim1]) {
//        			if (cnt) {
//        				--l;
//        				cout << "skip 1 bDim* with consequences, lower l: " << l << endl << std::flush;
//        			}
//        			else {
//        				cout << "skip 1 bDim* for FREE, they differ" << endl << std::flush;
//        			}
//        			++bDim0;
//        			++bDim1;
//        			++cnt;
//        		}
//
//        		cnt = 0;
//        		while (fiber[eDim0] != query[eDim1]) {
//        			if (cnt) {
//        				--r;
//        				cout << "skip 1 eDim* for FREE, they differ" << endl << std::flush;
//        			}
//        			else {
//        				cout << "skip 1 eDim* with consequences" << endl << std::flush;
//        			}
//        			--eDim0;
//        			--eDim1;
//        			++cnt;
//        		}
//        		cout << "eDim0, eDim1 after skipping:" << eDim0 << ", " << eDim1 << endl << std::flush;
//        		cout << "bDim0, bDim1 after skipping:" << bDim0 << ", " << bDim1 << endl << std::flush;
//
//        		// eDim* still include the last match, to report we add +1
//        		if (eDim0 + 1 - bDim0 >= minLength && offsetPairs[l + 1] < r && offsetPairs[l] < r) {
//        			int hash = ((bDim0 + 1) * (eDim1 + 1) * 10) - (bDim1 + 1) * (eDim0 + 1);
//        			// if new match
//        			if (!addedSeedHashes.count(hash)) {
//        				// add +1 to endDim*, it's how seeds / matches are reported
//        				extendedSeeds.push_back(TSeed(bDim0, bDim1, eDim0 + 1, eDim1 + 1));
//        				addedSeedHashes.insert(hash);
//        				offsetPairs[l] = r;
//#ifdef DEBUG
//        				cout << "l: " << l << ", r: " << r << endl;
//        				cout << "valid & new seed extension: " << extendedSeeds.back() << endl;
//#endif
//        			}
//        		}
//#ifdef DEBUG
//        		else {
//        			cout << "length of match is: " << eDim0 + 1 - bDim0 << " < " << minLength << endl;
//        			cout << "offsetPairs[l]: " << offsetPairs[l + 1] << ", r: " << r << endl;
//        		}
//#endif
//        	}
//        }
//
//        // find maximal seeds
////        for (; l >= 0 && r < rMmSize; --l) {
////            assert(abs(l - r) <= k);
////
////            bDim0 = getBeginDim0(seed) - lMmOffsets[l];
////            bDim1 = getBeginDim1(seed) - lMmOffsets[l];
////
////            // first skip of mismatch is "FREE"
////            if (fiber[bDim0] != query[bDim1]) {
////            	++bDim0;
////            	++bDim1;
////            	cout << "skip 1 for FREE" << endl << std::flush;
////            }
////
////            // skip beginning positions until match found
////            while (fiber[bDim0] != query[bDim1]) {
////            	cout << "skip 1 with consequences" << endl << std::flush;
////            	--l;
////            	++bDim0;
////            	++bDim1;
////
////            	// if we adjust the beginning positions, we can extend the right ones to get maximal seeds
////            	if (r < rMmSize - 1) {
////            		++r;
////            		cout << "skip r also " << endl << std::flush;
////            	}
////            }
////
////            cout << "bar r:" << r << endl << std::flush;
////            eDim0 = getEndDim0(seed) + rMmOffsets[r];
////            eDim1 = getEndDim1(seed) + rMmOffsets[r];
////            cout << "qux eDim0 eDim1:" << eDim0 << ", " << eDim1 << endl << std::flush;
////
////            while (fiber[eDim0] != query[eDim1]) {
////            	--eDim0;
////            	--eDim1;
////            }
////
////            ++eDim0;
////            ++eDim1;
////
////#ifdef DEBUG
////            cout << "l: " << l << ", r: " << r << endl;
////#endif
////            if (eDim0 - bDim0 >= minLength) {
////            	extendedSeeds.push_back(TSeed(bDim0, bDim1, eDim0, eDim1));
////            }
////#ifdef DEBUG
////            else {
////            	cout << "length of match is: " << eDim0 - bDim0 << " < " << minLength << endl;
////            }
////#endif
////            ++r;
////        }
//	}
