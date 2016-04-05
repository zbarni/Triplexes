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
#define QGRAM 3
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
	bool isOverlapping(TSeed const & a, TSeed const & b){ ENTER
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
								 ){ ENTER
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
#endif
						// extend seed to both sides as far as possible
						extendSeed(seed, tmp, getSequenceByNo(seqno,needle(pattern)), EXTEND_BOTH, scoreMatrix, scoreDropOff, UnGappedXDrop());

#ifdef TRIPLEX_DEBUG
						::std::cout << "seed_3:" << getBeginDim0(seed) << ":" << getEndDim0(seed) << " " << getBeginDim1(seed) << ":" << getEndDim1(seed) << ::std::endl;
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
						// create a new hit and append it to the gardeners hit list
						THit hit(queryid,
								 seqno,					// needle seq. number            
								 getBeginDim0(front(newset)),      // begin in haystack      
								 getBeginDim1(front(newset)),		// needle position
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
	        ::std::cout << "\n_collectSeeds(): new seed found: ******************************\n";
			::std::cout << "Q (inv = tfo):" << infix(finder) << ::std::endl; // The Infix of the match in the haystack.
			::std::cout << "T (inv = tts):" << infix(pattern, *finder.curHit) << ::std::endl;
			::std::cout << "H:" << (*finder.curHit).hstkPos << "-N" << (*finder.curHit).ndlSeqNo << ":P" << (*finder.curHit).ndlPos << ":D" << (*finder.curHit).diag << ::std::endl;
            //TODO@barni remove
			::std::cout << "needle itself@barni: " << getSequenceByNo((*finder.curHit).ndlSeqNo, needle(pattern)) << ::std::endl;
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
			std::cout << "match: " << match << "; mismatch: " << mismatch << "; scoredropoff: " << scoreDropOff << std::endl;
			std::cout << "tscore: " << (TScore)length(host(finder)) << std::endl;
//			exit(0);
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
	 *
	 */
	template<
	typename THitSet,
	typename TTts,
	typename TMotifSet,
	typename TSuffix,
	typename TSAIter,
	typename TError,
	typename TSize,
	typename TDrop
	>
	void verifyMatches(
			THitSet 	   &hitSet,
			TTts 			ttsFiber,
			TMotifSet 		tfoSet,
			TSuffix 		suffixQGram,
			const TSAIter 	itStartBucket,
			const TSAIter 	itEndBucket,
			TError const   &errorRate,
			int 		   &numLocations,
	        int 			endLocations[],
	        TSize const	   &minLength,
	        TDrop const    &xDrop) {
		typedef int								TScore;
		typedef Seed<Simple, DefaultSeedConfig>	TSeed;

	    int ePos;
	    int bPos;
	    int qPos;
	    int misM;
	    int maxSeedLength, maxSeedEndPos;
	    int k = ceil(errorRate * minLength);
	    TSAIter itSB, itEB;

	    // define a scoring scheme
	    TScore match 	= 1;
	    TScore mismatch = (TScore)_max((TScore) (-1.0/(errorRate+0.00000001)) + 1, -(TScore)length(ttsFiber));
	    Score<TScore> scoreMatrix(match, mismatch, std::numeric_limits<int>::max());
	    TScore scoreDropOff = (TScore) _max((TScore) xDrop * (-mismatch), minValue<TScore>()+1);
//	    std::cout << "MYERS match: " << match << "; mismatch: " << mismatch << "; scoredropoff: " << scoreDropOff << std::endl;

	    std::vector<std::pair<int, int> > seeds;

	    // iterate over all putative matches (end locations), find max seed and then extend
	    for (int loc = 0; loc < numLocations; loc++) {
	        ePos = endLocations[loc]; // end of current putative match
	        bPos = ePos - minLength;  // beginning of --||--
	        misM = 0;                 // nr of mismatches
	        qPos = minLength - 1;
	        maxSeedLength = 0;
	        maxSeedEndPos = 0;


	        // reset SA iterators
	        itSB = itStartBucket;
	        itEB = itEndBucket;

	        /////////////////////////////
	        // validate TODO check for # of consecutive mismatches
	        for (int pos = ePos; pos > bPos && misM <= k; --pos) {
	        	if (ttsFiber[pos] != suffixQGram[qPos--]) {
	        		++misM;
	        	}
	        }
	        // invalid alignment
	        if (misM > k) {
	        	continue;
	        }

	        // TODO barni make this more efficient
#ifndef DEBUG
	        std::cout << "TTSFIBER substring:\n";
	        for (int pos = bPos + 1; pos <= ePos; ++pos) {
	        	std::cout << ttsFiber[pos];
	        }
	        std::cout << endl;
	        std::cout << "TFOSUFFIX substring:\n";
	        for (int pos = 0; pos < minLength; ++pos) {
	        	std::cout << suffixQGram[pos];
	        }
	        std::cout << endl;
#endif
	        qPos = minLength - 1;
	        for (int pos = ePos; pos > bPos && qPos >= 0; ) {
	        	int matchLength = 0;

	        	while (ttsFiber[pos--] == suffixQGram[qPos--]) {
//	        		std::cout << ttsFiber[pos+1] << "--" << suffixQGram[qPos+1] << " ( "
//	        				<< (bool)((char)ttsFiber[pos+1] == (char)suffixQGram[qPos+1]) << "); ";
	        		matchLength++;
	        	}

	        	if (matchLength > maxSeedLength) {
	        		maxSeedLength = matchLength;
	        		maxSeedEndPos = pos + matchLength;
	        	}
	        }
	        std::cout << endl;
	        std::cout
					<< "MaxSeedLength: " << maxSeedLength << std::endl
					<< "MaxSeedEndPos: " << maxSeedEndPos << std::endl;

	        // found a valid match, extend for each tfo which had this suffix
	        printf("Suffixes for endPos #%d: %d\n", loc, ePos);
	        for (; itSB != itEB; ++itSB) {
	            unsigned seqId      = itSB->i1;
	            unsigned suffixPos  = itSB->i2; // start position
	#ifdef DEBUG
	            printf("seqId: %d\tsuffix: %d - %d\n", seqId, suffixPos, suffixPos + minLength - 1);
	#endif
	        }
	    }
	}

	/**
	 * One parameter is the q-gram index of the TFO set. This function iterates over each double
	 * stranded sequence in the haystack, and for each of these sequences iterates over each substring
	 * s of length minLength (== len(q-gram)). For each s we do Myers' approximate matching with maximal
	 * k (# allowed errors), and then verify / extend these matches (seeds).
	 */
	template<
	typename THaystack,		// haystack spec - double stranded sequences (tts)
	typename TQGramIndex,	// q-gram index - (tfo)
	typename TQuerySet,		// query set (needle) - single stranded sequences (tfo)
	typename TError,		// error rate
	typename TSize,			// minimum hit size
	typename TDrop,			// xdrop
	typename TSpec,			// specialization
	typename TId,			// sequence id
	typename TWorker
	>
	void plantMyers(
			Gardener<TId, TSpec>&gardener,
			THaystack			&haystack, 	// TTS set
			TQGramIndex			&index,
			TQuerySet			&queries,	// TFO set
			TError const		&errorRate,
			TSize const			&minLength,
			TDrop const			&xDrop,
			TWorker
	){
	    double t = sysTime();
		typedef typename Value<Gardener<TId, TSpec> >::Type		THitMap;
		typedef typename Value<THitMap>::Type					THitMapEntry;
		typedef typename Value<THitMapEntry,2>::Type			THitSetPointer;
		typedef typename Value<THitSetPointer>::Type			THitSet;

	    typedef typename Value<THaystack>::Type               	THaystackValue;
	    typedef typename Fibre<TQGramIndex, QGramSA>::Type      TSA;
	    typedef typename Fibre<TQGramIndex, QGramSADir>::Type   TSADir;
	    typedef typename Value<TSADir>::Type                    TSADirValue;
	    typedef typename Iterator<TSA, Standard>::Type          TSAIter;
	    typedef typename Iterator<TSADir, Standard>::Type       TSADirIter;

	    typedef typename Fibre<TQGramIndex, QGramCounts>::Type      TCounts;
	    typedef typename Fibre<TQGramIndex, QGramCountsDir>::Type   TCountsDir;
	    typedef typename Iterator<TCounts, Standard>::Type      	TIterCounts;
	    typedef typename Iterator<TCountsDir, Standard>::Type   	TIterCountsDir;
		typedef typename Suffix<typename GetSequenceByNo<TQGramIndex const>::Type >::Type	TSuffix;


	    int score;
	    int numLocations;
	    int *endLocations;
	    int *startLocations;
	    int alignmentLength;
	    unsigned char *alignment;

	    int cnt = 0;
	    int k 	= ceil(errorRate * minLength);

	    // init SA iterators for the tfoSet SA index
	    TIterCountsDir itCountsDir      = begin(indexCountsDir(index), Standard());
	    TIterCountsDir itCountsDirEnd   = end(indexCountsDir(index), Standard());
	    TIterCounts    itCountsBegin    = begin(indexCounts(index), Standard());

	    TSADirIter  itBucketDir         = begin(indexDir(index), Standard());
	    TSADirIter  itBucketDirEnd      = end(indexDir(index), Standard());
	    TSAIter     saBegin             = begin(indexSA(index), Standard());

	    TSADirValue 	bucketBegin = *itBucketDir;	// ptr to each unique qgram in SA of tfoSet (each bucket
	    											// may have several entries
	    THaystackValue 	ttsFiber; 					// current item in haystack (tts)
	    TSuffix 		suffixQGram;				// q-gram of appropriate length within suffix
	    THitSetPointer 	hitsPointer;				//

	    // iterate over each tts in the set of double stranded sequences (haystack)
	    for (int i = 0; i < length(haystack); ++i) {
	    	ttsFiber = haystack[i]; // next fiber in haystack
	    	unsigned char  	bitTarget[length(ttsFiber)];
	    	unsigned int   	targetLength = length(ttsFiber);

	    	// rewrite double stranded seq (tts) in a proper form for Myers (char -> index)
	    	for (int i = 0; i < length(ttsFiber); ++i) {
	    		bitTarget[i] = static_cast<unsigned int>(ttsFiber[i]);
	    	}

	    	// iterate over each substring of length 'minLength' == same as iterating over each
	    	// bucket (unique q-gram) in the q-gram index
	    	for (++itBucketDir; itBucketDir != itBucketDirEnd; ++itBucketDir)
	    	{
	    		TSADirValue bucketEnd = *itBucketDir; // end of this bucket == beginning of next one
	    		if (bucketBegin != bucketEnd)
	    		{
	    			unsigned char bitQGram[minLength];
	    			hitsPointer = new THitSet;

	    			TSAIter itBucketItem = saBegin + bucketBegin;
	    			TSAIter itEndBucket  = saBegin + bucketEnd;
	    			cnt ++;

	    			suffixQGram = suffix(indexText(index), *itBucketItem);

	    			// transform query into index for Myers
	    			for (int i = 0; i < minLength; ++i) {
	    				bitQGram[i] = static_cast<unsigned int>(suffixQGram[i]);
	    			}

#ifndef DEBUG
	    			// print query and target
	    			printf("--------------------------\n");
	    			std::cout << "Current suffix: " << suffixQGram << std::endl;
	    			std::cout << "Tts fiber: " << ttsFiber << std::endl;
	    			for (int i = 0; i < minLength; ++i) {
//	    				printf("%d;%c  ", bitQGram[i], (char)suffixQGram[i]);//(unsigned int)(char)(suf[i]));
	    				printf("%c",(char)suffixQGram[i]);//(unsigned int)(char)(suf[i]));
	    			}
	    			printf("\n");
//	    			for (int i = 0; i < targetLength; ++i) {
//	    				printf("%c", (char)bitTarget[i]);
//	    			}
//	    			printf("\n");
#endif

	    			std::cout << "error rate: " << errorRate << "; k: " << k << std::endl;

	    			// calculate Myers distance - this yields putative matches in endLocations
	    			edlibCalcEditDistance(bitQGram, minLength, bitTarget, targetLength, 4,
	    					k, EDLIB_MODE_HW, true, true, &score,
	    					&endLocations, &startLocations, &numLocations,
	    					&alignment, &alignmentLength);

	    			verifyMatches(*hitsPointer, ttsFiber, queries, suffixQGram, itBucketItem, itEndBucket,
	    					errorRate, numLocations, endLocations, minLength, xDrop);
	    		}
	    		bucketBegin = bucketEnd;
	    	}
	    }

	    std::cout << cnt << std::endl;
	    std::cout << "elapsed time: "<< sysTime() - t << " sec\n";
	    exit(-1);
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
			   TQuerySet			&queries,   // tts(original), now tfoSet
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
        std::cout << "qux\n";
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
