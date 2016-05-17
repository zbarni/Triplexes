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

#ifndef FBUSKE_APPS_TRIPLEXATOR_HEADER_MYERS_H
#define FBUSKE_APPS_TRIPLEXATOR_HEADER_MYERS_H

#include <limits>
#include "triplex_alphabet.h"
#include "helper.h"
#include "edlib.h"
#include <seqan/seeds2.h>  // Include module under test.
#include <seqan/sequence/adapt_std_list.h>
#include <seqan/misc/priority_type_base.h>
#include <seqan/misc/priority_type_heap.h>
#include <seqan/misc/misc_dequeue.h>

#define DEBUG
#define TOLERATED_ERROR 2
using namespace seqan;

namespace SEQAN_NAMESPACE_MAIN
{    
	template<
	typename THaystackFiber,
	typename TQuery,
	typename TSeed,
	typename TError,
	typename TVector,
	typename TOptions
	>
	void getMismatchOffsets(
			TSeed 			const 	&seed,
			THaystackFiber 	const 	&fiber,
			TQuery 			const 	&query,
			TError  		const	&errorRate,
			TVector					&lMmOffsets,
			TVector					&rMmOffsets,
			TVector					&lGuanines,
			TVector					&rGuanines,
			TOptions		const &options
	)
	{
		int posFiber;
		int posQuery;
		int mismatches;
		int consecutiveMismatches;
		unsigned int guanine;
		unsigned int localK = 1000;
		int k = ceil(errorRate * options.minLength);

		bool isLeftZeroIncluded = false;
		bool isRightEndIncluded = false;

		// go left
		guanine 				= 0;
		mismatches 				= 0;
		consecutiveMismatches 	= 0;
		posFiber = getBeginDim0(seed);
		posQuery = getBeginDim1(seed);

		while (posFiber >= 0 && posQuery >= 0 && mismatches <= localK + TOLERATED_ERROR
				&& consecutiveMismatches <= options.maxInterruptions) {
			if (fiber[posFiber] != query[posQuery]) {
				// update guanine counts (left from seed, not including those within seed)
				lGuanines.push_back(guanine);
				lMmOffsets.push_back(getBeginDim0(seed) - posFiber);
				++mismatches;
				if (posFiber == 0 || posQuery == 0) {
					isLeftZeroIncluded = true;
				}
				++consecutiveMismatches;
			}
			else {
				consecutiveMismatches = 0;

				if (isGuanine(fiber[posFiber])) {
					++guanine;
				}
			}
			localK = std::max((int)floor((getEndDim0(seed) - posFiber + 1) * errorRate), k);
			--posFiber;
			--posQuery;
		}

		unsigned int lastGuanine = (lGuanines.size()) ? lGuanines.back() : 0;
		// add corresponding end points if required
		if ((posFiber == -1 || posQuery == -1) && !isLeftZeroIncluded) {
			if ((posFiber == -1 && isGuanine(fiber[0])) ||
				(posQuery == -1 && isGuanine(query[0])))
			{
				lGuanines.push_back(lastGuanine + 1);
			}
			else {
				lGuanines.push_back(lastGuanine);
			}

			lMmOffsets.push_back(std::min((int)(getBeginDim0(seed)),(int)(getBeginDim1(seed))));
		}

		// go right
		guanine 				= 0;
		mismatches 				= 0;
		consecutiveMismatches 	= 0;
		posFiber = getEndDim0(seed);
		posQuery = getEndDim1(seed);

		int endFiber = length(fiber);
		int endQuery = length(query);

		// only go until the second last element, the last one is added later if needed (end* - 1)
		while (posFiber < endFiber && posQuery < endQuery && mismatches <= localK + TOLERATED_ERROR
				&& consecutiveMismatches <= options.maxInterruptions) {
			if (fiber[posFiber] != query[posQuery]) {
				// update guanine counts (left from seed, not including those within seed)
				rGuanines.push_back(guanine);
				rMmOffsets.push_back(posFiber - getEndDim0(seed));
				mismatches++;
				if (posFiber == endFiber - 1 || posQuery == endQuery - 1) {
					isRightEndIncluded = true;
				}
				++consecutiveMismatches;
			}
			else {
				consecutiveMismatches = 0;

				if (isGuanine(fiber[posFiber])) {
					++guanine;
				}
			}
			localK = std::max((int)floor((posFiber - getBeginDim0(seed) + 1) * errorRate), k);
			++posFiber;
			++posQuery;
		}

		lastGuanine = (rGuanines.size()) ? rGuanines.back() : 0;
		// add corresponding end points if required
		if ((posFiber == endFiber || posQuery == endQuery) && !isRightEndIncluded) {
			if ((posFiber == endFiber && isGuanine(fiber[endFiber - 1])) ||
				(posQuery == endQuery && isGuanine(query[endQuery - 1])))
			{
				rGuanines.push_back(lastGuanine + 1);
			}
			else {
				rGuanines.push_back(lastGuanine);
			}

			rMmOffsets.push_back(std::min((int)(length(query) - getEndDim1(seed) - 1), (int)(length(fiber) - getEndDim0(seed) - 1)));
		}
	}

	template<
	typename TSeedHashMap,
	typename TSeed,
	typename TSeedList,
	typename THaystackFiber,
	typename TQuery,
	typename TError,
	typename TOptions
	>
	void extendSimpleSeed(
			TSeedHashMap		  &addedSeedHashes,
			TSeed 			const &seed,
			TSeedList			  &extendedSeeds,
			THaystackFiber 	const &fiber,
			TQuery 			const &query,
			TError  		const &errorRate,
			int				const &k,
			unsigned int	const &seedGuanine,
			TOptions		const &options)
	{
		extendedSeeds.clear();
		typedef std::vector<unsigned int> TVector;
		// keep list of where mismatches occur, left and right from seed
		TVector lMmOffsets;
		TVector rMmOffsets;
		TVector lGuanines;
		TVector rGuanines;

		getMismatchOffsets(seed, fiber, query, errorRate, lMmOffsets, rMmOffsets, lGuanines, rGuanines, options);

		int bDim0;
		int bDim1;
		int eDim0;
		int eDim1;
		int mismatches;
		int localK 			= 1000; // should be inf but this is enough
		int guanineRate;

		// ============================================================================
		// START REAL SHIT
		// init left and right mismatch numbers
		int lMmSize = lMmOffsets.size();
		int rMmSize = rMmOffsets.size();
		int previousR = -1;

#ifdef DEBUG
		cout << "lMmOffsets: " << endl << "\t";
		for (int i = 0; i < lMmSize; i++) {
			cout << lMmOffsets[i] << " ";
		}
		cout << endl << "rMmOffsets: " << endl << "\t";
		for (int i = 0; i < rMmSize; i++) {
			cout << rMmOffsets[i] << " ";
		}
		cout << endl << "lGuanines: " << endl << "\t";
		for (int i = 0; i < lGuanines.size(); i++) {
			cout << lGuanines[i] << " ";
		}
		cout << endl << "rGuanines: " << endl << "\t";
		for (int i = 0; i < rGuanines.size(); i++) {
			cout << rGuanines[i] << " ";
		}
		cout << endl << endl << std::flush;
#endif

		for (int l = lMmSize - 1; l >= 0; --l) {
			for (int r = rMmSize - 1; r >= 0; --r) {
				bDim0 = getBeginDim0(seed) - lMmOffsets[l];
				bDim1 = getBeginDim1(seed) - lMmOffsets[l];
				eDim0 = getEndDim0(seed) + rMmOffsets[r];
				eDim1 = getEndDim1(seed) + rMmOffsets[r];
				int cnt = 0;

				// skip beginning positions until match found
				while (fiber[bDim0] != query[bDim1]) {
					if (cnt) {
						--l;
					}
					++bDim0;
					++bDim1;
					++cnt;
				}

				cnt = 0;
				while (fiber[eDim0] != query[eDim1]) {
					if (cnt) {
						--r;
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
				cout << "previousR: " << previousR  << std::flush << endl;
#endif

				// make sure #mismatches is valid
				if (localK < mismatches) {
#ifdef DEBUG
					cout << "localK < mismatches" << endl;
#endif
					continue;
				}

				if (r <= previousR) {
#ifdef DEBUG
					cout << "r <= previousR: " << r << " <= " << previousR << endl;
#endif
					continue;
				}

				guanineRate = lGuanines[l] + rGuanines[r] + seedGuanine;
				if (guanineRate < ceil(diag * options.minGuanineRate)) {
#ifdef DEBUG
					TSeed s (bDim0, bDim1, eDim0 + 1, eDim1 + 1);
					cout << "Discarding extended seed due to low guanine rate: " << endl;
					cout << "tts : " << infix(fiber, bDim0, eDim0 + 1) << endl;
					cout << "tfo : " << infix(query, bDim1, eDim1 + 1) << endl;
					cout << "Guanine rate: " << guanineRate << " vs " << ceil(diag * options.minGuanineRate) << endl;
					cout << "l: " << l << "; r: " << r << endl;
					cout << "seedGuanine: " << seedGuanine << endl;
					cout << "---- ***** ------" << endl;
#endif
//					continue;
				}
				else {
#ifdef DEBUG
					cout << "Guanine rate should be okay: " << guanineRate << " vs " << floor(diag * options.minGuanineRate) << endl;
					cout << "l: " << l << "; r: " << r << endl;
					cout << "seedGuanine: " << seedGuanine << endl;
#endif
				}

				// eDim* still include the last match, to report we add +1
				if (eDim0 + 1 - bDim0 >= options.minLength) {
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
#ifdef DEBUG
						cout << "valid & new seed extension: " << extendedSeeds.back() << endl;
						cout << "tts : " << infix(fiber, bDim0, eDim0 + 1) << endl;
						cout << "tfo : " << infix(query, bDim1, eDim1 + 1) << endl;
#endif
					}
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
			std::cout
					<< "== XXX >> Seed is already in map: " << seed << std::endl
					<< "\thash: " << hash << std::endl
					<< "\ttfoSeqNo: " << tfoSeqNo << std::endl
					<< "\thaystackFiberSeqNo: " << haystackFiberSeqNo  << std::flush << std::endl;
	#endif
			return false;
		}
		addedSeedHashes.insert(hash);
	#ifdef DEBUG
		std::cout << "===>>> New seed is added to map: " << seed << std::endl
				<< "\thash: " << hash << std::endl
				<< "\ttfoSeqNo: " << tfoSeqNo << std::endl
				<< "\thaystackFiberSeqNo: " << haystackFiberSeqNo  << std::flush << std::endl;
	#endif
		return true;
	}


	template<typename TChar>
	bool isGuanine(TChar const &c) {
		return c == 'G';
	}

	/**
	 *
	 */
	template<
	typename TTimes,
	typename TSeedMap,
	typename THitList,
	typename THaystack,
	typename TMergedHaystack,
	typename TSegment,
	typename TPosToFiberSeqNo,
	typename TMotifSet,
	typename TSuffix,
	typename TSAIter,
	typename TError,
	typename TOptions,
	typename THit,
	typename THitListKey
	>
	void verifyMatches(
			TTimes					&times,
			TSeedMap	   			&addedSeedHashMap,
			TSeedMap	   			&initMaxSeedHashMap,
			THitList				&hitList,
			THaystack 		const 	&haystack,
			TMergedHaystack const 	&mergedHaystack,
			TSegment 		const 	&segmentMap,
			TPosToFiberSeqNo const 	&posToFiberSeqNo,
			TMotifSet 		const	&needleSet,
			TSuffix 		const	&suffixQGram,
			TSAIter 		const	&itStartBucket,
			TSAIter 		const	&itEndBucket,
			TError 			const   &errorRate,
			int 		   	const	&numLocations,
			int 			const	endLocations[],
			TOptions		const	&options,
			THit,
			THitListKey)
	{
		typedef int								TScore;
		typedef Seed<Simple, DefaultSeedConfig>	TSeed;
		typedef std::vector< Seed<Simple, DefaultSeedConfig> > 			 TSeedList;
		typedef std::vector< Seed<Simple, DefaultSeedConfig> >::iterator TSeedListIterator;

		double t;
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
		int k = ceil(errorRate * options.minLength); // #mismatches allowed
		unsigned int guanine;
		TSAIter itSB, itEB;
		TSeedList extendedSeeds;

		// iterate over all putative matches (end locations), find max seed and then extend
		for (int match = 0; match < numLocations; match++) {
	#ifdef DEBUG
	//	    	cout << endl << "Iteration #" << match << " =====================================" << endl;
	#endif
			ePos = endLocations[match]; 	// end of current putative match in duplex (mergedHaystack)
			bPos = ePos - options.minLength + 1; 	// beginning of --||--
			if (bPos < 0) {
				continue;
			}

			// reset variables
			guanine	= 0;
			misM 	= 0;
			qPos 	= options.minLength - 1;
			maxSeedLength 		= 0;
			maxSeedMergedEnd	= 0;
			maxSeedQGramEnd 	= 0;
			maxSeedFiberEnd		= -1;
			t = sysTime();
			haystackFiberSeqNo 	= posToFiberSeqNo[bPos];
			times["gethaystackfiberno"] += sysTime() - t;
			assert(haystackFiberSeqNo >= 0);

			// make sure the found match doesn't span 2 different fibers in the original haystack
			if ((bPos - segmentMap[haystackFiberSeqNo]) + options.minLength > length(haystack[haystackFiberSeqNo])) {
	//	        if (haystackFiberSeqNo != posToFiberSeqNo[bPos + minLength] || bPos + minLength >= posToFiberSeqNo.size()) {
	#ifdef DEBUG
				cout << "Discarding Myers match because it spans over different fibers in original haystack." << std::flush << endl;
	#endif
				continue;
			}

			// reset SA iterators
			itSB = itStartBucket;
			itEB = itEndBucket;

			/////////////////////////////
			// validate TODO check for # of consecutive mismatches
			t = sysTime();
			int consecutiveMismatches = 0;
			for (int pos = ePos; pos > bPos && misM <= k + 1; --pos) {
				if (mergedHaystack[pos] != suffixQGram[qPos--]) {
					++misM;
					++consecutiveMismatches;
					if (consecutiveMismatches > options.maxInterruptions) {
						continue;
					}
				}
				else {
					consecutiveMismatches = 0;
				}
			}
			times["consmm"] += sysTime() - t;
			// invalid alignment. Note the k + 1, we want to be forgiving here because an extension might still
			// lower the error rate as wanted, so don't discard here just yet (will later, if it's not good).
			if (misM > k + 1) {
	#ifdef DEBUG
				cout << "Discarding Myers match because there are more than allowed mismatches."  << std::flush << endl;
	#endif
				continue;
			}

	#ifdef DEBUG
			cout << endl << "=====>>>>> Found valid match @ " << endl;
			cout << "Seed match is (Myers):" << endl << std::flush  << "\t";
			for (int pos = bPos; pos <= ePos; ++pos) {
				std::cout << mergedHaystack[pos];
			}
			std::cout << "\ttts" << endl << "\t";
			for (int pos = 0; pos < options.minLength; ++pos) {
				std::cout << suffixQGram[pos];
			}
			std::cout << "\ttfo"  << std::flush << endl;

	#endif
			// find the largest seed within the q-gram
			t 	 = sysTime();
			qPos = options.minLength - 1;
			for (int pos = ePos; pos > bPos && qPos >= 0; --pos, --qPos) {
				int tempGuanine = 0;
				int matchLength = 0;

				while (qPos >= 0 && mergedHaystack[pos] == suffixQGram[qPos]) {
					if (isGuanine(mergedHaystack[pos])) {
						++tempGuanine;
					}

					matchLength++;
					--pos;
					--qPos;
				}

				if (matchLength > maxSeedLength) {
					guanine				= tempGuanine;
					maxSeedLength 		= matchLength;
					maxSeedMergedEnd 	= pos  + matchLength;
					maxSeedQGramEnd		= qPos + matchLength;
					maxSeedFiberEnd		= maxSeedMergedEnd - segmentMap[haystackFiberSeqNo];
				}
			}

			if (maxSeedFiberEnd == -1 || !maxSeedLength) {
				cout << "BIG PROBLEM HERE, CAPTAIN! We have a seed of length 0?!!" << std::flush << endl;
			}
			times["maxseedfind"] += sysTime() - t;

	#ifdef DEBUG
			cout << endl << "Max. seed: " << endl  << std::flush << "\t";
			for (int pos = maxSeedMergedEnd - maxSeedLength + 1; pos <= maxSeedMergedEnd; ++pos) {
				cout << mergedHaystack[pos];
				assert(mergedHaystack[pos] == haystack[haystackFiberSeqNo][pos - segmentMap[haystackFiberSeqNo]]);
			}
			cout 	<< endl << "MaxSeed Length: " << maxSeedLength << std::endl
					<< "MaxSeed Position (mergedHaystack): " <<  maxSeedMergedEnd - maxSeedLength + 1
					<< " - " << maxSeedMergedEnd << std::endl << std::flush ;
	#endif

			// found a valid match, extend for each tfo which had this q-gram suffix
			for (; itSB != itEB; ++itSB) {
				int ndlSeqNo      	= itSB->i1;
				int qGramOffset  	= itSB->i2; // start position / offset of q-gram in tfo
				int qGramSeedBegin  = maxSeedQGramEnd - maxSeedLength + 1 + qGramOffset;
				int qGramSeedEnd 	= maxSeedQGramEnd + qGramOffset;

				THitListKey seqNoKey(haystackFiberSeqNo, ndlSeqNo);
				// check if seqNoKey key exists
				if (!initMaxSeedHashMap.count(seqNoKey)) {
					initMaxSeedHashMap[seqNoKey] = std::set<unsigned long long>();
				}
				if (!addedSeedHashMap.count(seqNoKey)) {
					addedSeedHashMap[seqNoKey] = std::set<unsigned long long>();
				}

				// calculate hash for maxSeed
				unsigned long long maxSeedHash = maxSeedFiberEnd - maxSeedLength + 1;
				maxSeedHash = (((((maxSeedHash << 16) + maxSeedFiberEnd) << 16) + qGramSeedBegin) << 16) + qGramSeedEnd;

				if (initMaxSeedHashMap[seqNoKey].count(maxSeedHash)) {
					continue;
				}

				initMaxSeedHashMap[seqNoKey].insert(maxSeedHash);

				TSeed seed(maxSeedFiberEnd - maxSeedLength + 1, qGramSeedBegin, maxSeedFiberEnd, qGramSeedEnd);
	#ifdef DEBUG
				std::cout << endl << "Original seed: " << seed << endl
						<< "Actual right ends (including this pos.): " << getEndDim0(seed) << ", " << getEndDim1(seed) << endl
						<< "seedFiber: " << infix(haystack[haystackFiberSeqNo], getBeginDim0(seed), getEndDim0(seed) + 1) << "\n"
						<< "seedQuery: " << infix(needleSet[ndlSeqNo], getBeginDim1(seed), getEndDim1(seed) + 1) << "\n" << std::flush ;
	#endif

				t = sysTime();
				extendSimpleSeed(addedSeedHashMap[seqNoKey], seed, extendedSeeds, haystack[haystackFiberSeqNo],
						needleSet[ndlSeqNo], errorRate, k, guanine, options);
				times["seedextend"] += sysTime() - t;

				for (TSeedListIterator seed = extendedSeeds.begin(); seed != extendedSeeds.end(); ++seed) {
	#ifdef DEBUG
					cout << "++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
					cout << "one seed after extension: " << *seed << endl;
					cout << "Match: " << endl;
					cout << "\tseedFiber: " << infix(haystack[haystackFiberSeqNo], getBeginDim0(*seed), getEndDim0(*seed)) << endl << std::flush ;
					cout << "\tseedQuery: " << infix(needleSet[ndlSeqNo], getBeginDim1(*seed), getEndDim1(*seed)) << endl;
					cout << "needle (tfo - isparallel): " << isParallel(needleSet[ndlSeqNo]) << "\t\t\t" << needleSet[ndlSeqNo] << endl;
					cout << "match length: " << getEndDim0(*seed) - getBeginDim0(*seed) << endl;
					cout << "tfo sequence id (index): \t\t" << ndlSeqNo << endl;
					cout << "+++++ ------- +++++++ ------ ++++++ will try to add new? seed" << endl << std::flush ;
	#endif
					// if new seed, add to seedMap and to hitSet
					t = sysTime();
					if (addIfNewSeed(haystackFiberSeqNo, ndlSeqNo, *seed, addedSeedHashMap[seqNoKey])) {
						times["addifnewseed"] += sysTime() - t;
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
	#ifdef DEBUG
						cout << "+++ adding new hit for matchKey<int,int> : " << matchKey.first << ", " << matchKey.second << endl << std::flush ;
	#endif
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
	typename TSegmentMap,
	typename TPosToFiberSeqNo
	>
	void mergeHaystack(
			THaystack const 	&haystack,
			TMergedHaystack 	&mergedHaystack,
			TSegmentMap  		&segmentMap,
			TPosToFiberSeqNo	&posToFiberSeqNo,
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
		std::cout 	<< "Merged haystack has length: " << totalLength << std::endl << std::flush ;
		std::cout 	<< "Fibers in merged haystack: " << std::endl << std::flush ;
	#endif

		// resize merged haystack to totalLength
		resize(mergedHaystack, totalLength, Exact());
		segmentMap.reserve(totalLength);
		bitTarget = new unsigned char[totalLength];
		totalLength = 0;

		for (int i = 0; i < length(haystack); ++i) {
			segmentMap.push_back(totalLength);
	#ifdef DEBUG
			std::cout 	<< "seq no: " << i << " (len = " << length(haystack[i]) << "): " << std::flush ;
	#endif
			// fill mergedHaystack and Myers vector
			for (int j = 0; j < length(haystack[i]); ++j) {
				mergedHaystack[totalLength + j] = haystack[i][j];
				// rewrite double stranded seq (tts) in a proper form for Myers (char -> index)
				bitTarget[totalLength + j] = charToBit(haystack[i][j]);//static_cast<unsigned int>(haystack[i][j]);

				// add starting position in mergedHaystack for current fiber
				posToFiberSeqNo.push_back(i);
			}
			totalLength += length(haystack[i]);
	#ifdef DEBUG
			std::cout 	<< haystack[i] << std::endl << std::flush ;
	#endif
		}

	#ifdef DEBUG
	//		std::cout << "Segment map / vector:\t";
	//		for (unsigned i = 0; i < segmentMap.size(); ++i) {
	//			cout << segmentMap[i] << " ";
	//		}
	//		cout << endl;
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
		typename THaystack,
		typename TNeedles,
		typename TError,
		typename THit
	>
	void mergeOverlappingHits(
			THitList 					&hitList,
			THitSetPointerMap 			&hitSetMap,
			THaystack 			const 	&haystack,
			TNeedles  			const	&needles,
			TError				const	&errorRate,
			int 				const	&minLength,
			THit)
	{
		typedef std::pair<int,int>										THitListKey;
		typedef typename 	THitList::iterator 							THitListIterator;
		typedef typename 	std::vector<THit*>::iterator 				THitIterator;
		typedef typename 	std::map<THitListKey, std::vector<THit*> >::iterator TDim0Iterator;

		int haystackFiberSeqNo;
		int ndlSeqNo;

		TDim0Iterator dim0It;
		TDim0Iterator nextDim0It;
	#ifdef DEBUG
		cout << endl << endl << "Merging overlaps" << endl;
	#endif

		for (THitListIterator it = hitList.begin(); it != hitList.end(); ++it) {
			haystackFiberSeqNo 	= (it->first).first;
			ndlSeqNo			= (it->first).second;
	#ifdef DEBUG
			cout << "Starting new tts-tfo pair" << endl;
	#endif

			for (dim0It = (it->second).begin(); dim0It != (it->second).end(); ++dim0It) {
				int currBegDim0 = dim0It->first.first;
				int currEndDim0 = dim0It->first.second;

				nextDim0It = dim0It;
				bool noMoreOverlap = false;
				for (++nextDim0It; nextDim0It != (it->second).end() && !noMoreOverlap; ++nextDim0It) {

					int nextBegDim0 = nextDim0It->first.first;
					int nextEndDim0 = nextDim0It->first.second;
	#ifdef DEBUG
					cout << endl << endl << "current matchKey<int,int> : " << currBegDim0 << ", " << currEndDim0 << endl << std::flush;
					cout << "next matchKey<int,int> : " << nextBegDim0 << ", " << nextEndDim0 << endl << std::flush;
	#endif
					// check if dim0 indices overlap
					if (nextBegDim0 >= currBegDim0 && nextBegDim0 < currEndDim0) {
	#ifdef DEBUG
						cout << "matches overlap, will iterate over all hit pairs..." << endl;
	#endif
						// iterate over all hit matches
						for (THitIterator currHitIt = dim0It->second.begin(); currHitIt != dim0It->second.end();) {
							bool incrementCurrHitIt = true;
							for (THitIterator nextHitIt = nextDim0It->second.begin(); nextHitIt != nextDim0It->second.end();) {
								THit *currHit = *currHitIt;
								THit *nextHit = *nextHitIt;
	#ifdef DEBUG
								cout 	<< "c Hit seed: = " << currBegDim0 << ", " << currBegDim0 + currHit->hitLength << ", - " << currHit->ndlPos << ", " << currHit->ndlPos + currHit->hitLength << endl
										<< "n Hit seed: = " << nextHit->hstkPos << ", " << nextHit->hstkPos + nextHit->hitLength << ", - " << nextHit->ndlPos << ", " << nextHit->ndlPos + nextHit->hitLength << endl
										<< "c Hit: hstkPos - ndlPos = " << currHit->hstkPos - currHit->ndlPos << endl
										<< "n Hit: hstkPos - ndlPos = " << nextHit->hstkPos - nextHit->ndlPos << endl;
	#endif
								if (currHit->hstkPos - currHit->ndlPos != nextHit->hstkPos - nextHit->ndlPos) {
									++nextHitIt;
	#ifdef DEBUG
									cout << "hits overlap but they are not aligned" << endl;
	#endif
									continue;
								}

								// case I: currHit contains nextHit
								if (currHit->ndlPos <= nextHit->ndlPos
										&& currHit->ndlPos + currHit->hitLength >= nextHit->ndlPos + nextHit->hitLength)
								{
	#ifdef DEBUG
									cout << "current hit contains next hit, delete next" << std::flush;
	#endif
									delete nextHit;
									nextHitIt = nextDim0It->second.erase(nextHitIt);
								}

								// case II: nextHit contains currHit
								else if (currHit->ndlPos >= nextHit->ndlPos
										&& currHit->ndlPos + currHit->hitLength <= nextHit->ndlPos + nextHit->hitLength
								) {
	#ifdef DEBUG
									cout << "next hit contains current hit, just delete current and skip to next current" << std::flush;
	#endif
									delete currHit;
									currHitIt = dim0It->second.erase(currHitIt);
									incrementCurrHitIt = false;
									break;
								}
								//  case III: hits overlap and are aligned, check if merge can be done
								else {
									int mismatches  = 0;
									int mergeEndPos = nextHit->hstkPos + nextHit->hitLength;
									int hstkPos = currHit->hstkPos;
									int ndlPos 	= currHit->ndlPos;
									int overlappedLength = mergeEndPos - currHit->hstkPos;
	#ifdef DEBUG
									cout << "trying to merge, overlapped segment len:" << overlappedLength << endl;
									if (mergeEndPos <= currHit->hstkPos + currHit->hitLength) {
										cout << "AJJJAJAJ wrong assumption while overlapping!" << endl;
									}
	#endif
									while (hstkPos < mergeEndPos) {
										if (haystack[haystackFiberSeqNo][hstkPos] != needles[ndlSeqNo][ndlPos]) {
											++mismatches;
										}
										++hstkPos;
										++ndlPos;
									}
									// should merge overlapping segments, errorRate isn't changed
									if (mismatches <= std::floor(errorRate * overlappedLength)) {
	#ifdef DEBUG
										cout << "okay, we can overlap these two..." << endl;
	#endif
										currHit->hitLength = overlappedLength;
										// delete next hit as it's now part of current hit
										delete nextHit;
										nextHitIt = nextDim0It->second.erase(nextHitIt);
									}
									else {
	#ifdef DEBUG
										cout << "damn mothafucka, we can NOT overlap these two because of the error rate?!" << endl;
										cout << "\tmismatches: " << mismatches << endl;
										cout << "\terror rate: " << errorRate << endl;
										cout << "\terrorRate * overlappedLength: " << std::floor(errorRate * overlappedLength) << endl;
	#endif
										++nextHitIt;
									}
								}
							}
							if (incrementCurrHitIt) {
								++currHitIt;
							}
						}
					}
					// no overlap --> add all hits and stop verifying current dim0It
					else {
	#ifdef DEBUG
						cout << "no overlap! stop after this" << endl;
	#endif
						noMoreOverlap = true;
					}
				}
	#ifdef DEBUG
				cout << "will add all remaining current hits, nextDim0It loop exited (or didn't enter)" << endl;
	#endif
				// add all remaining hits and stop verifying current dim0It
				for (THitIterator hit = dim0It->second.begin(); hit != dim0It->second.end(); ++hit) {
	#ifdef DEBUG
					cout 	<< "seed bDim0: " << (*hit)->hstkPos << endl
							<< "seed eDim0: " << (*hit)->hstkPos + (*hit)->hitLength << endl
							<< "seed bDim1: " << (*hit)->ndlPos << endl
							<< "seed eDim1: " << (*hit)->ndlPos + (*hit)->hitLength << endl;
	#endif
					add((*hitSetMap[haystackFiberSeqNo]), **hit);
				}
			} // end dim0It for loop
		} // end tts-tfo match for loop
	}


} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef FBUSKE_APPS_TRIPLEXATOR_HEADER_GARDENER_H
