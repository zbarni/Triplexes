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
#include "local_container.h"

//#define DEBUG
#define TOLERATED_ERROR 3
#define TOLERATED_SEED_ERROR 2 // temporary error to allow for potentially long matches to be explored
#define MAX_OFFSET 1 // 1 means they must be at least adjacent
using namespace seqan;

namespace SEQAN_NAMESPACE_MAIN
{    
	typedef struct SeedDimStruct {
		int bDim0;
		int bDim1;
		int eDim0;
		int eDim1;

		SeedDimStruct() {
			bDim0 = bDim1 = eDim0 = eDim1 = 0;
		}

		SeedDimStruct(int p1, int p2, int p3, int p4) {
			bDim0 = p1;
			bDim1 = p2;
			eDim0 = p3;
			eDim1 = p4;
		}
	} SeedDimStruct;

	ostream& operator << (ostream& os, const SeedDimStruct& s) {
		return os << "\n\tbDim0, eDim0: " << s.bDim0 << ", " << s.eDim0 << endl
					<< "\tbDim1, eDim1: " << s.bDim1 << ", " << s.eDim1;
	}


	/**
	 * Returns true if parameter char is guanine, false otherwise.
	 */
	template<typename TChar>
	bool isGuanine(TChar const &c) {
		return c == 'G';
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

	/**
	 *
	 */
	template<
	typename THaystackFiber,
	typename TQuery,
	typename TSeed,
	typename TError,
	typename TVector,
	typename TOptions
	>
	void calcMismatchOffsets(
			TSeed 			const 	&seed,
			THaystackFiber 	const 	&fiber,
			TQuery 			const 	&needle,
			TError  		const	&errorRate,
			TVector					&lMmOffsets,
			TVector					&rMmOffsets,
			TVector					&lGuanines,
			TVector					&rGuanines,
			TOptions		const 	&options)
	{
		int posFiber;
		int posQuery;
		int mismatches;
		int consecutiveMm;
		unsigned int guanine;
		unsigned int localK = 1000;
		int k = floor(errorRate * options.minLength);

		bool isLeftZeroIncluded = false;
		bool isRightEndIncluded = false;

		// go left
		guanine 		= 0;
		mismatches 		= 0;
		consecutiveMm 	= 0;
		posFiber = getBeginDim0(seed) - 1;
		posQuery = getBeginDim1(seed) - 1;

		while (posFiber >= 0 && posQuery >= 0 && mismatches <= localK + TOLERATED_ERROR
				&& consecutiveMm <= options.maxInterruptions) {
			if (fiber[posFiber] != needle[posQuery]) {
				// update guanine counts (left from seed, not including those within seed)
				lGuanines.push_back(guanine);
				lMmOffsets.push_back(getBeginDim0(seed) - posFiber);
				++mismatches;
				if (posFiber == 0 || posQuery == 0) {
					isLeftZeroIncluded = true;
				}
				++consecutiveMm;
			}
			else {
				consecutiveMm = 0;

				if (isGuanine(fiber[posFiber])) {
					++guanine;
				}
			}
			localK = std::max((int)floor((getEndDim0(seed) - posFiber + 1) * errorRate), k);
			--posFiber;
			--posQuery;
		}

		// add corresponding end points if required
		if ((posFiber == -1 || posQuery == -1) && !isLeftZeroIncluded) {
			lGuanines.push_back(guanine);
			lMmOffsets.push_back(std::min((int)(getBeginDim0(seed)),(int)(getBeginDim1(seed))));
		}

		// go right
		guanine 				= 0;
		mismatches 				= 0;
		consecutiveMm 	= 0;
		posFiber = getEndDim0(seed) + 1;
		posQuery = getEndDim1(seed) + 1;

		const int endFiber  = length(fiber);
		const int endNeedle = length(needle);

		// only go until the second last element, the last one is added later if needed (end* - 1)
		while (posFiber < endFiber && posQuery < endNeedle && mismatches <= localK + TOLERATED_ERROR
				&& consecutiveMm <= options.maxInterruptions) {
			if (fiber[posFiber] != needle[posQuery]) {
				// update guanine counts (left from seed, not including those within seed)
				rGuanines.push_back(guanine);
				rMmOffsets.push_back(posFiber - getEndDim0(seed));
				mismatches++;
				if (posFiber == endFiber - 1 || posQuery == endNeedle - 1) {
					isRightEndIncluded = true;
				}
				++consecutiveMm;
			}
			else {
				consecutiveMm = 0;

				if (isGuanine(fiber[posFiber])) {
					++guanine;
				}
			}
			localK = std::max((int)floor((posFiber - getBeginDim0(seed) + 1) * errorRate), k);
			++posFiber;
			++posQuery;
		}

		// add corresponding end points if required
		if ((posFiber == endFiber || posQuery == endNeedle) && !isRightEndIncluded) {
			rGuanines.push_back(guanine);
			rMmOffsets.push_back(std::min((int)(endNeedle - getEndDim1(seed) - 1), (int)(endFiber - getEndDim0(seed) - 1)));
		}

#ifdef DEBUG
		cout << "lMmOffsets: " << endl << "\t";
		for (int i = 0; i < lMmOffsets.size(); i++) {
			cout << lMmOffsets[i] << " ";
		}
		cout << endl << "rMmOffsets: " << endl << "\t";
		for (int i = 0; i < rMmOffsets.size(); i++) {
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
	}

	/**
	 * Adjust initially extended seed so that the Guanine rate is satisfied. Returns all subsegments
	 * of the initial seed which satisfy the Guanine rate constraint, in guaninePotentials array.
	 * Returns true if any potential was found, false otherwise.
	 * WARNING!: all seed dimensions are contained in the respective seeds, unlike the end output!
	 */
	template<
	typename THaystackFiber,
	typename TQuery,
	typename TGuaninePotentials,
	typename TGuanine,
	typename TOptions,
	typename TDim
	>
	bool adjustForGuanineRate(
			THaystackFiber	const	&fiber,
			TQuery			const 	&needle,
			SeedDimStruct	const	&extSeedDim,
			TGuaninePotentials 		&guaninePotentials,
			TGuanine		const	&initGuanines,
			TOptions		const	&options,
			TDim)
	{
		guaninePotentials.clear();

		// if guanineRate already satisfied
		if (initGuanines >= ceil((extSeedDim.eDim0 - extSeedDim.bDim0 + 1) * options.minGuanineRate)) {
#ifdef DEBUG
			cout << "Guanine rate satisfied at the beginning already." << endl << std::flush;
#endif
			// add current (initial) seed and return
			guaninePotentials.push_back(SeedDimStruct(extSeedDim.bDim0, extSeedDim.bDim1, extSeedDim.eDim0, extSeedDim.eDim1));
			return true;
		}

		TDim maxDiag 	= -1;
		TDim diag	 	= extSeedDim.eDim0 - extSeedDim.bDim0 + 1;
		TDim tmp_bDim0 	= extSeedDim.bDim0;
		TDim tmp_bDim1 	= extSeedDim.bDim1;

#ifdef DEBUG
		cout << "Guanine adjustment started.." << endl << std::flush;
		cout << "initial extSeedDim: " << endl << extSeedDim << endl << std::flush;
		cout << "\ttts : " << infix(fiber, extSeedDim.bDim0, extSeedDim.eDim0 + 1) << endl << std::flush;
		cout << "\ttfo : " << infix(needle, extSeedDim.bDim1, extSeedDim.eDim1 + 1) << endl << std::flush;
		cout << "\tdiagonal including end points [REAL] : " << diag << endl << std::flush;
		cout << "\tInit guanine rate: " << initGuanines << " vs " << ceil(diag * options.minGuanineRate) << endl << std::flush;
#endif

		while (fiber[tmp_bDim0] == needle[tmp_bDim1] && extSeedDim.eDim0 > tmp_bDim0)
		{
			TDim tmp_eDim0 	= extSeedDim.eDim0;
			TDim tmp_eDim1 	= extSeedDim.eDim1;
			TDim tmp_diag	= tmp_eDim0 - tmp_bDim0 + 1;

			// adjust right end until needed
			while (fiber[tmp_eDim0 - 1] == needle[tmp_eDim1 - 1]// as long as we have matches only
					&& !isGuanine(fiber[tmp_eDim0])				// stop if 'G' reached
					&& tmp_eDim0 > tmp_bDim0					// to stay in valid index range
					&& initGuanines < ceil(tmp_diag * options.minGuanineRate))
			{
				--tmp_eDim0;
				--tmp_eDim1;
				--tmp_diag;
			}

			// if we reached a mismatch it's not a valid adjustment
			if (fiber[tmp_eDim0] != needle[tmp_eDim1]) {
				++tmp_bDim0;
				++tmp_bDim1;
				--diag;
				continue;
			}

			// if valid adjustment, add to guaninePotentials
			if (initGuanines >= ceil(tmp_diag * options.minGuanineRate)	&& tmp_diag >= maxDiag) {
				// if it's a new max, discard earlier ones
				if (tmp_diag > maxDiag) {
					guaninePotentials.clear();
				}
				maxDiag = tmp_diag;
				guaninePotentials.push_back(SeedDimStruct(tmp_bDim0, tmp_bDim1, tmp_bDim0 + maxDiag - 1, tmp_bDim1 + maxDiag - 1));
#ifdef DEBUG
				cout << "\tIndices have been adjusted for guanine rate." << endl << std::flush;
				cout << "\tAdded new guaninePotential: " << guaninePotentials.back() << endl << std::flush;
#endif
			}

			// adjust left end (+1) and diag to find all matches
			++tmp_bDim0;
			++tmp_bDim1;
			--diag;

			// break if previous was Guanine, since we'd skip it in next step
			if (isGuanine(fiber[tmp_bDim0 - 1])) {
				break;
			}
		}

		if (!guaninePotentials.size() && initGuanines >= ceil((extSeedDim.eDim0 - extSeedDim.bDim0 + 1) * options.minGuanineRate)) {
			guaninePotentials.push_back(SeedDimStruct(extSeedDim.bDim0, extSeedDim.bDim1, extSeedDim.eDim0, extSeedDim.eDim1));
		}

		// true if any potential was found
		return guaninePotentials.size() > 0;
	}

	/**
	 * Shifts a seed window to the right by 1 positions.
	 */
	void shiftSeedToRight(SeedDimStruct &extSeedDim) {
		++extSeedDim.bDim0;
		++extSeedDim.eDim0;
		++extSeedDim.bDim1;
		++extSeedDim.eDim1;
	}

	/**
	 * Tries to minimize the window size so that maxLength is satisfied. Returns true is this is
	 * possible, and false otherwise. Also, the leftmost such window is returned in extSeedDim.
	 */
	template<
	typename THaystackFiber,
	typename TNeedle,
	typename TDiag,
	typename TGuanine,
	typename TOptions
	>
	bool resizeWindowToMaxLength(
			SeedDimStruct 			&extSeedDim,
			THaystackFiber 	const 	&fiber,
			TNeedle			const	&needle,
			TDiag					&diag,
			TGuanine				&guanineRate,
			TOptions 		const	&options)
	{
		// if no limit on max length, return here
		if (options.maxLength == -1) {
			return true;
		}

		// find max diag that satisfies maxLength by reducing the right end
		while (diag > options.maxLength && fiber[extSeedDim.eDim0 - 1] == needle[extSeedDim.eDim1 - 1]) {
			// adjust guanineRate if skipping a 'G'
			if (isGuanine(fiber[extSeedDim.eDim0])) {
				--guanineRate;
			}
			--extSeedDim.eDim0;
			--extSeedDim.eDim1;
			--diag;
#ifdef DEBUG
			cout 	<< "options.maxLength is violated, continue to reduce right end" << endl
					<< "\tdiag:" << diag << endl
					<< "\textSeedDim:" << endl << extSeedDim << endl << std::flush;
#endif
		}

		// couldn't find valid diag from the right (reached a mismatch before constraint could be satisfied), search from left;
		if (diag > options.maxLength) {
#ifdef DEBUG
			cout << "options.maxLength is DEFINITELY right-violated, try left; length: " << diag << endl << std::flush;
#endif
			// find max diag that satisfies maxLength by reducing the right end
			while (diag > options.maxLength && fiber[extSeedDim.bDim0 + 1] == needle[extSeedDim.bDim1 + 1]) {
				// adjust guanineRate if skipping a 'G'
				if (isGuanine(fiber[extSeedDim.bDim0])) {
					--guanineRate;
				}
				++extSeedDim.bDim0;
				++extSeedDim.bDim1;
				--diag;
#ifdef DEBUG
				cout << "options.maxLength is violated, continue to increase left end" << endl
					<< "\tdiag:" << diag << endl
					<< "\textSeedDim:" << endl << extSeedDim << endl << std::flush;
#endif
			}

			// couldn't find valid diag from the left
			if (diag > options.maxLength) {
				return false;
			}
		}
		return true;
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
			TQuery 			const &needle,
			TError  		const &errorRate,
			int				const &k,
			unsigned int	const &seedGuanine,
			TOptions		const &options)
	{
		extendedSeeds.clear();
		typedef std::vector<unsigned int> 	TVector;
		typedef std::vector<SeedDimStruct> 	TSeedDimVector;

		// keep list of where mismatches occur, left and right from seed
		TVector lMmOffsets;
		TVector rMmOffsets;
		TVector lGuanines;
		TVector rGuanines;

		calcMismatchOffsets(seed, fiber, needle, errorRate, lMmOffsets, rMmOffsets, lGuanines, rGuanines, options);

		SeedDimStruct extSeedDim;
		int mismatches;
		int localK = 1000; // should be inf but this is enough
		int guanineRate;

		// ============================================================================
		// START REAL SHIT
		// init left and right mismatch numbers
		int lMmSize  	= lMmOffsets.size();
		int rMmSize  	= rMmOffsets.size();
		int fiberSize 	= length(fiber);
		int needleSize	= length(needle);

		for (int l = lMmSize - 1; l >= 0; --l) {
#ifdef DEBUG
			cout << "Started / lowered l to: " << l << endl << std::flush;
#endif
			int cnt = 0;
			extSeedDim.bDim0 = getBeginDim0(seed) - lMmOffsets[l];
			extSeedDim.bDim1 = getBeginDim1(seed) - lMmOffsets[l];
			// skip beginning positions until match found
			while (fiber[extSeedDim.bDim0] != needle[extSeedDim.bDim1]) {
				if (cnt) {
					--l;
				}
				++extSeedDim.bDim0;
				++extSeedDim.bDim1;
				++cnt;
			}

			TSeedDimVector guaninePotentials;
			const int init_bDim0 = extSeedDim.bDim0;
			const int init_bDim1 = extSeedDim.bDim1;
			for (int r = rMmSize - 1; r >= 0; --r) {
				extSeedDim.eDim0 = getEndDim0(seed) + rMmOffsets[r];
				extSeedDim.eDim1 = getEndDim1(seed) + rMmOffsets[r];
				// reset bdim
				extSeedDim.bDim0 = init_bDim0;
				extSeedDim.bDim1 = init_bDim1;

				cnt = 0;
				while (fiber[extSeedDim.eDim0] != needle[extSeedDim.eDim1]) {
					if (cnt) {
						--r;
					}
					--extSeedDim.eDim0;
					--extSeedDim.eDim1;
					++cnt;
				}

#ifdef DEBUG
				cout << "Testing potential match: " << endl << extSeedDim << endl << std::flush;
#endif

				// if diag is < minLength, no reason to continue
				int diag 	= extSeedDim.eDim0 - extSeedDim.bDim0 + 1;
				if (diag < options.minLength) {
#ifdef DEBUG
					cout << "potential extension too small already, continue:" << diag << endl << std::flush;
					cout << "l: " << l << ", r: " << r << endl << std::flush;
#endif
					continue;
				}

				guanineRate = lGuanines[l] + rGuanines[r] + seedGuanine;
				// resize window if too big
				if (!resizeWindowToMaxLength(extSeedDim, fiber, needle, diag, guanineRate, options)) {
					continue;
				}

				localK		= std::max((int)floor(diag * errorRate), k);
				mismatches 	= l + r;
				// make sure #mismatches is valid
				if (localK < mismatches) {
#ifdef DEBUG
					cout << "localK < mismatches, continue resizing window (--r): " << localK << "; " << mismatches << endl << std::flush;
#endif
					continue;
				}
#ifdef DEBUG
				else {
					cout 	<< "mismatch constraint satisfied, proceed. "
							<< "(localK, k, mismatches)" << localK << "; " << k << "; " << mismatches << endl << std::flush;
				}
#endif
				// valid diag window, now shift this window to the left and check for valid seeds at each window position
				// error rate must be checked again since diag size could have become smaller
				while (extSeedDim.eDim0 < fiberSize && extSeedDim.eDim1 < needleSize
						&& fiber[extSeedDim.bDim0] == needle[extSeedDim.bDim1] && fiber[extSeedDim.eDim0] == needle[extSeedDim.eDim1])
				{
#ifdef DEBUG
					cout << "shifting diag window over valid positions: " << diag << endl << std::flush;
#endif
					// check guanine rate
					bool guanineOK 	= adjustForGuanineRate(fiber, needle, extSeedDim, guaninePotentials, guanineRate, options, int());
					if (guanineOK)
					{
						// iterate over each guanine potentials
						for (TSeedDimVector::iterator seedDim = guaninePotentials.begin(); seedDim != guaninePotentials.end(); ++seedDim)
						{
#ifdef DEBUG
							cout << "Guanine rate should be okay: " << guanineRate << " vs " << ceil((seedDim->eDim0 - seedDim->bDim0) * options.minGuanineRate) << endl << std::flush;
							cout << "\tl: " << l << "; r: " << r << endl << std::flush;
							cout << "\tnew? tts : " << infix(fiber, seedDim->bDim0, seedDim->eDim0 + 1) << endl << std::flush;
							cout << "\tnew? tfo : " << infix(needle, seedDim->bDim1, seedDim->eDim1 + 1) << endl << std::flush;
							cout << "\tAfter guanine rate: " << guanineRate << " vs " << ceil((seedDim->eDim0 - seedDim->bDim0) * options.minGuanineRate) << endl << std::flush;
#endif
							// new diag after guanine adjustment
							unsigned tmp_diag = seedDim->eDim0 + 1 - seedDim->bDim0;
							// check again if error rate is still satisfied
							localK = std::max((int)floor(tmp_diag * errorRate), k);
							if (localK < mismatches) {
//								cout << "However, error rate is now not satisfied any more!" << endl;
								continue;
							}

							// seedDim->eDim* still include the last match, to report we add +1
							if (seedDim->eDim0 + 1 - seedDim->bDim0 >= options.minLength) {
								if (seedDim->bDim0 >= 65536 || seedDim->bDim1 >= 65536 || seedDim->eDim0 >= 65536 || seedDim->eDim0 >= 65536) {
									cout << "YOU SUCK! find a new hash function." << endl << std::flush;
									exit(-1);
								}
								long long hash = seedDim->bDim0;
								hash = (((((hash << 16) + seedDim->eDim0) << 16) + seedDim->bDim1) << 16) + seedDim->eDim1;

								// if new match
								if (!addedSeedHashes.count(hash)) {
									// add +1 to endDim*, it's how seeds / matches are reported
									extendedSeeds.push_back(TSeed(seedDim->bDim0, seedDim->bDim1, seedDim->eDim0 + 1, seedDim->eDim1 + 1));
									addedSeedHashes.insert(hash);
#ifdef DEBUG
									cout << "valid & new seed extension: " << extendedSeeds.back() << endl << std::flush;
									cout << "\ttts : " << infix(fiber, seedDim->bDim0, seedDim->eDim0 + 1) << endl << std::flush;
									cout << "\ttfo : " << infix(needle, seedDim->bDim1, seedDim->eDim1 + 1) << endl << std::flush;
#endif
								}
							} // end if
						} // end for over guaninePotentials
					}
#ifdef DEBUG
					else {
						cout << "Discarding extended seed due to low guanine rate, adjustment failed." << endl << std::flush;
					}
#endif
					// shift seed to the right
					shiftSeedToRight(extSeedDim);
					// update temp guanine numbers within current / updated seed
					if (extSeedDim.eDim0 < fiberSize && isGuanine(fiber[extSeedDim.eDim0])) {
						++guanineRate;
					}
					if (extSeedDim.bDim0 > 0 && isGuanine(fiber[extSeedDim.bDim0 - 1])) {
						--guanineRate;
					}
				} // end while window shift left
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
			return false;
		}
		addedSeedHashes.insert(hash);
	#ifdef DEBUG
		std::cout << "===>>> New seed is added to map: " << seed << std::endl << std::flush
				<< "\thash: " << hash << std::endl
				<< "\ttfoSeqNo: " << tfoSeqNo << std::endl
				<< "\thaystackFiberSeqNo: " << haystackFiberSeqNo  << std::flush << std::endl << std::flush;
	#endif
		return true;
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
		int k = floor(errorRate * options.minLength); // #mismatches allowed
		unsigned int guanine;
		TSAIter itSB, itEB;
		TSeedList extendedSeeds;

		// iterate over all putative matches (end locations), find max seed and then extend
		for (int matchId = 0; matchId < numLocations; matchId++) {
			ePos = endLocations[matchId]; 		// end of current putative match in duplex (mergedHaystack)
			bPos = ePos - options.minLength + 1;// beginning of --||--
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
	#ifdef DEBUG
				cout << "Discarding Myers match because it spans over different fibers in original haystack." << std::flush << endl;
	#endif
				continue;
			}

			// reset SA iterators
			itSB = itStartBucket;
			itEB = itEndBucket;

			// check for consecutive mismatch constraint
			t = sysTime();
			int consMismatches = 0;
			bool consSatisfied = true;
			for (int pos = ePos; pos > bPos && misM <= k + 1; --pos) {
				if (mergedHaystack[pos] != suffixQGram[qPos--]) {
					++misM;
					++consMismatches;
					if (consMismatches > options.maxInterruptions) {
						consSatisfied = false;
						break;
					}
				}
				else {
					consMismatches = 0;
				}
			}
			times["consmm"] += sysTime() - t;

			// invalid alignment. Note the k + 1, we want to be forgiving here because an extension might still
			// lower the error rate as wanted, so don't discard here just yet (will later, if it's not good).
			if (!consSatisfied || misM > k + TOLERATED_SEED_ERROR) {
#ifdef DEBUG
				cout << "Discarding Myers match because there are more consecutive mismatches than allowed / there are more than allowed mismatches"  << std::flush << endl;
				cout << "k: " << k << "; misMatches: " << misM << "; consSatisfied?: " << consSatisfied << endl << std::flush;
#endif
				continue;
			}

			/////////////////////////////////////////////////////////////////////////
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
				continue;
//				cout << "BIG PROBLEM HERE, CAPTAIN! We have a seed of length 0?!!" << std::flush << endl;
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

//				TODO @barni fix this
				if (!addedSeedHashMap.count(seqNoKey)) {
					addedSeedHashMap[seqNoKey] = std::set<unsigned long long>();
				}

				TSeed maxSeed(maxSeedFiberEnd - maxSeedLength + 1, qGramSeedBegin, maxSeedFiberEnd, qGramSeedEnd);

	#ifdef DEBUG
				std::cout << endl << "MAXSEED: " << maxSeed << endl
						<< "Actual right ends (including this pos.): " << getEndDim0(maxSeed) << ", " << getEndDim1(maxSeed) << endl
						<< "seedFiber: " << infix(haystack[haystackFiberSeqNo], getBeginDim0(maxSeed), getEndDim0(maxSeed) + 1) << "\n"
						<< "seedQuery: " << infix(needleSet[ndlSeqNo], getBeginDim1(maxSeed), getEndDim1(maxSeed) + 1) << "\n" << std::flush ;
	#endif

				t = sysTime();
				extendSimpleSeed(addedSeedHashMap[seqNoKey], maxSeed, extendedSeeds, haystack[haystackFiberSeqNo],
						needleSet[ndlSeqNo], errorRate, k, guanine, options);
				times["seedextend"] += sysTime() - t;

				// WARNING! from now on, the (extended) seeds have a +1 at their end dimensions as
				// they should for reporting results. E.g., an 10 - 21 extended seed has an actual
				// triplex length of 20 and lies between 10 - 20.
				for (TSeedListIterator seed = extendedSeeds.begin(); seed != extendedSeeds.end(); ++seed) {
	#ifdef DEBUG
					cout << "++++++++++++++++++++++++++++++++++++++++++++++++" << endl << std::flush;
					cout << "one seed after extension: " << *seed << endl;
					cout << "Match: " << endl << std::flush;
					cout << "\tseedFiber: " << infix(haystack[haystackFiberSeqNo], getBeginDim0(*seed), getEndDim0(*seed)) << endl << std::flush;
					cout << "\tseedQuery: " << infix(needleSet[ndlSeqNo], getBeginDim1(*seed), getEndDim1(*seed)) << endl << std::flush;
					cout << "match length: " << getEndDim0(*seed) - getBeginDim0(*seed) << endl << std::flush;
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
	 * Given a haystack with individual fibers, creates a single 'merged haystack' by concatenating
	 * all fibers one after the other. A map of which fiber starts where in the merged haystack
	 * is kept in segmentMap param.
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
			unsigned long 		&totalHaystackLength)
	{

		// calculate total length of segments in haystack
		for (int i = 0; i < length(haystack); ++i) {
			totalHaystackLength += length(haystack[i]);
		}

	#ifdef DEBUG
		std::cout 	<< "Merged haystack has length: " << totalHaystackLength << std::endl << std::flush ;
		std::cout 	<< "Fibers in merged haystack: " << length(haystack) << std::endl << std::flush ;
	#endif

		// resize merged haystack to totalLength
		resize(mergedHaystack, totalHaystackLength, Exact());
		segmentMap.reserve(totalHaystackLength);
		bitTarget = new unsigned char[totalHaystackLength];
		totalHaystackLength = 0;

		for (int i = 0; i < length(haystack); ++i) {
			segmentMap.push_back(totalHaystackLength);
			// fill mergedHaystack and Myers vector
			for (int j = 0; j < length(haystack[i]); ++j) {
				mergedHaystack[totalHaystackLength + j] = haystack[i][j];
				// rewrite double stranded seq (tts) in a proper form for Myers (char -> index)
				bitTarget[totalHaystackLength + j] = charToBit(haystack[i][j]);

				// add starting position in mergedHaystack for current fiber
				posToFiberSeqNo.push_back(i);
			}
			totalHaystackLength += length(haystack[i]);
		}
	}

	/**
	 * Computes guanine 'G' rate of a triplex, only matching guanines are considered valid.
	 */
	template<
	typename TFiber,
	typename TNeedle,
	typename TPos,
	typename TSize
	>
	int calcGuanines (
			TFiber	const 	&fiber,
			TNeedle	const	&needle,
			TPos			fiberPos,
			TPos			ndlPos,
			TSize	const	&mergedLength)
	{
		int guanines = 0;
		for (int i = 0; i < mergedLength; ++i) {
			// must be guanine and fiber and needle position must match
			if (isGuanine(fiber[fiberPos + i]) && fiber[fiberPos + i] == needle[ndlPos + i]) {
				++guanines;
			}
		}
#ifdef DEBUG
		cout << "#guanines in overlapped segment: " << guanines << endl << std::flush;
#endif
		return guanines;
	}

	/**
	 *
	 */
	template<
		typename THitSetPointerMap,
		typename THitList,
		typename THaystack,
		typename TNeedles,
		typename TError,
		typename TOptions,
		typename THit
	>
	void mergeOverlappingHits (
			THitList 					&hitList,
			THitSetPointerMap 			&hitSetMap,
			THaystack 			const 	&haystack,
			TNeedles  			const	&needles,
			TError				const	&errorRate,
			TOptions			const	&options,
			THit)
	{
		typedef 			std::pair<int,int>				THitListKey;
		typedef typename 	THitList::iterator 				THitListIterator;
		typedef typename 	std::vector<THit*>::iterator 	THitIterator;
		typedef typename 	std::map<THitListKey, std::vector<THit*> >::iterator TDim0Iterator;

		int haystackFiberSeqNo;
		int ndlSeqNo;

		TDim0Iterator dim0It;
		TDim0Iterator nextDim0It;

#ifdef DEBUG
		// TODO @barni remove
		cout << "Dumping (added) hitList [added and not yet merged seeds]:" << endl;
		for (THitListIterator it = hitList.begin(); it != hitList.end(); ++it) {
			for (dim0It = (it->second).begin(); dim0It != (it->second).end(); ++dim0It) {
				int currBegDim0 = dim0It->first.first;
				int currEndDim0 = dim0It->first.second;
				for (THitIterator currHitIt = dim0It->second.begin(); currHitIt != dim0It->second.end(); ++currHitIt) {
					THit *currHit = *currHitIt;
					cout 	<< "hit (fib vs ndl): " << currBegDim0 << " - " << currBegDim0 + currHit->hitLength << " vs. " << currHit->ndlPos << " - " << currHit->ndlPos + currHit->hitLength << endl;
				}
			}
		}
		cout << "===" << endl << "Dumping ended." << endl;
#endif

		for (THitListIterator it = hitList.begin(); it != hitList.end(); ++it) {
			haystackFiberSeqNo 	= (it->first).first;
			ndlSeqNo			= (it->first).second;
#ifdef DEBUG
			cout << "Starting merge of new tts-tfo pair" << endl;
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
					cout << "haystackFiberSeqNo: " << haystackFiberSeqNo << ", ndlSeqNo: " << ndlSeqNo << endl << std::flush;
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
#ifdef DEBUG
									cout << "hits overlap but they are not aligned" << endl;
#endif
									++nextHitIt;
									continue;
								}

								int mergeEndPos = nextHit->hstkPos + nextHit->hitLength;
								int overlappedLength = mergeEndPos - currHit->hstkPos;
								if (overlappedLength > options.maxLength && options.maxLength != -1) {
#ifdef DEBUG
									cout << "overlapped length longer than maxLength, ignoring possible merge" << endl;
#endif
									++nextHitIt;
									continue;
								}

								// case I: currHit contains nextHit
								if (currHit->ndlPos <= nextHit->ndlPos && currHit->ndlPos + currHit->hitLength >= nextHit->ndlPos + nextHit->hitLength)
								{
#ifdef DEBUG
									cout << "current hit contains next hit, delete next" << endl << std::flush;
#endif
									delete nextHit;
									nextHitIt = nextDim0It->second.erase(nextHitIt);
								}

								// case II: nextHit contains currHit
								else if (currHit->ndlPos >= nextHit->ndlPos
										&& currHit->ndlPos + currHit->hitLength <= nextHit->ndlPos + nextHit->hitLength)
								{
#ifdef DEBUG
									cout << "next hit contains current hit, just delete current and skip to next current" << endl << std::flush;
#endif
									delete currHit;
									currHitIt = dim0It->second.erase(currHitIt);
									incrementCurrHitIt = false;
									break;
								}
								//  case III: hits overlap and are aligned, check if merge can be done
								else {
									int mismatches  = 0;
									int hstkPos = currHit->hstkPos;
									int ndlPos 	= currHit->ndlPos;
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
									int guanines = calcGuanines(haystack[haystackFiberSeqNo], needles[ndlSeqNo], currHit->hstkPos, currHit->ndlPos, overlappedLength);
									// should merge overlapping segments, errorRate isn't changed
									if (mismatches <= std::floor(errorRate * overlappedLength) &&
											guanines >= ceil(overlappedLength * options.minGuanineRate))
									{
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
										cout << "damn mothafucka, we can NOT overlap these two; error rate / guanine rate vioalted?!" << endl;
										cout << "\tmismatches: " << mismatches << endl;
										cout << "\terror rate: " << errorRate << endl;
										cout << "\tguanines: " << guanines << endl;
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
					delete *hit;
				}
			} // end dim0It for loop
		} // end tts-tfo match for loop
	}



	/***************************************************************************************************
	 ***************************************************************************************************
	 *										PALINDROMIC STUFF
	 ***************************************************************************************************
	 ***************************************************************************************************
	 */

	/**
	 * Check if matching positions on fiber and needle are not? too far apart - we want locality!
	 * The distance is calculated by computing the absolute starting positions of the current
	 * seed (+/-, A/P considered) on the genome.
	 */
	template<typename TSeed, typename TFiber, typename TNeedle>
	bool checkLocalOverlapConstraints(
			bool  	const	plusStrand,
			bool  	const	needleParallel,
			TSeed 	const 	&seed,
			TFiber 	const 	&fiber,
			TNeedle const 	&needle)
	{
		typedef typename Position<TFiber>::Type TPos;
		int offset; // can be negative as well

		TPos fiberAbsBPos 	= beginPosition(fiber);
		TPos needleAbsBPos	= beginPosition(needle);

		fiberAbsBPos  += (plusStrand) ? getBeginDim0(seed) : length(fiber) - getEndDim0(seed);
		needleAbsBPos += (needleParallel) ? getBeginDim1(seed) : length(needle) - getEndDim1(seed);

		offset = fiberAbsBPos - needleAbsBPos;
#ifdef DEBUG
		cout << "\nseed: " << seed << endl;
		cout << "fiberAbsBPos: " << fiberAbsBPos << endl;
		cout << "needleAbsBPos: " << needleAbsBPos << endl;
		cout << "overlap offset: " << offset << ", returning: " << (abs(offset) <= (int)(getEndDim0(seed) - getBeginDim0(seed))) << endl;
#endif
		// TODO change max offset here
		return abs(offset) <= (int)(getEndDim1(seed) - getBeginDim1(seed) + MAX_OFFSET - 1);
	}

	/**
	 * Receives a list of putative matches (endLocations) on fiber and finds all maximal
	 * triplexes starting from these matches. Finds and extends seeds according to constraints, and
	 * all hits are appended to hitList.
	 */
	template<
	typename TTimes,
	typename TSeedMap,
	typename THitList,
	typename TPair,
	typename TFiber,
	typename TNeedle,
	typename TSeed,
	typename TError,
	typename TOptions,
	typename THit
	>
	void verifyLocalTriplexes(
			TTimes					&times,
			bool 		const	plusStrand,
			TSeedMap			&addedSeedHashMap,
			THitList			&hitList,
			TPair		const 	&seqNoKey,
			TFiber		const	&fiber,
			TNeedle 	const 	&needle,
			TSeed		const 	&needleSearchWindow,
			int 		const	&numLocations,
			int 		const	endLocations[],
			TError		const	&errorRate,
			TOptions	const	&options,
			THit)
	{
		typedef int								TScore;
		typedef std::vector< TSeed > 			TSeedList;
		typedef typename TSeedList::iterator 	TSeedListIterator;

		TSeedList extendedSeeds;

		int k = floor(options.errorRate * options.minLength); // #mismatches allowed

		// iterate over all putative matches (end locations), find max seed and then extend
		for (int match = 0; match < numLocations; match++) {
#ifdef DEBUG
			cout << "'verifyLocalTrilexes'. Found " << numLocations
					<< " possible match positions. Now doing fiber ePos: " << endLocations[match] << endl;
#endif

			int ePos = endLocations[match]; 			// end of current putative match in fiber
			int bPos = ePos - options.minLength + 1; 	// beginning of --||--

			// starting position of match is invalid
			if (bPos < 0) {
				continue;
			}

			// reset variables
			unsigned int guanine= 0;
			int misM 			= 0;	// counter for number of mismatches
			int maxSeedLength	= 0;	// length of maximum seed

			// check for allowed consecutive mismatches
			bool consSatisfied 	= true;
			int consMismatches 	= 0;
			int fibCurPos 		= 0; 	// current position in fiber (for comparison)
			int ndlCurPos 		= 0; 	// current position in needle (for comparison)

			for (fibCurPos = ePos, ndlCurPos = getEndDim1(needleSearchWindow); fibCurPos > bPos && misM <= k + 1; --fibCurPos, --ndlCurPos) {
				if (fiber[fibCurPos] != needle[ndlCurPos]) {
					++misM;
					++consMismatches;

					if (consMismatches > options.maxInterruptions) {
						consSatisfied = false;
						break;
					}
				}
				else {
					consMismatches = 0;
				}
			}

			// invalid alignment. Note the k + 1, we want to be forgiving here because an extension might still
			// lower the error rate as wanted, so don't discard here just yet (will later, if it's not good).
			if (!consSatisfied || misM > k + TOLERATED_SEED_ERROR) {
#ifdef DEBUG
				cout << "Discarding Myers match because there are more consecutive mismatches than allowed. // "
						<< "Discarding Myers match because there are more than allowed mismatches."  << std::flush << endl;
				cout << "k: " << k << "; misMatches: " << misM << "; consSatisfied?: " << consSatisfied << endl << std::flush;
#endif
				continue;
			}

			/////////////////////////////////////////////////////////////////////////
			// find the largest seed within the q-gram
			// current position in needle is end of seed window
			TSeed 	maxSeed;
			int 	needleWindowBPos = getBeginDim1(needleSearchWindow);
			for (fibCurPos = ePos, ndlCurPos = getEndDim1(needleSearchWindow); fibCurPos > bPos && ndlCurPos >= needleWindowBPos; --fibCurPos, --ndlCurPos)
			{
				int tempGuanine = 0;
				int matchLength = 0;
				while (ndlCurPos >= needleWindowBPos && fiber[fibCurPos] == needle[ndlCurPos]) {
					if (isGuanine(fiber[fibCurPos])) {
						++tempGuanine;
					}
					matchLength++;
					--fibCurPos;
					--ndlCurPos;
				}

				if (matchLength > maxSeedLength) {
					guanine	= tempGuanine;

					setBeginDim0(maxSeed, fibCurPos + 1);
					setBeginDim1(maxSeed, ndlCurPos + 1);
					setEndDim0(maxSeed, fibCurPos + matchLength);
					setEndDim1(maxSeed, ndlCurPos  + matchLength);

					maxSeedLength = matchLength;
				}
			}

			// make sure we have a valid maximum seed
			if (!maxSeedLength) {
				continue;
//				cout << "BIG PROBLEM HERE, CAPTAIN! We have a seed of length 0?!!" << std::flush << endl;
			}

			if (!addedSeedHashMap.count(seqNoKey)) {
				addedSeedHashMap[seqNoKey] = std::set<unsigned long long>();
			}
#ifdef DEBUG
			cout << "MAXSEED:" << maxSeed << endl << "\t"
					<< infix(fiber, getBeginDim0(maxSeed), getEndDim0(maxSeed) + 1) << endl << "\t"
					<< infix(needle,getBeginDim1(maxSeed), getEndDim1(maxSeed) + 1) << endl;
#endif

			// extend valid maximum seed
			double tt = sysTime();
			extendSimpleSeed(addedSeedHashMap[seqNoKey], maxSeed, extendedSeeds,
					fiber, needle, errorRate, k, guanine, options);
			times["seedextend"] += sysTime() - tt;

			for (TSeedListIterator seed = extendedSeeds.begin(); seed != extendedSeeds.end(); ++seed) {
#ifdef DEBUG
				cout << "one seed after extension: " << *seed << endl;
				cout << "Match: " << endl << std::flush;
				cout << "\tseedFiber: " << infix(fiber, getBeginDim0(*seed), getEndDim0(*seed)) << endl << std::flush ;
				cout << "\tseedQuery: " << infix(needle, getBeginDim1(*seed), getEndDim1(*seed)) << endl << std::flush;
				cout << "match length: " << getEndDim0(*seed) - getBeginDim0(*seed) << endl << std::flush;
				cout << "+++++ ------- +++++++ ------ ++++++ will try to add new? seed" << endl << std::flush ;
#endif
				// check if extended seed / match fits overlap (shift) constraints
				if (!checkLocalOverlapConstraints(plusStrand, isParallel(needle), *seed, fiber, needle))
				{
#ifdef DEBUG
					cout << "Overlap / offset constraint violated, discarding extended seed." << endl;
#endif
					continue;
				}

				// if new seed, add to seedMap and to hitSet
				if (addIfNewSeed(seqNoKey.first, seqNoKey.second, *seed, addedSeedHashMap[seqNoKey])) {
					THit *hit = new THit(
							seqNoKey.first,			// fiber seq no
							seqNoKey.second,		// needle seq no
							getBeginDim0(*seed),	// fiber beginning
							getBeginDim1(*seed),	// needle beginning
							-1, 					// the diagonal - NOT IMPORTANT
							0,
							getEndDim0(*seed)-getBeginDim0(*seed));

					TPair matchKey(getBeginDim0(*seed), getEndDim0(*seed));
					((hitList[seqNoKey])[matchKey]).push_back(hit);
				}
			}
		}
	}

	/**
	 * This function receives a fiber and a needle as input parameters, and computes all the matches
	 * between the two satisfying some maximum shift constraint.
	 */
	template<
	typename TTimes,
	typename TSeedHashMap,
	typename THitList,
	typename TFiber,
	typename TNeedle,
	typename TBitArray,
	typename TSeqNo,
	typename TOptions,
	typename TSeed,
	typename THit
	>
	void computeLocalTriplexes(
			TTimes				&times,
			TSeedHashMap 		&addedSeedHashMap,
			THitList			&hitList,
			TFiber 		const 	&fiber,
			TNeedle 	const 	&needle,
			TBitArray 			&bitFiber,			// bit representation of fiber, will be shifted
			TSeqNo		const	&fiberSeqNo,
			TSeqNo		const	&needleSeqNo,
			unsigned 			k,
			unsigned 			alphabetSize,
			bool		const 	&plusStrand,
			TOptions 	const 	&options,
			TSeed,
			THit)
	{
		typedef std::pair<TSeqNo, TSeqNo> TPair;
		TPair seqNoKey(fiberSeqNo, needleSeqNo);

		int score;
		int numLocations;
		int *endLocations;
		int *startLocations;
		int alignmentLength;
		unsigned char *alignment;

		// compute error rate
		double eR = options.errorRate;
		if (options.maximalError >= 0){
			eR = min(options.errorRate, max(double(options.maximalError)/options.minLength, 0.0));
		}

#ifdef DEBUG
		cout << "Overlap found: " << beginPosition(fiber) << " - " << endPosition(fiber) << " vs "
				<< beginPosition(needle) << " - " << endPosition(needle) << endl
				<< endl << "fiber (" << (plusStrand ? "+" : "-") << "): " << fiber << endl
				<< "nedle (" << isParallel(needle) << " = " << getMotif(needle) << "): " << needle << endl;
#endif
		unsigned fiberLength   = length(fiber);
		unsigned needleLength  = length(needle);
		unsigned char *bitNeedle = new unsigned char [needleLength];
		unsigned char * const bitNeedleBase = bitNeedle; // store original address for deletion later

		// init seed window and number of shifts
		TSeed needleSearchWindow(0, 0, 0, options.minLength - 1);

#ifdef DEBUG
		cout << "plusStrand?: " << plusStrand << endl;
		cout << "needleParallel?: " << isParallel(needle) << endl;
#endif

		// transform query into index for Myers
		for (int i = 0; i < needleLength; ++i) {
			bitNeedle[i] = charToBit((needle)[i]);
		}

		unsigned shifts = needleLength - options.minLength;
		for (unsigned i = 0; i <= shifts; ++i) {
#ifdef DEBUG
			cout << endl << "seedWindow (shifted): " << needleSearchWindow << endl;
			cout << "nedle window (search space): " << infix(needle, getBeginDim1(needleSearchWindow), getEndDim1(needleSearchWindow) + 1) << endl;
#endif

			//calculate Myers distance - this yields putative matches in endLocations
			numLocations = 0;
			double t = sysTime();
			edlibCalcEditDistance(bitNeedle + getBeginDim1(needleSearchWindow),
					options.minLength,
					bitFiber, //
					fiberLength, // size of current fiber window
					alphabetSize, k + TOLERATED_SEED_ERROR, EDLIB_MODE_HW, false, false, &score,
					&endLocations, &startLocations, &numLocations,
					&alignment, &alignmentLength);
    		times["myers"] += sysTime() - t;

    		t = sysTime();
			verifyLocalTriplexes(times, plusStrand, addedSeedHashMap, hitList, seqNoKey, fiber, needle,
					needleSearchWindow,	numLocations, endLocations, eR, options, THit());
    		times["verify"] += sysTime() - t;

			free(endLocations);

			// shift window to the right
			setBeginDim1(needleSearchWindow, getBeginDim1(needleSearchWindow) + 1);
			setEndDim1(needleSearchWindow, getEndDim1(needleSearchWindow) + 1);
		}

		delete [] bitNeedleBase;
	}

} //namespace SEQAN_NAMESPACE_MAIN

#endif
