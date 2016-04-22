#include <iostream>
#include <sstream>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/index.h>    
#include <seqan/find.h>    
#include <seqan/seeds.h>    
#include <seqan/file.h>   // Required to print strings in tests.
#include "myers/edlib.h"
#include <cmath>

#define QGRAM 19
//#define DEBUG

using namespace seqan;

template<
typename TTts, 
         typename TMotifSet,
         typename TSuffix,
         typename TSAIter>

         void verifyMatches(
                 TTts tts, 
                 TMotifSet tfoSet,
                 TSuffix suffix, 
                 const TSAIter itStartBucket,
                 const TSAIter itEndBucket,
                 int k,
                 int numLocations, 
                 int endLocations[]) {

             int ePos;
             int bPos;
             int qPos;
             int misM;
             TSAIter itSB, itEB;

             // iterate over all putative matches (end locations)
             // and start searching from back to start
             for (int loc = 0; loc < numLocations; loc++) {
                 ePos = endLocations[loc]; // end of current putative match
                 bPos = ePos - QGRAM;    // beginning of --||--
                 misM = 0;               // nr of mismatches
                 qPos = QGRAM - 1;

                 // reset SA iterators
                 itSB = itStartBucket;
                 itEB = itEndBucket;

                 for (int pos = ePos; pos > bPos && misM <= k; --pos) {
                     if (tts[pos] != suffix[qPos--]) {
                         ++misM;
                     }
                 }

                 // invalid alignment
                 if (misM > k) {
                     continue;
                 }

                 // found a valid match, extend for each tfo which had this suffix
                 //        printf("Suffixes for endPos #%d: %d\n", loc, ePos);
                 for (; itSB != itEB; ++itSB) {
                     //            unsigned seqId      = itSB->i1;
                     //            unsigned suffixPos  = itSB->i2;
#ifdef DEBUG
                     printf("seqId: %d\tsuffix: %d - %d\n", seqId, suffixPos, suffixPos + QGRAM - 1);
#endif
                 }
             }
         }

////////////////////////////////////////////
template<typename TMotifSet, typename TttsSet>
void align(TMotifSet tfoSet, TttsSet ttsSet) {

    double t = sysTime();

    typedef Index<TMotifSet, IndexQGram<UngappedShape<QGRAM>, OpenAddressing> > TQGramIndex;
    typedef typename Value<TttsSet>::Type               TTts;
    typedef Shape<Dna5, UngappedShape<QGRAM> >          TShape;
    typedef typename Fibre<TQGramIndex, QGramSA>::Type      TSA;
    typedef typename Fibre<TQGramIndex, QGramSADir>::Type   TSADir;
    typedef typename Value<TSADir>::Type                    TSADirValue;
    typedef typename Iterator<TSA, Standard>::Type          TSAIter;
    typedef typename Iterator<TSADir, Standard>::Type       TSADirIter;

    typedef typename Fibre<TQGramIndex, QGramCounts>::Type       TCounts;
    typedef typename Fibre<TQGramIndex, QGramCountsDir>::Type    TCountsDir;
    typedef typename Iterator<TCounts, Standard>::Type      TIterCounts;
    typedef typename Iterator<TCountsDir, Standard>::Type   TIterCountsDir;
    typedef typename seqan::Suffix<seqan::String<Dna5 > >::Type TSuffix;

    int score;
    int numLocations;
    int *endLocations;
    int *startLocations;
    int alignmentLength;
    unsigned char *alignment; 

    TTts tts = ttsSet[0]; // for now, only 1 row / sequence in tts
    unsigned char target[length(tts)];
    unsigned int targetLength = length(tts);

    for (int i = 0; i < length(tts); ++i) {
        target[i] = static_cast<unsigned int>(tts[i]);
    }

    TShape      shape;
    TQGramIndex index(tfoSet);  
    indexRequire(index, QGramCounts());
    indexRequire(index, QGramSADir());

    int cnt = 0;

    // initialize distance matrix
    //    int seqNum = countSequences(index);

    // init SA iterators
    TIterCountsDir itCountsDir      = begin(indexCountsDir(index), Standard());
    TIterCountsDir itCountsDirEnd   = end(indexCountsDir(index), Standard());
    TIterCounts itCountsBegin       = begin(indexCounts(index), Standard());

    TSADirIter  itSADir             = begin(indexDir(index), Standard());
    TSADirIter  itSADirEnd          = end(indexDir(index), Standard());
    TSAIter     saBegin             = begin(indexSA(index), Standard()); 

    std::cout << "index size: " << length(index) << std::endl;
    TSADirValue bucketBegin = *itSADir;
    for (++itSADir; itSADir != itSADirEnd; ++itSADir)
    {
        TSADirValue bucketEnd = *itSADir; // end of this bucket == beginning of next one
        if (bucketBegin != bucketEnd)
        {
            TSAIter itBucket    = saBegin + bucketBegin;
            TSAIter itEndBucket = saBegin + bucketEnd;
            cnt ++;

            TSuffix suf = suffix(indexText(index), *itBucket);

            unsigned char query[QGRAM];
            unsigned int  queryLength = QGRAM;

            // transform query 
            for (int i = 0; i < QGRAM; ++i) {
                query[i] = static_cast<unsigned int>(suf[i]);
            }

#ifdef DEBUG
            // print query and target
            printf("--------------------------\n");
            for (int k = 0; k < queryLength; ++k) { 
                printf("%d;%c  ", query[k], (char)suf[k]);//(unsigned int)(char)(suf[k]));
            }
            printf("\n");
            for (int k = 0; k < targetLength; ++k) { 
                printf("%d ", target[k]);
            }
            printf("\n");
#endif

            int k = 1;
            edlibCalcEditDistance(query, queryLength, target, targetLength, 4,
                    k, EDLIB_MODE_HW, true, true, &score,
                    &endLocations, &startLocations, &numLocations,
                    &alignment, &alignmentLength);

            verifyMatches(tts, tfoSet, suf, itBucket, itEndBucket, 
                    k, numLocations, endLocations);
        }
        bucketBegin = bucketEnd;
    }
    std::cout << cnt << std::endl; 
    std::cout << "elapsed time: "<< sysTime() - t << " sec\n";
}
/*
   Y-Xfoo bar qux asd X X Z X X
   ggagaggagagaggg
 */
void seed() {
    CharString seqH = "GGGAGAGAGGAGAGG";
    CharString seqV = "GGGAGAGGGGAGAGG";
    Score<int, Simple> scoreMatrix(1, -8, std::numeric_limits<int>::max());

    Seed<Simple> seed(4, 2, 9, 7);
    std::cout << "original\n"
        << "seedH: " << infix(seqH, beginPositionH(seed),
                endPositionH(seed)) << "\n"
        << "seedV: " << infix(seqV, beginPositionV(seed),
                endPositionV(seed)) << "\n";

    extendSeed(seed, seqH, seqV, EXTEND_BOTH, scoreMatrix, 16, UnGappedXDrop());
    std::cout << "result\n"
        << "seedH: " << infix(seqH, beginPositionH(seed),
                endPositionH(seed)) << "\n"
        << "seedV: " << infix(seqV, beginPositionV(seed),
                endPositionV(seed)) << "\n"
        << seed << std::endl;

}

void call(char *a) {
    a = new char[10];
}

int main(int argc, char *argv[])
{
    seed();
    return 0;
    std::ifstream fTfo (argv[1]);
    std::ifstream fTts (argv[2]);

    typedef Dna5String       Ttts;
    typedef Dna5String       Ttfo;
    typedef StringSet<Ttfo>  TMotifSet;
    typedef StringSet<Ttts>  TttsSet;

    std::string s;
    Ttts        rTTS;
    Ttfo        rTFO;
    TMotifSet   tfoSet;
    TttsSet     ttsSet;

    ////////////////////////////////////////////
    // read files
    while (fTfo >> s) {
        rTFO = static_cast<Ttfo>(s);
        appendValue(tfoSet, rTFO);
    }

    std::string concatenated;
    while (fTts >> s) {
        concatenated.append(s);
    }
    appendValue(ttsSet, static_cast<Ttts>(concatenated));

    fTfo.close();
    fTts.close();

    align(tfoSet, ttsSet);

    return 0;
}

