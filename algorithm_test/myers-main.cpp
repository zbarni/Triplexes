#include <iostream>
#include <sstream>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/index.h>    
#include <seqan/find.h>    
#include <seqan/file.h>   // Required to print strings in tests.
#include "myers/edlib.h"
#include <cmath>

#define QGRAM 3

using namespace seqan;

void printAlignment(const unsigned char* query, const int queryLength,
                    const unsigned char* target, const int targetLength,
                    const unsigned char* alignment, const int alignmentLength,
                    const int position, const int modeCode, const char* idxToLetter) {
    int tIdx = -1;
    int qIdx = -1;
    if (modeCode == EDLIB_MODE_HW) {
        tIdx = position;
        for (int i = 0; i < alignmentLength; i++) {
            if (alignment[i] != 1)
                tIdx--;
        }
    }
    for (int start = 0; start < alignmentLength; start += 50) {
        // target
        printf("T: ");
        int startTIdx;
        for (int j = start; j < start + 50 && j < alignmentLength; j++) {
            if (alignment[j] == 1)
                printf("_");
            else
                printf("%c", idxToLetter[target[++tIdx]]);
            if (j == start)
                startTIdx = tIdx;
        }
        printf(" (%d - %d)\n", std::max(startTIdx, 0), tIdx);
        // query
        printf("Q: ");
        int startQIdx = qIdx;
        for (int j = start; j < start + 50 && j < alignmentLength; j++) {
            if (alignment[j] == 2)
                printf("_");
            else
                printf("%c", idxToLetter[query[++qIdx]]);
            if (j == start)
                startQIdx = qIdx;
        }
        printf(" (%d - %d)\n\n", std::max(startQIdx, 0), qIdx);
    }
}

////////////////////////////////////////////
template<typename TMotifSet, typename TttsSet>
void align(TMotifSet tfoSet, TttsSet ttsSet) {

    double t = sysTime();
    typedef typename Value<TMotifSet>::Type             TTfo;
    typedef typename Value<TttsSet>::Type               TTts;
    typedef Shape<Dna5, UngappedShape<QGRAM> >          TShape;
    typedef Index<TMotifSet, IndexQGram<UngappedShape<QGRAM>, OpenAddressing> > TQGramIndex;
//    typedef typename Value<TShape>::Type                TShapeValue;
    typedef typename Fibre<TQGramIndex, QGramSA>::Type      TSA;
    typedef typename Fibre<TQGramIndex, QGramSADir>::Type   TSADir;
    typedef typename Value<TSADir>::Type                    TSADirValue;
//    typedef typename Value<TSA>::Type                       TSAValue;
    typedef typename Iterator<TSA, Standard>::Type          TSAIter;
    typedef typename Iterator<TSADir, Standard>::Type       TSADirIter;

    typedef typename Fibre<TQGramIndex, QGramCounts>::Type       TCounts;
    typedef typename Fibre<TQGramIndex, QGramCountsDir>::Type    TCountsDir;
    typedef typename Value<TCountsDir>::Type                TDirValue;
    typedef typename Iterator<TCounts, Standard>::Type      TIterCounts;
    typedef typename Iterator<TCountsDir, Standard>::Type   TIterCountsDir;
    
    char idxToLetter[4] = {'A', 'C', 'G', 'T'};
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
    int seqNum = countSequences(index);

    TIterCountsDir itCountsDir      = begin(indexCountsDir(index), Standard());
    TIterCountsDir itCountsDirEnd   = end(indexCountsDir(index), Standard());
    TIterCounts itCountsBegin       = begin(indexCounts(index), Standard());

    TSADirIter  itSADir             = begin(indexDir(index), Standard());
    TSADirIter  itSADirEnd          = end(indexDir(index), Standard());
    TSAIter     saBegin             = begin(indexSA(index), Standard()); 


    TSADirValue bucketBegin = *itSADir;
    for (++itSADir; itSADir != itSADirEnd; ++itSADir)
    {
        TSADirValue bucketEnd = *itSADir; // end of this bucket == beginning of next one
        if (bucketBegin != bucketEnd)
        {
            TSAIter itA = saBegin + bucketBegin;
            TSAIter itEnd = saBegin + bucketEnd;
//            std::cout << *itA << std::endl;
            cnt ++;
    
            seqan::Suffix<seqan::String<Dna5 > >::Type suf = suffix(indexText(index), *itA);

            unsigned char query[QGRAM];
            unsigned int  queryLength = QGRAM;

            // transform query 
            for (int i = 0; i < QGRAM; ++i) {
                query[i] = static_cast<unsigned int>(suf[i]);
            }

//            // print query and target
//            for (int k = 0; k < queryLength; ++k) { 
//                printf("%d;%d  ", query[k], (unsigned int)(char)(suf[k]));
//            }
//            printf("\n");
//            for (int k = 0; k < targetLength; ++k) { 
//                printf("%d ", target[k]);
//            }

            edlibCalcEditDistance(query, queryLength, target, targetLength, 4,
                    1, EDLIB_MODE_HW, true, true, &score,
                    &endLocations, &startLocations, &numLocations,
                    &alignment, &alignmentLength);

            if (alignment) {
                std::cout << "--- NEW ALIGNMENT: " << std::endl;
                std::cout << "alignment length: " << alignmentLength << std::endl;
                std::cout << "query length: " << queryLength << std::endl;
                printAlignment(query, queryLength, target, length(target),
                        alignment, alignmentLength,
                        *(endLocations), EDLIB_MODE_HW, idxToLetter);
                printf("Best score: %d\n", score);
            }
        }
        bucketBegin = bucketEnd;
    }
    std::cout << cnt << std::endl; 
    std::cout << "elapsed time: "<< sysTime() - t << " sec\n";
}

int main(int argc, char *argv[])
{
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

