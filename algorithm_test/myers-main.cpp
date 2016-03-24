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
    typedef typename Value<TShape>::Type                TShapeValue;
    typedef typename Fibre<TQGramIndex, QGramSA>::Type      TSA;
    typedef typename Fibre<TQGramIndex, QGramSADir>::Type   TSADir;
    typedef typename Value<TSADir>::Type                    TSADirValue;
    typedef typename Value<TSA>::Type                       TSAValue;
    typedef typename Iterator<TSA, Standard>::Type          TSAIter;
    typedef typename Iterator<TSADir, Standard>::Type       TSADirIter;

    typedef typename Fibre<TQGramIndex, QGramCounts>::Type       TCounts;
    typedef typename Fibre<TQGramIndex, QGramCountsDir>::Type    TCountsDir;
    typedef typename Value<TCountsDir>::Type                TDirValue;
    typedef typename Iterator<TCounts, Standard>::Type      TIterCounts;
    typedef typename Iterator<TCountsDir, Standard>::Type   TIterCountsDir;
    
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

    TSADirIter itSADir              = begin(indexDir(index), Standard());
    TSADirIter itSADirEnd           = end(indexDir(index), Standard());
    TSAIter saBegin                 = begin(indexSA(index), Standard());

    // for each bucket count common q-grams for each sequence pair
//    TDirValue bucketBegin = *itCountsDir;
//    for (++itCountsDir; itCountsDir != itCountsDirEnd; ++itCountsDir)
//    {
//        TDirValue bucketEnd = *itCountsDir;
//
//        // q-gram must occur in at least 2 different sequences
//        if (bucketBegin != bucketEnd)
//        {
//            TIterCounts itA = itCountsBegin + bucketBegin;
//            TIterCounts itEnd = itCountsBegin + bucketEnd;
//            for (; itA != itEnd; ++itA)
//            {
////                ;
//                std::cout << *itA << std::endl;
//                TSAIter occ = saBegin + indexDir(index)[bucketBegin];
//                TSAIter occEnd = saBegin + indexDir(index)[bucketBegin + 1];
//
//                for(; occ != occEnd; ++occ)
//                {
//                    Pair<unsigned> ndlPos;
//                    posLocalize(ndlPos, *occ, stringSetLimits(index));
//                    std::cout << ndlPos << std::endl;
//                    std::cout << "bucket: " << bucketBegin << std::endl;
//                }
//            }
//        }
//        bucketBegin = bucketEnd;
//    }

    TSADirValue bucketBegin = *itSADir;
    for (++itSADir; itSADir != itSADirEnd; ++itSADir)
    {
        TSADirValue bucketEnd = *itSADir;
        if (bucketBegin != bucketEnd)
        {
            TSAIter itA = saBegin + bucketBegin;
            TSAIter itEnd = saBegin + bucketEnd;
            std::cout << *itA << std::endl;
            for (; itA != itEnd; ++itA)
            {
//                std::cout << *itA << std::endl;
//                printf("shit\n");
            }
        }
        bucketBegin = bucketEnd;
    }

//    for (unsigned i = 0; i < length(indexSA(index)); ++i)
//    {
//        unsigned textPos = (saAt(i, index).i2 == 0) ? length(index) - 1 : saAt(i, index).i2 - 1;
//        std::cout << textAt(textPos, index) << "\t" << suffix(indexText(index), textPos) << std::endl;
//        std::cout << infix(indexText(index), i, i + 3) << std::endl;
//    }
    
exit(0);














//    char idxToLetter[4] = {'A', 'C', 'G', 'T'};
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

    // for each tfo in set
    for (int iter = 0; iter < length(tfoSet); ++iter) {
        TTfo tfo = tfoSet[iter];
        unsigned char query[length(tfo)];
        unsigned int queryLength = length(tfo);

        // transform query 
        for (int i = 0; i < length(tfo); ++i) {
            query[i] = static_cast<unsigned int>(tfo[i]);
        }

//        // print query and target
//        for (int k = 0; k < queryLength; ++k) { 
//            printf("%d;%d  ", query[k], (unsigned int)(char)(tfo[k]));
//        }
//        printf("\n");
//        for (int k = 0; k < targetLength; ++k) { 
//            printf("%d ", target[k]);
//        }

        edlibCalcEditDistance(query, queryLength, target, targetLength, 4,
                3, EDLIB_MODE_HW, true, true, &score,
                &endLocations, &startLocations, &numLocations,
                &alignment, &alignmentLength);

//        if (alignment) {
//        std::cout << "alignment length: " << alignmentLength << std::endl;
//        std::cout << "query length: " << queryLength << std::endl;
//            printAlignment(query, length(tfo), target, length(target),
//                    alignment, alignmentLength,
//                    *(endLocations), EDLIB_MODE_HW, idxToLetter);
//            printf("Best score: %d\n", score);
//        }
    }
    
    std::cout << "elapsed time: "<< sysTime() - t << " sec\n";
}

template<typename TMotifSet, typename TttsSet>
void readFasta(std::string& fTfo, std::string& fTts, TMotifSet& tfoSet, TttsSet& ttsSet) {
//    StringSet<CharString> ids, tfoIds;
//
//    SeqFileIn seqFileIn(fTts.c_str());
//    SeqFileIn tfoFileIn(fTfo.c_str());
//
//    // Reads all remaining records.
//    readRecords(ids, ttsSet, seqFileIn);
//    readRecords(tfoIds, tfoSet, tfoFileIn);
}

int main(int argc, char *argv[])
{
    std::ifstream fTfo (argv[1]);
    std::ifstream fTts (argv[2]);
//    std::string fTfo = argv[1];
//    std::string fTts = argv[2];

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
    //    readFasta(fTfo, fTts, tfoSet, ttsSet);
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

