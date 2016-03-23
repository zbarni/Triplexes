#include <iostream>
#include <sstream>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/index.h>    
#include <seqan/find.h>    
#include <seqan/file.h>   // Required to print strings in tests.
#include <seqan/seq_io.h>
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
    typedef typename Value<TMotifSet>::Type         TTfo;
    typedef typename Value<TttsSet>::Type           TTts;

    char idxToLetter[4] = {'A', 'C', 'G', 'T'};
    int score1, numLocations1;
    int *endLocations1, *startLocations1;
    unsigned char *alignment; int alignmentLength;

    TTts tts = ttsSet[0];
    unsigned char target[length(tts)];
    unsigned int targetLength = length(tts);
    for (int i = 0; i < length(tts); ++i) {
        target[i] = static_cast<unsigned int>(tts[i]);
    }

//    unsigned int charToIdx[128];
//    charToIdx['A'] = 0;
//    charToIdx['C'] = 1;
//    charToIdx['T'] = 2;
//    charToIdx['G'] = 3;

    // for each tfo in set
    for (int iter = 0; iter < length(tfoSet); ++iter) {
        TTfo tfo = tfoSet[iter];
        unsigned char query[length(tfo)];
        unsigned int queryLength = length(tfo);

        // transform query and target
        for (int i = 0; i < length(tfo); ++i) {
//            query[i] = charToIdx[(unsigned int)(char)(tfo[i])];
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
                2, EDLIB_MODE_HW, true, true, &score1,
                &endLocations1, &startLocations1, &numLocations1,
                &alignment, &alignmentLength);

//        std::cout << "alignment length: " << alignmentLength << std::endl;
//        std::cout << "query length: " << queryLength << std::endl;
//        if (alignment) {
//            printAlignment(query, length(tfo), target, length(target),
//                    alignment, alignmentLength,
//                    *(endLocations1), EDLIB_MODE_HW, idxToLetter);
//        }
    }
    
    std::cout << "elapsed time: "<< sysTime() - t << " sec\n";
}

template<typename TMotifSet, typename TttsSet>
void readFasta(std::string& fTfo, std::string& fTts, TMotifSet& tfoSet, TttsSet& ttsSet) {
    StringSet<CharString> ids, tfoIds;

    SeqFileIn seqFileIn(fTts.c_str());
    SeqFileIn tfoFileIn(fTfo.c_str());

    // Reads all remaining records.
    readRecords(ids, ttsSet, seqFileIn);
    readRecords(tfoIds, tfoSet, tfoFileIn);
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

