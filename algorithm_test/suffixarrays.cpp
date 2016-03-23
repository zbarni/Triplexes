#include <iostream>
#include <sstream>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/index.h>    
#include <seqan/find.h>    
#include <seqan/file.h>   // Required to print strings in tests.
#include <seqan/seq_io.h>

#define QGRAM 3

using namespace seqan;

void printResults(std::string fn, double t, long long cnt) {

    std::cout   << std::endl
        << "Function:\t\t" << fn << std::endl
        << "Total time:\t\t" << t << std::endl
        << "Found Counter:\t" << cnt << std::endl;
}

////////////////////////////////////////////
// haystack: tts is preprocessed into a finder 
template<typename TMotifSet, typename TttsSet>
void createSuffixData(TMotifSet tfoSet, TttsSet ttsSet) {

    double t = sysTime();
    typedef typename Value<TMotifSet>::Type         TTfo;
    typedef typename Value<TttsSet>::Type           TTts;

    TTfo symbol = static_cast<TTfo>(std::string("Y"));

    for (int i = 0; i < length(tfoSet); ++i) {
        TTfo tfo = tfoSet[i];
        for (int j = 0; j < length(ttsSet); ++j) {
            TTts tts = ttsSet[j];
            append(tfo,symbol);
            append(tfo,tts);
            Index<TTfo, IndexEsa<> > *esaIndex = new Index<TTfo, IndexEsa<> >(tfo);
            
            indexRequire(*esaIndex, FibreLcp());

//            for (unsigned i = 0; i < length(indexLcp(*esaIndex)); ++i)
//            {
//                std::cout << lcpAt(i, *esaIndex) << " ";
//            }
//            std::cout << std::endl << tfo << std::endl;
            delete esaIndex;
        }
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
//        std::cout << rTFO << std::endl;
        appendValue(tfoSet, rTFO);
    }

    std::string concatenated;
    while (fTts >> s) {
        concatenated.append(s);
//        rTTS = static_cast<Ttts>(s);
//        appendValue(ttsSet, rTTS);
    }
    appendValue(ttsSet, static_cast<Ttts>(concatenated));
    
    fTfo.close();
    fTts.close();
    createSuffixData(tfoSet, ttsSet);

    return 0;
}

