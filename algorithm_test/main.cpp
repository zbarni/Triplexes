#include <iostream>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/index.h>    
#include <seqan/find.h>    
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
void runSimpleIndexSearch(TMotifSet tfoSet, TttsSet ttsSet) {
    SEQAN_PROTIMESTART(timeT);
    double t = sysTime();

    typedef typename Value<TMotifSet>::Type         TTfo;
    typedef Shape<Dna5, UngappedShape<QGRAM> >      TShape;
    typedef Index<TttsSet, IndexQGram<UngappedShape<QGRAM>, OpenAddressing> > TQGramIndex;

    long long cnt = 0;
    TShape shape;
    TQGramIndex index(ttsSet);  

    for (int i = 0; i < length(tfoSet); ++i) {
        TTfo tfo = tfoSet[i];
        for (int j = 0; j <= length(tfo) - QGRAM; ++j) {
            TTfo qGram = infix(tfo, j, j + QGRAM);
            Finder< TQGramIndex > finder(index);

            while(find(finder, qGram)) {
                cnt ++;
            }
        }
    }
    printResults("simpleIndexSearch", sysTime() - t, cnt);
}


////////////////////////////////////////////
// new one: search for 
template<typename TMotifSet, typename TttsSet>
void runIndexSearch(TMotifSet tfoSet, TttsSet ttsSet) {
}

////////////////////////////////////////////
// original: tfo is preprocessed (pattern)
//           tts (haystack) is not processed in any way
template<typename THaystack, typename TPattern>
void runOnlineSearch(THaystack haystack, TPattern pattern) {
    SEQAN_PROTIMESTART(timeT);
    double t = sysTime();

    typedef Index<TPattern, IndexQGram<UngappedShape<QGRAM>, OpenAddressing> > TQGramIndex;
    typedef Shape<Dna5, UngappedShape<QGRAM> >          TShape;
    typedef typename Value<TShape>::Type                TShapeValue;
    typedef typename Fibre<TQGramIndex, QGramSA>::Type  TSA;
    typedef typename Iterator<TSA, Standard>::Type      TSAIter;

    long long cnt = 0;
    TShape shape;

    TQGramIndex index(pattern);  

    indexRequire(index, QGramSADir());
    for (int i = 0; i < length(haystack); ++i) {
        hashInit(shape, begin(haystack[i]));
        for (unsigned j = 0; j < length(haystack[i]) - length(shape) + 1; ++j) {
            TShapeValue hash = hashNext(shape, begin(haystack[i]) + j);

            TSAIter saBegin = begin(indexSA(index), Standard());
            TSAIter occ     = saBegin + indexDir(index)[getBucket(index.bucketMap, hash)];
            TSAIter occEnd  = saBegin + indexDir(index)[getBucket(index.bucketMap, hash) + 1];

            for(; occ != occEnd; ++occ) {
                cnt++;
            }
        }
    }
    printResults("onlineSearch", sysTime() - t, cnt);
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

    while (fTts >> s) {
        rTTS = static_cast<Ttts>(s);
        appendValue(ttsSet, rTTS);
    }
    
    fTfo.close();
    fTts.close();

//    std::cout << "haystack: ttsSet\n";
//    runOnlineSearch(ttsSet, tfoSet);
    std::cout << "haystack: tfoSet\n";
    runOnlineSearch(tfoSet, ttsSet);

    std::cout << "index, ttsSet into Finder\n";
    runSimpleIndexSearch(tfoSet, ttsSet);
    return 0;
}

