#include <iostream>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/index.h>    

using namespace seqan;

int main(int argc, char *argv[])
{
    std::ifstream fTfo (argv[1]);
    std::ifstream fTts (argv[2]);

    typedef Dna5String       Ttts;
    typedef Dna5String       Ttfo;
    typedef StringSet<Ttfo>  TMotifSet;
    typedef StringSet<Ttts>  TttsSet;
    typedef Index<TttsSet, IndexQGram<UngappedShape<9>, OpenAddressing> > TQGramIndex;
    typedef Pattern<TQGramIndex> TPattern;
        
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

    ////////////////////////////////////////////
    // create index
    TQGramIndex index(ttsSet);

    resize(indexShape(index), 9);
    // create pattern   
//    TPattern pattern(index, UngappedShape<9> );
    return 0;
    
//    typedef Index<DnaString, IndexQGram<UngappedShape<4>, OpenAddressing > > TIndex;
//    TIndex indexQ("CATGATTACATA");
//    Finder<TIndex> finder(indexQ);
//    stringToShape(indexShape(indexQ), "1111");
//    
//    //    hash(indexShape(index), "ATCA");
//    //    for (unsigned i = 0; i < length(getOccurrences(index, indexShape(index))); ++i)
//    //        std::cout << getOccurrences(index, indexShape(index))[i] << std::endl;
//
//    CharString needle = "CATGA";
//
//    Pattern<CharString> pattern(needle);
//    while (find(finder, pattern))
//        std::cout << '[' << beginPosition(finder) << ',' << endPosition(finder) << ")\t" << infix(finder) << std::endl;
//    return 0;
}

