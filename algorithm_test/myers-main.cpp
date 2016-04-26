#include <iostream>
#include <sstream>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/index.h>    
#include <seqan/find.h>    
#include <seqan/seeds.h>    
#include <seqan/file.h>   // Required to print strings in tests.
#include <cmath>

#define QGRAM 19
//#define DEBUG

using namespace seqan;

template<typename VECTOR>
void printV(VECTOR const &v) {
    std::cout << std::endl;
    for (int i = 0; i < v.size(); ++i) {
        std::cout << v[i] << " ";
    }
    std::cout << std::endl;
}
template<typename SEQ, typename SEED>
void printSeed(SEQ const &seqH, SEQ const &seqV, SEED const &seed) {
//    std::cout << beginPositionH(seed) << " "
//              << beginPositionV(seed) << " "
//              << endPositionH(seed) << " "
//              << endPositionV(seed);

    std::cout << seed << std::endl
        <<  infix(seqH, beginPositionH(seed),
                endPositionH(seed)) << "\n"
        <<  infix(seqV, beginPositionV(seed),
                endPositionV(seed)) << "\n"
        << std::endl;
}

template<
    typename TSeed,
    typename THaystackFiber,
    typename TQuery>
void extendSimpleSeed(
    TSeed                 &seed,
    THaystackFiber  const &fiber,
    TQuery          const &query,
    int             const &k) {
        int posFiber;
        int posQuery;
        int mismatch;

        // keep list of where mismatches occur, left and right from seed
        std::vector<int> mismatchOffsetLeft;
        std::vector<int> mismatchOffsetRight;

        // go left
        mismatch = 0;
        posFiber = beginPositionH(seed);
        posQuery = beginPositionV(seed);
        while (posFiber >= 0 && posQuery >= 0 && mismatch <= k) {
            if (fiber[posFiber] != query[posQuery]) {
                mismatchOffsetLeft.push_back(beginPositionH(seed) - posFiber);
                mismatch++;
            }
            posFiber--;
            posQuery--;
        }

        // go right
        mismatch = 0;
        posFiber = endPositionH(seed);
        posQuery = endPositionV(seed);
        int endFiber = length(fiber);
        int endQuery = length(query);
        while (posFiber < endFiber && posQuery < endQuery && mismatch <= k + 1) {
            if (fiber[posFiber] != query[posQuery]) {
                mismatchOffsetRight.push_back(posFiber - endPositionH(seed));
                mismatch++;
            }
            posFiber++;
            posQuery++;
        }

        // find maximal extension
        int leftMismatches  = mismatchOffsetLeft.size();
        int rightMismatches = mismatchOffsetRight.size();

        // III if #mismatches less than k, whole segments match
        if (leftMismatches + rightMismatches <= k) {
//            assert(false);
            // update left position
            int bH = beginPositionH(seed);
            setBeginPositionH(seed, std::max(0, (int)(beginPositionH(seed) - beginPositionV(seed))));
            setBeginPositionV(seed, std::max(0, (int)(beginPositionV(seed) - bH)));

            // update right positions
            int eH = endPositionH(seed);
            setEndPositionH(seed, std::min((int)(length(fiber)), (int)(endPositionH(seed) + length(query) - endPositionV(seed))));
            setEndPositionV(seed, std::min((int)(length(query)), (int)(endPositionV(seed) + length(fiber) - eH)));
            printSeed(fiber, query, seed);
            return;
        }

        int limitL          = -1;
        int limitR          = -1;

        // add corresponding end points if required
        if (leftMismatches <= k) {
            mismatchOffsetLeft.push_back(std::min((int)(beginPositionH(seed)),(int)(beginPositionV(seed))));
            limitL = leftMismatches;
            ++leftMismatches;
        }
        if (rightMismatches <= k) {
            mismatchOffsetRight.push_back(std::min((int)(length(query) - endPositionV(seed) - 1), 
                    (int)(length(fiber) - endPositionH(seed) - 1)));
            limitR = rightMismatches;
            ++rightMismatches;
        }

//        std::cout << std::endl << "MismatchOffsetLeft:";
//        printV(mismatchOffsetLeft);
//        std::cout << std::endl << "MismatchOffsetRight:";
//        printV(mismatchOffsetRight);
//        std::cout << "\n-----------------------\n";

        int diagonal;           // size of matching segment between fiber and query
        int maxDiagonal  = 0;
        int seedDiagonal = endPositionH(seed) - beginPositionH(seed);
        int l, r;
        int lAdd, rAdd;
        TSeed tempSeed   = seed;

        // init left and right pointers
        l = leftMismatches - 1;
        r = (leftMismatches <= k) ? k - leftMismatches + 1 : 0;
            
//        std::cout << "limitL, limitR: " << limitL << ", " << limitR << std::endl;

        for (; l >= 0 && r < rightMismatches; --l) {
            assert(abs(l - r) <= k);
            diagonal = mismatchOffsetRight[r] + mismatchOffsetLeft[l] + seedDiagonal;
            if (l != limitL && r != limitR) {
                --diagonal;
            }

//            std::cout << "l #" << l << " : " << mismatchOffsetLeft[l]
//                      << std::endl << "r #" << r << " : " << mismatchOffsetRight[r]
//                      << std::endl << "diagonal: " << diagonal << std::endl;

            if (diagonal > maxDiagonal) {
                maxDiagonal = diagonal;
                lAdd = (l == limitL) ? 0 : 1;
                rAdd = (r == limitR) ? -1 : 0; // -1 only because the right end of seeds in seqan have an offset of 1
                // update left position
                setBeginPositionH(tempSeed, beginPositionH(seed) - mismatchOffsetLeft[l] + lAdd);
                setBeginPositionV(tempSeed, beginPositionV(seed) - mismatchOffsetLeft[l] + lAdd);
                // update right positions
                setEndPositionH(tempSeed, endPositionH(seed) + mismatchOffsetRight[r] - rAdd);
                setEndPositionV(tempSeed, endPositionV(seed) + mismatchOffsetRight[r] - rAdd);
//                std::cout << "MaxSeed is now: " << tempSeed << std::endl;
            }
            r++;
        }
        seed = tempSeed;
        printSeed(fiber, query, seed);
}

void seed() {
    //                 0  3 5  8 0   4               30
    CharString seqH = "GGGXGXG X AGAGG A YZ ABC A GAGG";
    CharString seqV =   "GAGGG G AGAGG A YZ ABC A GAGG";
    //                   0       8   2               28

    Seed<Simple> seed(10, 8, 15, 13);
//    std::cout << "Original seed:\n"
//        << "seedH: " << infix(seqH, beginPositionH(seed),
//                endPositionH(seed)) << "\n"
//        << "seedV: " << infix(seqV, beginPositionV(seed),
//                endPositionV(seed)) << "\n"
//        << seed;
    extendSimpleSeed(seed, seqH, seqV, 2);
}

int main(int argc, char *argv[])
{
    seed();
    return 0;
}

/*
// unit tests
 ################################################################################################
 ## k == 2 ##
    //                 0       8 0   4               30
    CharString seqH = "GGGAGGG X AGAGG X YZ ABC X GAGG";
    CharString seqV =      "GG G AGAGG A YZ ABC X GAGG";
    ///                     0    5   9 
 >>>
 Seed<Simple, TConfig>(5, 0, 30, 25, lower diag = 5, upper diag = 5)
 seedH: GG X AGAGG X YZ ABC X GAG
 seedV: GG G AGAGG A YZ ABC X GAG

 ################################################################################################
 ## k == 2 ##
    //                 0  3 5  8 0   4               30
    CharString seqH = "GGGXGXG X AGAGG X YZ ABC X GXGG";
    CharString seqV =   "GAGGG G AGAGG A YZ ABC A GAGG";
    //                   0       8   2               28
diagonal: 19
Seed<Simple, TConfig>(6, 4, 25, 23, lower diag = 2, upper diag = 2)
seedH: G X AGAGG X YZ ABC 
seedV: G G AGAGG A YZ ABC 


 ################################################################################################
 ## k == 3 ##
 //                 0  3 5  8 0   4               30
 CharString seqH = "GGGXGXG X AGAGG X YZ ABX X GAGG";
 CharString seqV =   "GAGGG G AGAGG A YZ ABC A GAGG";
 //                   0       8   2               28

diagonal: 22
MaxSeed is now: Seed<Simple, TConfig>(9, 7, 31, 29, lower diag = 2, upper diag = 2)
Seed<Simple, TConfig>(9, 7, 31, 29, lower diag = 2, upper diag = 2)
seedH:  AGAGG X YZ ABX X GAGG
seedV:  AGAGG A YZ ABC A GAGG

 ################################################################################################
 ## k == 3 ##
    //                 0  3 5  8 0   4               30
    CharString seqH = "GGGXGXG X AGAGG A YZ ABC A GAGG";
    CharString seqV =   "GAGGG G AGAGG A YZ ABC A GAGG";
    //                   0       8   2               28

 Seed<Simple, TConfig>(2, 0, 31, 29, lower diag = 2, upper diag = 2)
 GXGXG X AGAGG A YZ ABC A GAGG
 GAGGG G AGAGG A YZ ABC A GAGG
*/
