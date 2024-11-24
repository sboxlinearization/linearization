#include <bits/stdc++.h>
#include <assert.h>
#include <omp.h>

using namespace std;

// #define FULL_STAT 0

// typedef vector<pair<uint64_t, uint64_t>> STAT;
// one node
// sbatch -N 1 --ntasks-per-node 1 --ntasks-per-socket 1 -c 128 -t 1:00:00 --export=SBOX=aes,LB=18 ./batch.sh
// multiple nodes
// sbatch --tasks 8 -N 8 --ntasks-per-node 1 --ntasks-per-socket 1 -c 128 -t 48:00:00 --export=SBOX=apn8c1,LB=16 ./batch.sh
// sbatch --tasks 8 -N 8 --ntasks-per-node 1 --ntasks-per-socket 1 -c 128 -t 48:00:00 --export=SBOX=five8,LB=17 ./batch.sh
// sbatch --tasks 2 -N 2 --ntasks-per-node 1 --ntasks-per-socket 1 -c 128 -t 48:00:00 --export=SBOX=aes,LB=24 ./batch.sh
// sbatch --tasks 4 -N 4 --ntasks-per-node 1 --ntasks-per-socket 1 -c 128 -t 48:00:00 --export=SBOX=five9,LB=32 ./batch.sh
static vector<uint64_t> DERIVATIVES;

#ifdef SBOX_CUSTOM
    static int BITS;
    static bool IS_QUADRATIC = false;
    #ifdef SBOX_CUSTOM_LARGE
        typedef uint16_t word;
        static word S[65536];
    #else
        typedef uint8_t word;
        static word S[256];
    #endif
#endif

#ifdef SBOX_AES
    // AES S-box = inverse8
    constexpr int BITS = 8;
    constexpr bool IS_QUADRATIC = false;
    typedef uint8_t word;
    static const word S[] = {99, 124, 119, 123, 242, 107, 111, 197, 48, 1, 103, 43, 254, 215, 171, 118, 202, 130, 201, 125, 250, 89, 71, 240, 173, 212, 162, 175, 156, 164, 114, 192, 183, 253, 147, 38, 54, 63, 247, 204, 52, 165, 229, 241, 113, 216, 49, 21, 4, 199, 35, 195, 24, 150, 5, 154, 7, 18, 128, 226, 235, 39, 178, 117, 9, 131, 44, 26, 27, 110, 90, 160, 82, 59, 214, 179, 41, 227, 47, 132, 83, 209, 0, 237, 32, 252, 177, 91, 106, 203, 190, 57, 74, 76, 88, 207, 208, 239, 170, 251, 67, 77, 51, 133, 69, 249, 2, 127, 80, 60, 159, 168, 81, 163, 64, 143, 146, 157, 56, 245, 188, 182, 218, 33, 16, 255, 243, 210, 205, 12, 19, 236, 95, 151, 68, 23, 196, 167, 126, 61, 100, 93, 25, 115, 96, 129, 79, 220, 34, 42, 144, 136, 70, 238, 184, 20, 222, 94, 11, 219, 224, 50, 58, 10, 73, 6, 36, 92, 194, 211, 172, 98, 145, 149, 228, 121, 231, 200, 55, 109, 141, 213, 78, 169, 108, 86, 244, 234, 101, 122, 174, 8, 186, 120, 37, 46, 28, 166, 180, 198, 232, 221, 116, 31, 75, 189, 139, 138, 112, 62, 181, 102, 72, 3, 246, 14, 97, 53, 87, 185, 134, 193, 29, 158, 225, 248, 152, 17, 105, 217, 142, 148, 155, 30, 135, 233, 206, 85, 40, 223, 140, 161, 137, 13, 191, 230, 66, 104, 65, 153, 45, 15, 176, 84, 187, 22};
#endif

#ifdef SBOX_RANDOM8
    // Random S-box
    constexpr int BITS = 8;
    constexpr bool IS_QUADRATIC = false;
    typedef uint8_t word;
    static const word S[] = {198, 239, 134, 40, 240, 173, 176, 170, 94, 56, 106, 67, 251, 116, 239, 19, 84, 224, 47, 87, 241, 137, 79, 26, 2, 164, 162, 127, 111, 220, 72, 112, 241, 219, 107, 45, 230, 119, 103, 11, 170, 149, 165, 234, 98, 58, 246, 166, 158, 194, 191, 88, 224, 165, 150, 7, 22, 147, 112, 182, 99, 29, 9, 164, 7, 83, 245, 165, 129, 227, 95, 183, 140, 163, 34, 226, 50, 250, 29, 158, 21, 150, 15, 213, 105, 104, 154, 96, 248, 177, 209, 121, 19, 235, 175, 120, 14, 170, 25, 182, 229, 79, 23, 47, 52, 115, 54, 111, 22, 83, 44, 168, 103, 192, 222, 95, 134, 112, 74, 80, 173, 62, 49, 70, 27, 164, 202, 76, 231, 165, 57, 49, 222, 248, 74, 0, 75, 221, 209, 104, 55, 164, 182, 204, 222, 14, 0, 15, 208, 211, 121, 243, 234, 231, 98, 64, 182, 101, 120, 210, 214, 53, 69, 156, 190, 251, 252, 113, 69, 36, 212, 91, 91, 249, 83, 160, 37, 126, 4, 71, 157, 44, 1, 239, 80, 116, 110, 223, 63, 170, 103, 226, 16, 212, 0, 120, 242, 59, 182, 138, 13, 255, 1, 99, 150, 66, 150, 119, 195, 104, 4, 194, 193, 7, 174, 76, 7, 163, 1, 209, 171, 136, 145, 104, 184, 227, 244, 67, 41, 157, 171, 189, 101, 101, 141, 40, 246, 193, 247, 32, 248, 93, 103, 92, 34, 135, 230, 49, 221, 167, 147, 237, 12, 252, 16, 196};
#endif


#ifdef SBOX_INVERSE8
    // sage: monomial_function(8, 2**8-2)
    constexpr int BITS = 8;
    constexpr bool IS_QUADRATIC = false;
    typedef uint8_t word;
    static const word S[] = {0, 1, 142, 244, 71, 167, 122, 186, 173, 157, 221, 152, 61, 170, 93, 150, 216, 114, 192, 88, 224, 62, 76, 102, 144, 222, 85, 128, 160, 131, 75, 42, 108, 237, 57, 81, 96, 86, 44, 138, 112, 208, 31, 74, 38, 139, 51, 110, 72, 137, 111, 46, 164, 195, 64, 94, 80, 34, 207, 169, 171, 12, 21, 225, 54, 95, 248, 213, 146, 78, 166, 4, 48, 136, 43, 30, 22, 103, 69, 147, 56, 35, 104, 140, 129, 26, 37, 97, 19, 193, 203, 99, 151, 14, 55, 65, 36, 87, 202, 91, 185, 196, 23, 77, 82, 141, 239, 179, 32, 236, 47, 50, 40, 209, 17, 217, 233, 251, 218, 121, 219, 119, 6, 187, 132, 205, 254, 252, 27, 84, 161, 29, 124, 204, 228, 176, 73, 49, 39, 45, 83, 105, 2, 245, 24, 223, 68, 79, 155, 188, 15, 92, 11, 220, 189, 148, 172, 9, 199, 162, 28, 130, 159, 198, 52, 194, 70, 5, 206, 59, 13, 60, 156, 8, 190, 183, 135, 229, 238, 107, 235, 242, 191, 175, 197, 100, 7, 123, 149, 154, 174, 182, 18, 89, 165, 53, 101, 184, 163, 158, 210, 247, 98, 90, 133, 125, 168, 58, 41, 113, 200, 246, 249, 67, 215, 214, 16, 115, 118, 120, 153, 10, 25, 145, 20, 63, 230, 240, 134, 177, 226, 241, 250, 116, 243, 180, 109, 33, 178, 106, 227, 231, 181, 234, 3, 143, 211, 201, 66, 212, 232, 117, 127, 255, 126, 253};
#endif

#ifdef SBOX_CUBE8
    // sage: monomial_function(8, 2**8-2)
    constexpr int BITS = 8;
    constexpr bool IS_QUADRATIC = true;
    typedef uint8_t word;
    static const word S[] = {0, 1, 8, 15, 64, 85, 120, 107, 58, 115, 146, 221, 231, 186, 127, 36, 205, 193, 191, 181, 228, 252, 166, 184, 107, 47, 185, 251, 223, 143, 61, 107, 38, 115, 70, 21, 145, 208, 193, 134, 115, 110, 179, 168, 89, 80, 169, 166, 127, 39, 101, 59, 161, 237, 139, 193, 182, 166, 12, 26, 245, 241, 127, 125, 45, 161, 191, 53, 10, 146, 168, 54, 252, 56, 206, 12, 70, 150, 68, 146, 191, 62, 87, 208, 241, 100, 41, 186, 242, 59, 186, 117, 33, 252, 89, 130, 223, 7, 37, 251, 15, 195, 197, 15, 97, 241, 59, 173, 44, 168, 70, 196, 217, 12, 89, 138, 96, 161, 208, 23, 251, 102, 219, 64, 223, 86, 207, 64, 117, 231, 97, 245, 145, 23, 181, 53, 80, 138, 228, 56, 41, 231, 173, 101, 179, 44, 221, 68, 62, 181, 96, 237, 10, 221, 196, 21, 26, 217, 228, 33, 145, 87, 237, 45, 130, 80, 206, 26, 219, 85, 7, 143, 85, 207, 185, 37, 195, 8, 197, 8, 185, 102, 143, 86, 21, 150, 179, 54, 242, 101, 100, 245, 182, 169, 56, 33, 53, 62, 139, 134, 120, 47, 86, 7, 102, 37, 120, 61, 47, 61, 219, 207, 197, 195, 1, 1, 125, 39, 41, 117, 10, 68, 110, 38, 134, 205, 96, 45, 242, 173, 36, 125, 39, 36, 97, 100, 206, 217, 184, 169, 139, 205, 23, 87, 150, 196, 58, 110, 182, 184, 138, 130, 54, 44, 58, 38};
#endif
#ifdef SBOX_FIVE8
    // sage: monomial_function(8, 2**8-2)
    constexpr int BITS = 8;
    constexpr bool IS_QUADRATIC = true;
    typedef uint8_t word;
    static const word S[] = {0, 1, 32, 51, 116, 108, 46, 36, 38, 226, 1, 215, 169, 116, 244, 59, 180, 233, 17, 94, 32, 100, 255, 169, 132, 28, 38, 172, 235, 106, 51, 160, 3, 150, 108, 235, 26, 150, 15, 145, 116, 36, 28, 94, 150, 223, 132, 223, 77, 132, 167, 124, 180, 100, 36, 230, 44, 32, 193, 223, 46, 59, 185, 190, 96, 174, 55, 235, 1, 214, 44, 233, 103, 108, 55, 46, 253, 239, 215, 215, 38, 180, 244, 116, 167, 44, 15, 150, 55, 96, 226, 167, 77, 3, 226, 190, 85, 15, 77, 5, 89, 26, 59, 106, 3, 156, 28, 145, 244, 114, 145, 5, 233, 239, 116, 96, 5, 26, 226, 239, 169, 106, 51, 226, 190, 100, 94, 150, 156, 5, 100, 239, 174, 46, 44, 190, 32, 124, 223, 145, 233, 172, 108, 59, 124, 185, 1, 214, 174, 114, 169, 103, 214, 214, 172, 190, 255, 230, 255, 244, 180, 185, 3, 28, 235, 255, 38, 32, 89, 145, 233, 51, 253, 44, 55, 244, 174, 255, 156, 223, 17, 89, 89, 3, 85, 193, 96, 230, 17, 156, 94, 193, 114, 36, 253, 185, 85, 26, 160, 253, 239, 124, 103, 230, 51, 185, 193, 89, 96, 106, 106, 114, 167, 180, 215, 214, 235, 36, 230, 59, 215, 1, 160, 100, 108, 174, 172, 124, 38, 253, 156, 85, 160, 167, 103, 114, 17, 15, 172, 160, 132, 26, 193, 77, 46, 169, 17, 132, 94, 5, 28, 85, 15, 77, 55, 103};
#endif

#ifdef SBOX_APN8_C1
    // APN Candidate 1 for VecLin=15
    constexpr int BITS = 8;
    constexpr bool IS_QUADRATIC = true;
    typedef uint8_t word;
    static const word S[] = {0, 1, 0, 6, 0, 12, 2, 9, 0, 24, 50, 45, 4, 17, 52, 38, 0, 48, 87, 96, 100, 89, 49, 11, 8, 33, 109, 67, 104, 76, 15, 44, 0, 96, 107, 12, 173, 192, 196, 174, 200, 177, 145, 239, 97, 21, 58, 73, 16, 65, 44, 122, 217, 133, 231, 188, 208, 152, 222, 145, 29, 88, 17, 83, 0, 192, 76, 139, 213, 24, 155, 81, 91, 130, 37, 251, 138, 94, 246, 37, 146, 99, 137, 127, 35, 223, 58, 193, 193, 41, 232, 7, 116, 145, 95, 189, 32, 129, 7, 161, 88, 244, 125, 214, 179, 11, 166, 25, 207, 122, 216, 106, 162, 50, 210, 69, 190, 35, 204, 86, 57, 176, 123, 245, 33, 165, 97, 226, 0, 130, 128, 5, 152, 23, 26, 146, 171, 48, 25, 133, 55, 161, 135, 22, 181, 6, 98, 214, 73, 247, 156, 37, 22, 188, 243, 94, 238, 73, 9, 169, 38, 197, 205, 41, 19, 253, 250, 19, 69, 191, 156, 97, 116, 131, 175, 95, 131, 81, 63, 234, 210, 13, 108, 180, 232, 35, 102, 170, 189, 123, 49, 240, 64, 3, 140, 200, 13, 67, 195, 138, 176, 234, 78, 19, 249, 174, 5, 85, 103, 21, 252, 137, 78, 49, 215, 175, 159, 244, 54, 90, 178, 212, 25, 120, 70, 100, 225, 196, 166, 137, 3, 43, 126, 69, 235, 215, 154, 172, 13, 60, 113, 98, 129, 149, 245, 235, 7, 30, 65, 75, 131, 142, 193, 198, 1, 1};
#endif


#ifdef SBOX_INVERSE7
    // sage: monomial_function(7, 2**7-2)
    constexpr int BITS = 7;
    constexpr bool IS_QUADRATIC = false;
    typedef uint8_t word;
    static const word S[] = {0, 1, 65, 126, 97, 42, 63, 54, 113, 18, 21, 74, 94, 79, 27, 103, 121, 92, 9, 112, 75, 10, 37, 68, 47, 35, 102, 14, 76, 71, 114, 110, 125, 91, 46, 25, 69, 22, 56, 61, 100, 119, 5, 96, 83, 123, 34, 24, 86, 116, 80, 52, 51, 81, 7, 62, 38, 60, 98, 73, 57, 39, 55, 6, 127, 2, 108, 84, 23, 36, 77, 29, 99, 59, 11, 20, 28, 70, 95, 13, 50, 53, 122, 44, 67, 109, 48, 117, 104, 107, 124, 33, 17, 120, 12, 78, 43, 4, 58, 72, 40, 118, 26, 15, 88, 106, 105, 89, 66, 85, 31, 115, 19, 8, 30, 111, 49, 87, 101, 41, 93, 16, 82, 45, 90, 32, 3, 64};
#endif

#ifdef SBOX_CUBE7
    // sage: monomial_function(7, 2**7-2)
    constexpr int BITS = 7;
    constexpr bool IS_QUADRATIC = true;
    typedef uint8_t word;
    static const word S[] = {0, 1, 8, 15, 64, 85, 120, 107, 12, 69, 39, 104, 73, 20, 82, 9, 96, 119, 36, 53, 62, 61, 74, 79, 68, 27, 35, 122, 31, 84, 72, 5, 10, 51, 49, 14, 38, 11, 45, 6, 117, 4, 109, 26, 92, 57, 116, 23, 44, 3, 91, 114, 30, 37, 89, 100, 123, 28, 47, 78, 76, 63, 40, 93, 80, 113, 29, 58, 13, 56, 112, 67, 54, 95, 88, 55, 110, 19, 48, 75, 33, 22, 32, 17, 98, 65, 83, 118, 111, 16, 77, 52, 41, 66, 59, 86, 102, 127, 24, 7, 87, 90, 25, 18, 115, 34, 46, 121, 71, 2, 42, 105, 81, 94, 99, 106, 126, 101, 124, 97, 108, 43, 125, 60, 70, 21, 103, 50};
#endif

#ifdef SBOX_FIVE9
    // sage: monomial_function(8, 2**8-2)
    constexpr int BITS = 9;
    constexpr bool IS_QUADRATIC = true;
    typedef uint16_t word;
    static const word S[] = {0,1,32,51,34,295,83,324,98,227,466,321,53,432,468,67,38,243,142,73,397,92,372,183,147,198,427,492,333,28,36,359,226,229,159,138,328,75,356,117,56,191,469,320,487,100,91,458,249,42,12,205,218,269,126,443,244,167,401,464,162,501,406,211,174,158,78,108,377,77,456,238,84,228,292,390,502,66,215,369,307,471,347,429,365,141,340,166,286,378,230,144,309,85,156,494,479,489,354,326,384,178,364,76,413,299,176,20,439,5,203,363,127,157,74,186,425,79,461,57,234,136,335,319,329,47,189,457,362,200,345,233,388,290,486,338,71,357,484,212,220,250,302,282,213,419,110,266,434,448,344,312,303,217,4,480,61,207,327,423,323,231,301,155,37,133,26,168,214,498,296,30,453,485,106,88,193,433,39,325,46,90,153,255,387,115,245,279,281,493,62,216,45,446,222,351,310,417,404,273,408,139,251,506,246,225,452,449,41,366,82,263,187,248,145,192,331,140,160,373,428,367,22,199,407,2,313,190,260,405,507,376,154,399,420,163,124,109,275,272,430,239,392,219,180,241,195,148,116,437,450,17,283,478,252,43,54,414,460,118,101,201,462,368,280,48,370,72,318,274,261,315,164,472,470,184,382,262,93,55,349,161,447,81,242,10,65,171,125,467,474,102,422,268,80,232,491,197,476,224,69,111,35,27,210,424,509,149,128,254,510,402,403,105,300,196,436,330,346,438,308,173,14,389,146,15,505,374,130,411,40,291,337,332,426,421,29,336,431,240,50,123,465,394,380,177,350,129,294,495,341,398,236,371,395,6,194,89,500,381,482,253,277,24,441,418,287,278,504,179,23,334,95,16,481,444,33,490,94,391,499,316,477,256,386,393,107,114,285,18,165,440,483,360,410,259,265,134,289,188,137,86,488,293,415,68,175,358,63,96,206,131,348,7,508,181,258,271,182,169,21,284,496,235,475,342,511,352,185,304,204,343,52,237,264,451,170,375,455,8,58,99,150,221,209,396,44,355,361,339,64,104,3,317,379,87,400,298,297,385,143,305,103,459,473,311,120,132,314,208,202,306,503,409,454,442,353,11,257,121,122,70,270,288,152,416,445,151,59,135,223,113,172,276,25,435,247,31,267,497,412,112,49,463,97,9,13,119,383,19,322,60};
#endif

#ifdef SBOX_CUSTOM
    static int NUM;
#else
    constexpr int NUM = (1 << BITS);
#endif
// static word SCALAR[NUM][NUM];

static int LB;

struct STAT {
    uint64_t max_depth;
    vector<uint64_t> nodes;
    vector<uint64_t> leaves;
    uint64_t prefix_skips;

    #ifdef FULL_STAT
    vector<vector<uint64_t>> nodes_by_size;
    vector<vector<uint64_t>> leaves_by_size;
    #endif

    STAT(uint64_t _max_depth) : max_depth(_max_depth), prefix_skips(0) {
        nodes.resize(_max_depth + 1);
        leaves.resize(_max_depth + 1);
        #ifdef FULL_STAT
        nodes_by_size.assign(_max_depth + 1, vector<uint64_t>(NUM+1));
        leaves_by_size.assign(_max_depth + 1, vector<uint64_t>(NUM+1));
        #endif
    }
    void add(const STAT &cur_stat) {
        for(int depth = 0; depth <= max_depth; depth++) {
            nodes[depth] += cur_stat.nodes[depth];
            leaves[depth] += cur_stat.leaves[depth];
            #ifdef FULL_STAT
            for(int x = 0; x <= NUM; x++) {
                nodes_by_size[depth][x] += cur_stat.nodes_by_size[depth][x];
                leaves_by_size[depth][x] += cur_stat.leaves_by_size[depth][x];
            }
            #endif
        }
        prefix_skips += cur_stat.prefix_skips;
    }
    void report() const {
        uint64_t total_nodes = 0;
        uint64_t total_leaves = 0;
        for(int depth = 0; depth <= max_depth; depth++) {
            total_nodes += nodes[depth];
            total_leaves += leaves[depth];
            printf("depth %d:  %12lu = 2^%6.2lf nodes,  %12lu = 2^%6.2lf leaves,  %12lu = 2^%6.2f non-leaves (%.2f %%)\n",
                depth,
                nodes[depth], log(nodes[depth])/log(2),
                leaves[depth], log(leaves[depth])/log(2),
                (nodes[depth]-leaves[depth]), log(nodes[depth]-leaves[depth])/log(2),
                (nodes[depth]-leaves[depth]) * 100.0 / nodes[depth]
            );
            #ifdef FULL_STAT
            printf("  ");
            uint64_t total = 0;
            for(int x = 0; x <= NUM/4; x++) {
                auto num = nodes_by_size[depth][x];
                total += num;
                if (num == 0) {
                    printf("- ");
                }
                else {
                    printf("%5.2lf ", log(num)/log(2));
                }
                if (x % 16 == 15) printf("|\n  ");
            }
            printf("\n");
            printf("FULL STAT total: %lu = %.2lf\n", total, log(total)/log(2));
            #endif
        }
        printf("  total:  %12lu = 2^%6.2lf nodes,  %12lu = 2^%6.2lf leaves\n",
            total_nodes, log(total_nodes)/log(2),
            total_leaves, log(total_leaves)/log(2)
        );
        printf("  prefix   skips: %lu = 2^%6.2lf\n", prefix_skips, log(prefix_skips)/log(2));
        printf("  prefix checked: %lu = 2^%6.2lf\n", (1 << (2*BITS+2)) - prefix_skips, log((1 << (2*BITS+2)) - prefix_skips)/log(2));
    }
};



static void recurse(int depth, const vector<word> &X, STAT &stat) {
    stat.nodes[depth]++;
    #ifdef FULL_STAT
    stat.nodes_by_size[depth][X.size()]++;
    #endif
    if (X.size() < LB) {
        stat.leaves[depth]++;
        #ifdef FULL_STAT
        stat.leaves_by_size[depth][X.size()]++;
        #endif
        // N_LEAVES++;
        // if (depth > MAX_DEPTH) {
        //     printf("reached depth %d: %lu\n", depth, X.size());
        //     // for (auto v: X)
        //     //     printf("%d ", v);
        //     // printf("\n");
        //     MAX_DEPTH = depth;
        // }
        return;
    }
    if (depth == BITS) {
        #pragma omp critical
        {
            printf("SOLUTION depth %d: %lu : ", depth, X.size());
            for (auto v: X)
                printf("%d ", v);
            printf("\n");
            fflush(stdout);
        }
        return;
    }

    uint64_t ret = 0;
    for (uint64_t mask = 0; mask < NUM; mask++) {
        vector<word> Xnew[2];
        // split X into two sets according to the mask evaluation
        for (auto x: X) {
            // NOTE: LSB-based index at output (depth 0 = output LSB)
            uint8_t ybit = (S[x] >> depth) & 1;
            uint8_t xbit = __builtin_popcount(x & mask) & 1; // SCALAR[x][mask]
            Xnew[ybit ^ xbit].push_back(x);
        }
        recurse(depth+1, Xnew[0], stat);
        recurse(depth+1, Xnew[1], stat);
    }
    return;
}
static void recurse_start(
    word mask1, uint8_t const1, word mask2, uint8_t const2,
    STAT &stat
) {

    vector<word> X;
    for (uint64_t x = 0; x < NUM; x++) {
        uint8_t ybit1 = (S[x] >> 0) & 1;
        uint8_t xbit1 = __builtin_popcount(x & mask1) & 1; // SCALAR[x][mask]
        if ( (xbit1 ^ ybit1) != const1 )
            continue;
        uint8_t ybit2 = (S[x] >> 1) & 1;
        uint8_t xbit2 = __builtin_popcount(x & mask2) & 1; // SCALAR[x][mask]
        if ( (xbit2 ^ ybit2) != const2 )
            continue;

        X.push_back(x);
    }
    recurse(2, X, stat);
}

static void anf_in_place(vector<word> & anf, int n) {
    for (int i = 0; i < n; i++) {
        // xor up
        for (uint64_t j = 0; j < (1 << n); j++) {
            uint64_t j2 = j ^ (1 << i);
            if (j < j2) {
                anf[j2] ^= anf[j];
            }
        }
    }
}

static void verify_quadratic() {
    bool is_quadratic = true;
    vector<word> anf(S, S+NUM);
    anf_in_place(anf, BITS);
    for(uint64_t e = 0; e < NUM; e++) {
        if (__builtin_popcount(e) > 2) {
            is_quadratic &= !anf[e];
        }
    }
    #ifdef SBOX_CUSTOM
        IS_QUADRATIC = is_quadratic;
    #endif
    printf("S-box is quadratic? %d %d\n", is_quadratic, IS_QUADRATIC);
    fflush(stdout);
    assert(is_quadratic == IS_QUADRATIC);
    return;
}

static void compute_derivatives_quadratic() {
    set<uint64_t> der_set;
    for (uint64_t delta = 0; delta < NUM; delta++) {
        vector<word> anf(S, S+NUM);
        for (uint64_t x = 0; x < NUM; x++) {
            anf[x] = S[x] ^ S[x ^ delta];
        }
        anf_in_place(anf, BITS);

        // note: output bits indexed starting from LSB
        // input bits are same as masks (LSB=LSB, MSB=MSB)
        uint64_t mask1 = 0;
        uint64_t mask2 = 0;
        uint64_t const1 = anf[0] & 1;
        uint64_t const2 = (anf[0] & 2) >> 1;
        for (int i = 0; i < BITS; i++) {
            uint64_t coef = anf[1 << i];
            uint64_t bit1 = coef & 1;
            uint64_t bit2 = (coef & 2) >> 1;
            mask1 |= bit1 << i;
            mask2 |= bit2 << i;
        }
        uint64_t derivative = 0;
        // derivatives with constant part?
        // derivative |= const2;
        // derivative |= mask2 << 1;
        // derivative |= (const1 << (BITS + 1));
        // derivative |= (mask1 << (BITS + 2));

        // only the linear part matters now
        derivative |= mask2;
        derivative |= (mask1 << BITS);
        
        DERIVATIVES.push_back(derivative);
        der_set.insert(derivative);

        // testing
        for (uint64_t x = 0; x < NUM; x++) {
            uint64_t y = S[x] ^ S[x ^ delta];
            uint64_t val1 = (__builtin_popcount(x & mask1) & 1) ^ const1;
            uint64_t val2 = (__builtin_popcount(x & mask2) & 1) ^ const2;
            assert((y&1) == val1);
            assert( ((y&2) >> 1) == val2);
        }
    }
    printf("Computed %lu derivatives (%lu unique)\n", DERIVATIVES.size(), der_set.size());
    fflush(stdout);

    // keep unique ones
    DERIVATIVES.clear();
    for(auto der: der_set) {
        DERIVATIVES.push_back(der);
    }
    uint64_t n_der = DERIVATIVES.size();
    assert( (n_der & (n_der - 1)) == 0 ); // power of 2 - derivatives form a linear space
}

int main(int argc, char * argv[]) {
    #ifdef SBOX_CUSTOM
        if (argc != 5) {
            printf("Usage: %s <n> <lower bound> <batch-index> <batch-number>\n", argv[0]);
            return -1;
        }
        printf("Input S-box (space separated integers):\n");
        uint64_t n = atoi(argv[1]);
        printf("n = %ld\n", n);
        BITS = n;
        NUM = 1 << n;
        #ifdef SBOX_CUSTOM_LARGE
            assert(1 < n && n <= 16);
        #else
            assert(1 < n && n <= 8);
        #endif
        for(uint64_t x = 0; x < NUM; x++) {
            int y;
            assert(1 == scanf("%d", &y));
            S[x] = y;
            printf("%d ", y);
        }
        printf("\n");
    #else
        if (argc != 4) {
            printf("Usage: %s <lower bound> <batch-index> <batch-number>\n", argv[0]);
            return -1;
        }
        assert(sizeof(S)/sizeof(S[0]) == NUM);
    #endif

    verify_quadratic();
    if (IS_QUADRATIC) {
        compute_derivatives_quadratic();
    }

    #ifndef SBOX_CUSTOM
        LB = atoi(argv[1]);
        int batch_index = atoi(argv[2]);
        int batch_number = atoi(argv[3]);
    #else
        LB = atoi(argv[2]);
        int batch_index = atoi(argv[3]);
        int batch_number = atoi(argv[4]);
    #endif
    
    printf("RUNNING WITH LB %d, BATCH %d/%d\n", LB, batch_index, batch_number);
    fflush(stdout);
   
    assert(1 <= batch_number);
    assert(0 <= batch_index && batch_index < batch_number);

    // unused
    // for (uint64_t i = 0; i < NUM; i++) {
    // for (uint64_t j = 0; j < NUM; j++) {
    //     SCALAR[i][j] = __builtin_popcount(i & j) & 1;
    // }}

    STAT full_stat(BITS+1);

    // omp_set_dynamic(0);
    // omp_set_num_threads(128);

    // #pragma omp parallel for reduction (+:ret)
    #pragma omp parallel for
    for(uint64_t mask1 = 0; mask1 < NUM; mask1++) {
    #pragma omp parallel for
    for(uint64_t mask2 = 0; mask2 < NUM; mask2++) {
        uint64_t total_index = mask1 * NUM + mask2;
        if (total_index % batch_number != batch_index) {
            continue;
        }
        if (mask2 < batch_number) {
            #pragma omp critical
            {
                printf("checking depth 0 masks (%02lx %02lx)\n", mask1, mask2);
                fflush(stdout);
            }
        }

        STAT cur_stat(BITS+1);
        for(uint64_t const1 = 0; const1 < 2; const1++) {
        for(uint64_t const2 = 0; const2 < 2; const2++) {
            uint64_t full_prefix = 0;
            // derivatives with constants?
            // full_prefix |= const2;
            // full_prefix |= mask2 << 1;
            // full_prefix |= (const1 << (BITS + 1));
            // full_prefix |= (mask1 << (BITS + 2));

            // for now, only linear derivatives
            full_prefix |= mask2;
            full_prefix |= (mask1 << BITS);
            bool skip = false;
            for (auto der: DERIVATIVES) {
                if ((der ^ full_prefix) < full_prefix) {
                    // printf("skip %08lx %08x: %08x < this\n", full_prefix, der, full_prefix ^ der);
                    skip = true;
                    cur_stat.prefix_skips++;
                    break;
                }
            }
            if (!skip) {
                recurse_start(mask1, const1, mask2, const2, cur_stat);
            }
        }}

        #pragma omp critical
        full_stat.add(cur_stat);
    }}
    printf("\n");

    full_stat.report();

    fflush(stdout);
    return 0;
}