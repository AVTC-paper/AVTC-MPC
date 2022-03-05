#ifndef AVTC_H_
#define AVTC_H_

#include "stdbool.h"
#include <stdio.h>
#include <gmp.h>
#include "stdlib.h"
#include <string.h>
#include <sys/random.h>
#include "omp.h"
#include <sys/time.h>
#include <openssl/rand.h>
#include <openssl/sha.h>


#define ARITH_MAXDEPTH 7040
#define SECRET_LENGTH 128
#define ITERATIONS 26
#define RANDOM_TAPE_SEED_SIZE 256


#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))


size_t fieldSize;
mpz_t field;
mpz_t blindingFactorUpperbound;
mpz_t plaintextShareUpperBound;
mpz_t N;
mpz_t Nsquare;
gmp_randstate_t globalRand;


typedef struct {
    mpz_t incMsgs0[ARITH_MAXDEPTH];
    mpz_t incMsgs1[ARITH_MAXDEPTH];
} OLEOutputs;

typedef struct {
    mpz_t x0;
    mpz_t x1;
    mpz_t r0;
    mpz_t r1;
    mpz_t y0;
    mpz_t y1;
} OLEInput;

typedef struct {
    OLEInput msgs[ARITH_MAXDEPTH];
} OLEInputs;

typedef struct {
    char b_0;
    char b_1;
    mpz_t r_0;
    mpz_t r_1;
} ConversionShare;

typedef struct {
    ConversionShare shares[SECRET_LENGTH];
} ConversionShares;

typedef struct {
    int reconstructionIndex;
    mpz_t shares[1665];
} ReconstructionShares;

typedef struct {
    char blindedBit;
    char blindingBit0;
    char blindingBit1;
    mpz_t blindedShare;
    mpz_t blindingShare0;
    mpz_t blindingShare1;
} LSBExtractionShares;

typedef struct {
    int lsbExtractionDepth;
    LSBExtractionShares lsbExtractionShares[ARITH_MAXDEPTH];
    int distributionDepth;
    mpz_t sharesFromParty0[4609];
    mpz_t sharesFromParty1[4609];
} AllLSBExtractionShares;

typedef struct {
    int *depth;
    OLEOutputs *fromOLE[3];
    OLEInputs *toOLE[3];
    int *conversionIndex;
    ConversionShares *conversions[3];
    ReconstructionShares *reconstructions[3];
    AllLSBExtractionShares *extractionShares[3];
} Views;

typedef struct {
    mpz_t r1;
    mpz_t r2;
} r_ij;

typedef struct {
    int invocationNum;
    int party;
    const unsigned char *seed;
    gmp_randstate_t r;
} RandomTape;


typedef struct {
    mpz_t shares[3];
} SS;

typedef struct {
    mpz_t b, t;
    char lsb;
    SS b_shares;
    SS t_shares;
} PrecomputedShares;

typedef struct {
    mpz_t b;
    mpz_t t;
} BTShares;


typedef struct {
    OLEInputs *oleInputs;
    OLEOutputs *oleOutputs;
    const unsigned char *randomnessSeed;
    BTShares precomputedShares[SECRET_LENGTH + 1];
    ConversionShares *conversionShares;
    ReconstructionShares *reconstructions;
    AllLSBExtractionShares *lsbExtractionShares;
} CommittedView;

typedef struct {
    CommittedView committedViews[3];
} CommittedViews;

typedef struct {
    int party;
    BTShares precomputedShares[SECRET_LENGTH + 1];
    AllLSBExtractionShares *lsbExtractionShares;
    ConversionShares *conversionShares;
    OLEOutputs *oleOutputs;
    ReconstructionShares *reconstructions;
    RandomTape *tape;
} OpenedView;

typedef struct {
    OLEInput msgs[ARITH_MAXDEPTH];
} ComputedToOLE;

typedef struct {
    int *depth;
    OpenedView views[2];
} OpenedViewTuple;

typedef struct {
    SS next;
    SS lsb;
} LayerOutput;

typedef struct {
    PrecomputedShares *elements;
    int length;
} ModReducePrecomp;

typedef struct {
    char shares[SECRET_LENGTH];
} BooleanShares;

typedef struct {
    BooleanShares booleanShares[ITERATIONS];
    CommittedViews simulations[ITERATIONS];
} Algorithm1Output;


#define ySize 736


typedef struct {
    unsigned char x[64];
    uint32_t y[ySize];
} V;

typedef struct {
    uint32_t yp[3][8];
    unsigned char h[3][32];
} A;


typedef struct {
    A as[ITERATIONS];
    unsigned char rs[ITERATIONS][3][4];
    unsigned char keys[ITERATIONS][3][16];
    V localViews[ITERATIONS][3];
} sha256TotalViews;

typedef struct {
    unsigned char ke[16];
    unsigned char ke1[16];
    V ve;
    V ve1;
    unsigned char re[4];
    unsigned char re1[4];
} Z;

sha256TotalViews runSHA256Proof(unsigned char shares[ITERATIONS][3][SECRET_LENGTH / 8]);

void simulateVerifierViewSelection(int *es);

void randBlindingFactor(RandomTape *tape, mpz_t out, mpz_t blindingFactorUpperbound);

int nextRandomTapeEntry(mpz_t out, RandomTape *tape, mpz_t upperBound);

void randFieldElement(mpz_t out, gmp_randstate_t rand);

int partyIndex(int i, int j);

int fromIndex(int i, int j);

void subModField(mpz_t out, mpz_t x, mpz_t y);

void mpc_Mul(mpz_t x[3], mpz_t y[3], mpz_t z[3], RandomTape *tapes[3], Views views);

SS secretShareNoTape(mpz_t w);

RandomTape *newRandomTape();

Views createViews();

void reconstructArithStandalone(SS ss, mpz_t x);

RandomTape *tapeFromSeed(const unsigned char *rawSeed);

void testMul();

void testXor();

void testBooleanizeLSB();

SS mpc_xor(SS x, SS y, RandomTape *tapes[3], Views views);

int extractLSBByReconstruction(Views views, SS v, mpz_t vPlaintext, mpz_t lsb);

SS zero();

BooleanShares
verifyAlgorithm1(SS z, const char *S, OpenedView *v1, OpenedView *v2, ComputedToOLE *ole1, ComputedToOLE *ole2);

OpenedViewTuple *openView(Algorithm1Output *out, int i, int nonRevealedParty);

Z *packSHA256A(int *es, sha256TotalViews sha256Views);

int verifySHA256(A a, int e, Z z, unsigned char inputs[2][SECRET_LENGTH / 8]);

char booleanizeLSB(SS in, RandomTape *tapes[3], Views views);

#endif /* AVTC_H_ */