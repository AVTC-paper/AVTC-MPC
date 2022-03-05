#include "avtc.h"

// WARNING:
// This code is for educational purposes only.
// Do not use this in any commercial or production environment.

void simulateVerifierViewSelection(int* es) {
    for (int i=0; i < ITERATIONS; i++) {
        unsigned char randByte[1];
        if (getrandom(randByte, 1, GRND_RANDOM) == -1) {
            fprintf(stderr, "not enough randomness");
            exit(2);
        }
        es[i] = randByte[0] % 3;
    }
}

int findSkippedParty(OpenedView *v1, OpenedView *v2) {
    int skippedParty = -1;
    for (int i=0; i < 3; i++) {
        if (i == v1->party || i == v2->party) {
            continue;
        }
        return i;
    }

    return -1;
}

void verifyAddScalar(mpz_t in[3], mpz_t out[3], mpz_t n, int skippedIndex) {
    for (int i=0; i <= 2; i++) {
        if (i == skippedIndex) {
            continue;
        }
        mpz_set(out[i], in[i]);
    }

    for (int i=0; i <= 2; i++) {
        if (i == skippedIndex) {
            continue;
        }
        if (i == 0) {
            mpz_add(out[i], out[i], n);
            mpz_mod(out[i], out[i], field);
        }
    }
}

void assertEq(mpz_t a, mpz_t b, const char *errMsg, int depth, int currParty, int jIndex) {
    if (mpz_cmp(a, b) == 0) {
        return;
    }
    gmp_fprintf(stderr, errMsg, depth, currParty, jIndex, a, b);
    fflush(stderr);
    int n = 9/0;
    exit(2);
}

void verifyAdd(mpz_t x[3], mpz_t y[3], mpz_t z[3], int skippedIndex) {
    for (int i=0; i <= 2; i++) {
        if (i == skippedIndex) {
            continue;
        }
        mpz_add(z[i], x[i], y[i]);
        mpz_mod(z[i], z[i], field);
    }
}

void verifySub(mpz_t x[3], mpz_t y[3], mpz_t z[3], int skippedIndex) {
    mpz_t negY[3];
    for (int i=0; i <= 2; i++) {
        if (i == skippedIndex) {
            continue;
        }
        mpz_init(negY[i]);
        mpz_sub(negY[i], field, y[i]);
    }

    return verifyAdd(x, negY, z, skippedIndex);
}

void verifyMulScalar(mpz_t in[3], mpz_t out[3], mpz_t n, int skippedIndex) {
    for (int i = 0; i < 3; i++) {
        if (i == skippedIndex) {
            continue;
        }
        mpz_mul(out[i], in[i], n);
        mpz_mod(out[i], out[i], field);
    }
}

// x, y and z are all 3 sized extractionShares, but one of them is un-initialized because
// we only have knowledge of two out of three of the extractionShares.
void verifyMul(mpz_t x[3], mpz_t y[3], mpz_t z[3], OpenedView *v1, OpenedView *v2, ComputedToOLE *ole1, ComputedToOLE *ole2, int *circuitDepth) {
    int depth = *circuitDepth;

    bool secondParty = false;
    for (int i = 0; i <= 2; i++) { // i is the current party index
        // Skip the third party, we are not verifying it
        int currentParty;
        if (i != v1->party && i != v2->party) {
            continue;
        }
        RandomTape *tape;
        ComputedToOLE *toOLE;
        OLEOutputs *sentFromOLE;

        if (secondParty) {
            tape = v2->tape;
            toOLE = ole2;
            sentFromOLE = v1->oleOutputs;
        } else {
            tape = v1->tape;
            toOLE = ole1;
            sentFromOLE = v2->oleOutputs;
        }
        nextRandomTapeEntry(toOLE->msgs[depth].r0, tape, field);
        nextRandomTapeEntry(toOLE->msgs[depth].r1, tape, field);

        for (int j = 0; j <= 1; j++) { // j is either 0 or 1 to indicate the two other parties
            mpz_t XiYj;
            mpz_init(XiYj);

            int jIndex = partyIndex(i, j);
            int from = fromIndex(i, jIndex);


            if (j == 0) {
                mpz_init_set(toOLE->msgs[depth].x0, x[i]);
            } else {
                mpz_init_set(toOLE->msgs[depth].x1, x[i]);
            }

            if (j == 0) {
                mpz_init_set(toOLE->msgs[depth].y0, y[i]);
            } else {
                mpz_init_set(toOLE->msgs[depth].y1, y[i]);
            }

            // Skip the third party, we are not simulating interaction with it
            if (jIndex != v1->party && jIndex != v2->party) {
                continue;
            }

            mpz_mul(XiYj, x[i], y[jIndex]);

            mpz_t Sij;
            mpz_init(Sij);
            if (j == 0) {
                mpz_add(Sij, XiYj, toOLE->msgs[depth].r0);
            } else {
                mpz_add(Sij, XiYj, toOLE->msgs[depth].r1);
            }
            mpz_mod(Sij, Sij, field);

            if (from == 0) {
                assertEq(sentFromOLE->incMsgs0[depth], Sij, "depth: %d, i: %d, jIndex: %d, S_{i,j} received from 0 does not match expected (%Zd, %Zd)\n", depth, i, jIndex);
            } else {
                assertEq(sentFromOLE->incMsgs1[depth], Sij, "depth: %d, i: %d, jIndex: %d, S_{i,j} received from 1 does not match expected (%Zd, %Zd)\n", depth, i, jIndex);
            }

        }

        secondParty = true;
    }

    secondParty = false;
    for (int i = 0; i <= 2; i++) {
        // Skip the third party, we are not verifying it
        int currentParty;
        if (i != v1->party && i != v2->party) {
            continue;
        }
        ComputedToOLE *toOLE;
        OpenedView *openedView;

        if (secondParty) {
            openedView = v2;
            toOLE = ole2;
        } else {
            openedView = v1;
            toOLE = ole1;
        }

        mpz_t XiYi;
        mpz_init(XiYi);
        mpz_mul(XiYi, x[i], y[i]);


        mpz_t Ci;
        mpz_init(Ci);


        mpz_add(Ci, XiYi, openedView->oleOutputs->incMsgs0[depth]);
        mpz_add(Ci, Ci, openedView->oleOutputs->incMsgs1[depth]);
        mpz_mod(Ci, Ci, field);

        subModField(Ci, Ci, toOLE->msgs[depth].r0);
        subModField(Ci, Ci, toOLE->msgs[depth].r1);

        mpz_set(z[i], Ci);

        secondParty = true;

    }

    (*circuitDepth)++;
}


void testMul() {

    mpz_t x;
    mpz_init_set_si(x, 20);

    mpz_t y;
    mpz_init_set_si(y, 30);

    SS xs = secretShareNoTape(x);
    SS ys = secretShareNoTape(y);

    SS zs;
    mpz_init(zs.shares[0]);
    mpz_init(zs.shares[1]);
    mpz_init(zs.shares[2]);

    RandomTape *tapes[3];
    for (int i = 0; i < 3; i++) {
        tapes[i] = newRandomTape();
    }

    Views views = createViews();

    mpc_Mul(xs.shares, ys.shares, zs.shares, tapes, views);

    mpz_t z;
    mpz_init(z);
    reconstructArithStandalone(zs, z);
    long res = mpz_get_si(z);
    if (600 != res) {
        printf("Expected 600 but got %d", res);
        exit(2);
    }

    *views.depth = 0;

    int skippedParty = 1;

    OpenedView openedViews[2];
    int index = 0;
    for (int i=0; i <= 2; i++) {
        // Reveal only parties 0 and 1
        if (i == skippedParty) {
            continue;
        }

        openedViews[index].party = i;
        openedViews[index].tape = tapeFromSeed(tapes[i]->seed);
        openedViews[index].oleOutputs = views.fromOLE[i];
        index++;
    }

    // Backup zs[0]
    mpz_t skippedZ;
    mpz_init_set(skippedZ, zs.shares[skippedParty]);

    mpz_clear(xs.shares[skippedParty]);
    mpz_clear(ys.shares[skippedParty]);

    mpz_clear(zs.shares[0]);
    mpz_clear(zs.shares[1]);
    mpz_clear(zs.shares[2]);
    mpz_init(zs.shares[0]);
    mpz_init(zs.shares[1]);
    mpz_init(zs.shares[2]);

    int depth = 0;
    ComputedToOLE toOLE[2];

    verifyMul(xs.shares, ys.shares, zs.shares, &openedViews[0], &openedViews[1], &toOLE[0], &toOLE[1], &depth);

    mpz_init_set(zs.shares[skippedParty], skippedZ);

    reconstructArithStandalone(zs, z);
    res = mpz_get_si(z);
    if (600 != res) {
        printf("Expected 600 but got %d", res);
        exit(2);
    }

    int parties[3][2] = {
            {1, 2},
            {0, 2},
            {0, 1}
    };

    for (int i=0; i < depth; i++) {
        for (int p = 0; p <= 1; p++) {
            int party = parties[skippedParty][p];
            OLEInputs *in = views.toOLE[party];
            if (mpz_cmp(in->msgs[i].x0, toOLE[p].msgs[i].x0) != 0) {
                gmp_printf("i: %d, party: %d, p: %d, p.x0: %Zd\n x0: %Zd\n", i, party, p, in->msgs[i].x0, toOLE[p].msgs[i].x0);
                exit(2);
            }
            if (mpz_cmp(in->msgs[i].x1, toOLE[p].msgs[i].x1) != 0) {
                gmp_printf("i: %d, p: %d, p.x1: %Zd\n x1: %Zd\n", i, p, in->msgs[i].x1, toOLE[p].msgs[i].x1);
                exit(2);
            }
        }
    }

    printf("PASS");
}

SS verifyNeg(SS x, int skippedIndex) {
    SS res;

    for (int i=0; i < 3; i++) {
        if (i == skippedIndex) {
            continue;
        }
        mpz_init(res.shares[i]);
        mpz_sub(res.shares[i], field, x.shares[i]);
    }
    return res;
}

void initSharesInRevealedParties(mpz_t x[3], int skippedIndex) {
    for (int i=0; i < 3; i++) {
        if (i == skippedIndex) {
            continue;
        }
        mpz_init(x[i]);
    }
}



SS verifyXOR(SS x, SS y, OpenedView *v1, OpenedView *v2, ComputedToOLE *ole1, ComputedToOLE *ole2, int *circuitDepth) {
    int skippedIndex = findSkippedParty(v1, v2);

    SS xplusy;
    initSharesInRevealedParties(xplusy.shares, skippedIndex);
    verifyAdd(x.shares, y.shares, xplusy.shares, skippedIndex);

    mpz_t xy[3];
    initSharesInRevealedParties(xy, skippedIndex);

    verifyMul(x.shares, y.shares, xy, v1, v2, ole1, ole2, circuitDepth);


    SS twoXY;
    initSharesInRevealedParties(twoXY.shares, skippedIndex);
    verifyAdd(xy, xy, twoXY.shares, skippedIndex); // 2xy

    SS twoXYNeg = verifyNeg(twoXY, skippedIndex);

    SS z;
    initSharesInRevealedParties(z.shares, skippedIndex);
    verifyAdd(twoXYNeg.shares, xplusy.shares, z.shares, skippedIndex);

    return z;
}

void testXor() {
    for (int x = 0; x <= 1; x++) {
        for (int y = 0; y <= 1; y++) {
            mpz_t n;
            mpz_init_set_si(n, x);

            mpz_t m;
            mpz_init_set_si(m, y);

            SS ns = secretShareNoTape(n);
            SS ms = secretShareNoTape(m);

            RandomTape *tapes[3];
            for (int i = 0; i < 3; i++) {
                tapes[i] = newRandomTape();
                tapes[i]->party = i;
            }

            Views views = createViews();

            SS result = mpc_xor(ns, ms, tapes, views);

            mpz_t out;
            mpz_init(out);
            reconstructArithStandalone(result, out);

            int actual = mpz_get_si(out);
            int expected = x ^y;

            if (actual != expected) {
                fprintf("expected %, got %d\n", expected, actual);
                fflush(stderr);
                exit(2);
            }

            OpenedView openedViews[2];
            int index = 0;
            for (int i=0; i <= 2; i++) {
                // Reveal only parties 0 and 1
                if (i == 2) {
                    continue;
                }

                openedViews[index].party = i;
                openedViews[index].tape = tapeFromSeed(tapes[i]->seed);
                openedViews[index].oleOutputs = views.fromOLE[i];
                index++;
            }

            // Backup zs[2]
            mpz_t z2;
            mpz_init_set(z2, result.shares[2]);

            mpz_clear(result.shares[0]);
            mpz_clear(result.shares[1]);
            mpz_clear(result.shares[2]);
            mpz_init(result.shares[0]);
            mpz_init(result.shares[1]);
            mpz_init(result.shares[2]);

            int depth = 0;
            ComputedToOLE toOLE0;
            ComputedToOLE toOLE2;

            mpz_clear(ns.shares[2]);
            mpz_clear(ms.shares[2]);

            result = verifyXOR(ns, ms, &openedViews[0], &openedViews[1], &toOLE0, &toOLE2, &depth);

            mpz_init_set(result.shares[2], z2);

            mpz_clear(out);
            mpz_init(out);
            reconstructArithStandalone(result, out);
            actual = mpz_get_si(out);
            expected = x ^y;

            if (actual != expected) {
                fprintf("expected %, got %d\n", expected, actual);
                fflush(stderr);
                exit(2);
            }
        }
    }
}

void assertBlindingFactorInRange(mpz_t x) {
    if (mpz_cmp(x, blindingFactorUpperbound) < 0 ) {
        return;
    }
    fprintf(stderr, "blinding factor too large");
    fflush(stderr);
    exit(2);
}

void assertPlaintextShareInRange(mpz_t x) {
    if (mpz_cmp(x, plaintextShareUpperBound) < 0 ) {
        return;
    }
    fprintf(stderr, "plaintext share too large (%d bits)", mpz_sizeinbase(x, 2));
    fflush(stderr);
    exit(2);
}

SS verifySecretShare(mpz_t s, OpenedView *v) {
    mpz_t x1;
    mpz_t x2;
    mpz_t x3;
    mpz_init(x1);
    mpz_init(x2);
    mpz_init(x3);

    randFieldElement(x1, v->tape->r);
    randFieldElement(x2, v->tape->r);
    // (x1+x2+x3)%field = s
    mpz_t sum;
    mpz_init(sum);
    mpz_add(sum, x1, x2);
    mpz_sub(sum, s, sum);     // sum = s - (x1+x2)
    while (mpz_sgn(sum) == -1) {
        mpz_add(sum, sum, field);
    }
    mpz_set(x3, sum);

    SS ss;
    mpz_init(ss.shares[0]);
    mpz_init(ss.shares[1]);
    mpz_init(ss.shares[2]);
    mpz_set(ss.shares[0], x1);
    mpz_set(ss.shares[1], x2);
    mpz_set(ss.shares[2], x3);

    return ss;
}

SS verifySecretShareBit(int x, OpenedView *v) {
    if (x > 1 || x < 0) {
        fprintf(stderr, "x is %d but should be 0 or 1\n", x);
        fflush(stderr);
        exit(2);
    }
    mpz_t n;
    mpz_init_set_si(n, x);
    return verifySecretShare(n, v);
}

OpenedView *findReconstructionShareView(OpenedView *v1, OpenedView *v2) {
    int skippedParty = findSkippedParty(v1, v2);
    int partyStorageIndex = (skippedParty+1) % 3;
    if (v1->party == partyStorageIndex) {
        return v1;
    }
    return v2;
}

int verifyReconstruct(SS ss, mpz_t out, OpenedView *v1, OpenedView *v2) {
    int skippedParty = findSkippedParty(v1, v2);
    OpenedView *v = findReconstructionShareView(v1, v2);
    int reconsIndex = v->reconstructions->reconstructionIndex;
    mpz_init(ss.shares[skippedParty]);
    mpz_set(ss.shares[skippedParty], v->reconstructions->shares[reconsIndex]);
    v->reconstructions->reconstructionIndex++;
    reconstructArithStandalone(ss, out);
    return reconsIndex;
}

void verifyZeroAssertion(SS ss, const char *caller, OpenedView *v1, OpenedView *v2) {
    mpz_t sum;
    mpz_init(sum);

    int reconsIndex = verifyReconstruct(ss, sum, v1, v2);
    if (mpz_sgn(sum) != 0) {
        fprintf(stderr, "reconstruction of '%s' has %d bits, index is: %d\n", caller, mpz_sizeinbase(sum, 2),
                reconsIndex);
        gmp_printf("Shares:\n  %Zd\n  %Zd\n  %Zd\n", ss.shares[0], ss.shares[1], ss.shares[2]);
        fflush(stderr);
        int divideByZero = 9/0;
        exit(2);
    }
}

SS verifyComputeLSB(SS v, OpenedView *v1, OpenedView *v2, ComputedToOLE *ole1, ComputedToOLE *ole2, int *circuitDepth) {
    int skippedParty = findSkippedParty(v1, v2);

    // Input validation:
    int depth = v1->lsbExtractionShares->lsbExtractionDepth;

    // Check that blinded values and blinding factors are within appropriate range
    assertPlaintextShareInRange(v1->lsbExtractionShares->lsbExtractionShares[depth].blindedShare);
    assertPlaintextShareInRange(v2->lsbExtractionShares->lsbExtractionShares[depth].blindedShare);
    assertBlindingFactorInRange(v1->lsbExtractionShares->lsbExtractionShares[depth].blindingShare0);
    assertBlindingFactorInRange(v1->lsbExtractionShares->lsbExtractionShares[depth].blindingShare1);
    assertBlindingFactorInRange(v2->lsbExtractionShares->lsbExtractionShares[depth].blindingShare0);
    assertBlindingFactorInRange(v2->lsbExtractionShares->lsbExtractionShares[depth].blindingShare1);


    // Check that blinded bits are indeed bits
    char bitSharesAndBitHidingFactors[2][3] = {
            {v1->lsbExtractionShares->lsbExtractionShares[depth].blindingBit1, v1->lsbExtractionShares->lsbExtractionShares[depth].blindingBit0, v1->lsbExtractionShares->lsbExtractionShares[depth].blindedBit},
            {v2->lsbExtractionShares->lsbExtractionShares[depth].blindingBit1, v2->lsbExtractionShares->lsbExtractionShares[depth].blindingBit0, v2->lsbExtractionShares->lsbExtractionShares[depth].blindedBit},
    };

    for (int i=0; i<2; i++) {
        for (int j=0; j < 3; j++) {
            if (bitSharesAndBitHidingFactors[i][j] > 1 || bitSharesAndBitHidingFactors[i][j] < 0) {
                fprintf(stderr, "pre-shared bit is not a bit\n");
                exit(2);
            }
        }
    }

    // Secret sharing
    SS blindedPlaintextShares[3];
    SS blindingFactorsShares1[3];
    SS blindingFactorsShares2[3];
    SS blindedBitsShares[3];
    SS blindingBitsShares1[3];
    SS blindingBitsShares2[3];

    // Parties secret share their local input

    bool secondParty = false;
    for (int i = 0; i < 3; i++) {
        if (i == skippedParty) {
            mpz_init(blindedPlaintextShares[skippedParty].shares[0]);
            mpz_init(blindedPlaintextShares[skippedParty].shares[1]);
            mpz_init(blindedPlaintextShares[skippedParty].shares[2]);
            mpz_init(blindingFactorsShares1[skippedParty].shares[0]);
            mpz_init(blindingFactorsShares1[skippedParty].shares[1]);
            mpz_init(blindingFactorsShares1[skippedParty].shares[2]);
            mpz_init(blindingFactorsShares2[skippedParty].shares[0]);
            mpz_init(blindingFactorsShares2[skippedParty].shares[1]);
            mpz_init(blindingFactorsShares2[skippedParty].shares[2]);
            mpz_init(blindedBitsShares[skippedParty].shares[0]);
            mpz_init(blindedBitsShares[skippedParty].shares[1]);
            mpz_init(blindedBitsShares[skippedParty].shares[2]);
            mpz_init(blindingBitsShares1[skippedParty].shares[0]);
            mpz_init(blindingBitsShares1[skippedParty].shares[1]);
            mpz_init(blindingBitsShares1[skippedParty].shares[2]);
            mpz_init(blindingBitsShares2[skippedParty].shares[0]);
            mpz_init(blindingBitsShares2[skippedParty].shares[1]);
            mpz_init(blindingBitsShares2[skippedParty].shares[2]);
            continue;
        }
        OpenedView *v;
        if (secondParty) {
            v = v2;
        } else {
            v = v1;
        }
        blindedPlaintextShares[i] = verifySecretShare(v->lsbExtractionShares->lsbExtractionShares[depth].blindedShare, v);
        blindingFactorsShares1[i] = verifySecretShare(v->lsbExtractionShares->lsbExtractionShares[depth].blindingShare0, v);
        blindingFactorsShares2[i] = verifySecretShare(v->lsbExtractionShares->lsbExtractionShares[depth].blindingShare1, v);
        blindedBitsShares[i] = verifySecretShareBit(v->lsbExtractionShares->lsbExtractionShares[depth].blindedBit, v);
        blindingBitsShares1[i] = verifySecretShareBit(v->lsbExtractionShares->lsbExtractionShares[depth].blindingBit0, v);
        blindingBitsShares2[i] = verifySecretShareBit(v->lsbExtractionShares->lsbExtractionShares[depth].blindingBit1, v);
        v->lsbExtractionShares->lsbExtractionDepth++;

        secondParty = true;
    }

    // Populate secret shares distributed from the non-revealed party
    bool otherPartyIs0[3][3] = {
            {false, true, false},
            {true, false, false},
            {true, false, false},
    };

    secondParty = false;
    for (int i = 0; i < 3; i++) {
        if (i == skippedParty) {
            continue;
        }
        OpenedView *v;
        if (secondParty) {
            v = v2;
        } else {
            v = v1;
        }
        if (otherPartyIs0[v->party][skippedParty]) {
            mpz_set(blindedPlaintextShares[skippedParty].shares[v->party],
                    v->lsbExtractionShares->sharesFromParty0[v->lsbExtractionShares->distributionDepth]);
        } else {
            mpz_set(blindedPlaintextShares[skippedParty].shares[v->party],
                    v->lsbExtractionShares->sharesFromParty1[v->lsbExtractionShares->distributionDepth]);
        }
        v->lsbExtractionShares->distributionDepth++;
        if (otherPartyIs0[v->party][skippedParty]) {
            mpz_set(blindingFactorsShares1[skippedParty].shares[v->party], v->lsbExtractionShares->sharesFromParty0[v->lsbExtractionShares->distributionDepth]);
        } else {
            mpz_set(blindingFactorsShares1[skippedParty].shares[v->party], v->lsbExtractionShares->sharesFromParty1[v->lsbExtractionShares->distributionDepth]);
        }
        v->lsbExtractionShares->distributionDepth++;
        if (otherPartyIs0[v->party][skippedParty]) {
            mpz_set(blindingFactorsShares2[skippedParty].shares[v->party], v->lsbExtractionShares->sharesFromParty0[v->lsbExtractionShares->distributionDepth]);
        } else {
            mpz_set(blindingFactorsShares2[skippedParty].shares[v->party], v->lsbExtractionShares->sharesFromParty1[v->lsbExtractionShares->distributionDepth]);
        }
        v->lsbExtractionShares->distributionDepth++;
        if (otherPartyIs0[v->party][skippedParty]) {
            mpz_set(blindedBitsShares[skippedParty].shares[v->party], v->lsbExtractionShares->sharesFromParty0[v->lsbExtractionShares->distributionDepth]);
        } else {
            mpz_set(blindedBitsShares[skippedParty].shares[v->party], v->lsbExtractionShares->sharesFromParty1[v->lsbExtractionShares->distributionDepth]);
        }
        v->lsbExtractionShares->distributionDepth++;
        if (otherPartyIs0[v->party][skippedParty]) {
            mpz_set(blindingBitsShares1[skippedParty].shares[v->party], v->lsbExtractionShares->sharesFromParty0[v->lsbExtractionShares->distributionDepth]);
        } else {
            mpz_set(blindingBitsShares1[skippedParty].shares[v->party], v->lsbExtractionShares->sharesFromParty1[v->lsbExtractionShares->distributionDepth]);
        }
        v->lsbExtractionShares->distributionDepth++;
        if (otherPartyIs0[v->party][skippedParty]) {
            mpz_set(blindingBitsShares2[skippedParty].shares[v->party], v->lsbExtractionShares->sharesFromParty0[v->lsbExtractionShares->distributionDepth]);
        } else {
            mpz_set(blindingBitsShares2[skippedParty].shares[v->party], v->lsbExtractionShares->sharesFromParty1[v->lsbExtractionShares->distributionDepth]);
        }
        v->lsbExtractionShares->distributionDepth++;
        secondParty = true;
    }


    SS v0 = zero();
    for (int i = 0; i < 3; i++) {
        v0 = verifyXOR(v0, blindedBitsShares[i], v1, v2, ole1, ole2, circuitDepth);
        v0 = verifyXOR(v0, blindingBitsShares1[i], v1, v2, ole1, ole2, circuitDepth);
        v0 = verifyXOR(v0, blindingBitsShares2[i], v1, v2, ole1, ole2, circuitDepth);
    }

    SS sumBlinded = zero();
    for (int i = 0; i < 3; i++) {
        verifyAdd(sumBlinded.shares, blindedPlaintextShares[i].shares, sumBlinded.shares, skippedParty);
    }

    SS sumBlindingFactors = zero();
    for (int i = 0; i < 3; i++) {
        verifyAdd(sumBlindingFactors.shares, blindingFactorsShares1[i].shares, sumBlindingFactors.shares, skippedParty);
        verifyAdd(sumBlindingFactors.shares, blindingFactorsShares2[i].shares, sumBlindingFactors.shares, skippedParty);
    }

    SS shouldBeZero = zero();
    verifySub(v.shares, sumBlinded.shares, shouldBeZero.shares, skippedParty);
    verifyAdd(shouldBeZero.shares, sumBlindingFactors.shares, shouldBeZero.shares, skippedParty);
    verifySub(shouldBeZero.shares, v0.shares, shouldBeZero.shares, skippedParty);

    verifyZeroAssertion(shouldBeZero, "computeLSB", v1, v2);
    return v0;
}

SS verifyXORPlaintext(mpz_t x[3], int y, OpenedView *v1, OpenedView *v2) {
    int skippedParty = findSkippedParty(v1, v2);
    SS zs = zero();

    // (x xor y) = x+y - 2*x*y
    // if y == 1:  (x xor y) = x+1 - 2*x
    // if y == 0:  (x xor y) = x
    if (y == 0) {
        for (int i = 0; i < 3; i++) {
            if (i == skippedParty) {
                continue;
            }
            mpz_set(zs.shares[i], x[i]);
        }
        return zs;
    }

    if (y != 1) {
        fprintf(stderr, "y is %d but expected 1\n", y);
        exit(2);
    }

    mpz_t one;
    mpz_init_set_si(one, 1);

    verifyAddScalar(x, zs.shares, one, skippedParty);

    SS twoX;
    mpz_init(twoX.shares[0]);
    mpz_init(twoX.shares[1]);
    mpz_init(twoX.shares[2]);
    verifyAdd(x, x, twoX.shares, skippedParty);

    SS minusTwoX = verifyNeg(twoX, skippedParty);
    verifyAdd(zs.shares, minusTwoX.shares, zs.shares, skippedParty);
    return zs;
}

char verifyBooleanizeLSB(SS in, OpenedView *v1, OpenedView *v2, ComputedToOLE *ole1, ComputedToOLE *ole2, int *circuitDepth, int *convIndex) {
    int skippedParty = findSkippedParty(v1, v2);

    int conversionIndex = *convIndex;
    (*convIndex)++;

    mpz_t ourShare[3];
    char ourBit[3];
    char bitsSentToParties[3][2];

    bool secondParty = false;
    for (int i = 0; i < 3; i++) {
        if (i == skippedParty) {
            continue;
        }
        OpenedView *v;
        if (secondParty) {
            v = v2;
        } else {
            v = v1;
        }
        // Create blinding factors
        mpz_t r_i;
        randBlindingFactor(v->tape, r_i, blindingFactorUpperbound);

        int lsb = mpz_tstbit(r_i, 0);

        // Secret share the blinding factors in both domains;
        SS rs = verifySecretShare(r_i, v);
        mpz_init_set(ourShare[i], rs.shares[i]);
        for (int j=0; j<3; j++) {
            if (skippedParty == j) {
                continue;
            }
        }

        mpz_t tmp;
        randBlindingFactor(v->tape, tmp, blindingFactorUpperbound);
        int b0 = mpz_tstbit(tmp, 0); // random bit for party 0
        int b1 = mpz_tstbit(tmp, 1); // random bit for party 1
        ourBit[i] = b0 ^ b1 ^ lsb; // We get the XOR of the rest of the bits
        bitsSentToParties[i][0] = b0;
        bitsSentToParties[i][1] = b1;

        secondParty = true;
    }

    char lsbPi[3];
    secondParty = false;
    for (int i=0; i < 3; i++) {
        if (i == skippedParty) {
            continue;
        }

        OpenedView *v;
        if (secondParty) {
            v = v2;
        } else {
            v = v1;
        }
        lsbPi[i] = ourBit[i] ^ v->conversionShares->shares[conversionIndex].b_0 ^ v->conversionShares->shares[conversionIndex].b_1;
        secondParty = true;
    }

    //printf("lsbP0: %d, lsbP2: %d\n", lsbPi[0], lsbPi[2]);

    // Sum all arithmetic blinding factors and create a sharing 'r'
    SS r;
    secondParty = false;
    for (int i = 0; i < 3; i++) {
        if (i == skippedParty) {
            continue;
        }

        OpenedView *v;
        if (secondParty) {
            v = v2;
        } else {
            v = v1;
        }
        mpz_init(r.shares[i]);
        mpz_add(r.shares[i], ourShare[i], v->conversionShares->shares[conversionIndex].r_0);
        mpz_mod(r.shares[i], r.shares[i], field);
        mpz_add(r.shares[i], r.shares[i], v->conversionShares->shares[conversionIndex].r_1);
        mpz_mod(r.shares[i], r.shares[i], field);

        secondParty = true;
    }

    SS x;
    mpz_init(x.shares[0]);
    mpz_init(x.shares[1]);
    mpz_init(x.shares[2]);

    // Blind arithmetic secret with blinding factor
    verifyAdd(r.shares, in.shares, x.shares, skippedParty);

    // Reconstruct and reveal blinded secret
    mpz_t xPlaintext;
    mpz_init(xPlaintext);

    verifyReconstruct(x, xPlaintext, v1, v2);

    int lsb = mpz_tstbit(xPlaintext, 0);
    // Without loss of generality, add LSB to party 0's share
    if (v1->party == 0) {
        lsbPi[v1->party] = lsbPi[v1->party] ^ lsb;
    }

    if (v2->party == 0) {
        lsbPi[v2->party] = lsbPi[v2->party] ^ lsb;
    }

    // Encode output by putting party i's sharing in bit i.
    char out = 0;
    for (int i=0; i < 3; i++) {
        if (i == skippedParty) {
            continue;
        }
        out = out | (lsbPi[i] << i);
    }
    return out;
}

void verifyLessThanHalfOrder(SS v, OpenedView *v1, OpenedView *v2, ComputedToOLE *ole1, ComputedToOLE *ole2, int *circuitDepth) {
    int skippedParty = findSkippedParty(v1, v2);
    SS w;
    mpz_init(w.shares[0]);
    mpz_init(w.shares[1]);
    mpz_init(w.shares[2]);
    verifyAdd(v.shares, v.shares, w.shares, skippedParty);
    SS shouldBeZero = verifyComputeLSB(w, v1, v2, ole1, ole2, circuitDepth);
    verifyZeroAssertion(shouldBeZero, "lessThanHalfOrder", v1, v2);
}

void verifyLowerThan(SS v, mpz_t u, OpenedView *v1, OpenedView *v2, ComputedToOLE *ole1, ComputedToOLE *ole2, int *circuitDepth) {
    int skippedParty = findSkippedParty(v1, v2);
    verifyLessThanHalfOrder(v, v1, v2, ole1, ole2, circuitDepth);
    SS w = verifyNeg(v, skippedParty);

    SS x = zero();
    verifyAddScalar(w.shares, x.shares, u, skippedParty);
    verifyLessThanHalfOrder(x, v1, v2, ole1, ole2, circuitDepth);
}

LayerOutput verifySquareModNExtractLSB(SS x_i, OpenedView *v1, OpenedView *v2, ComputedToOLE *ole1, ComputedToOLE *ole2, int *circuitDepth, int bitIndex) {
    int skippedParty = findSkippedParty(v1, v2);
    SS a;
    mpz_init(a.shares[0]);
    mpz_init(a.shares[1]);
    mpz_init(a.shares[2]);

    verifyMul(x_i.shares, x_i.shares, a.shares,v1, v2, ole1, ole2, circuitDepth); // a = x^2

    SS y;
    mpz_init(y.shares[0]);
    mpz_init(y.shares[1]);
    mpz_init(y.shares[2]);

    bool secondParty = false;

    SS b;
    mpz_init(b.shares[0]);
    mpz_init(b.shares[1]);
    mpz_init(b.shares[2]);

    SS t;
    mpz_init(t.shares[0]);
    mpz_init(t.shares[1]);
    mpz_init(t.shares[2]);

    for (int i=0; i<3; i++) {
        if (i == skippedParty) {
            continue;
        }
        OpenedView *v;
        if (secondParty) {
            v = v2;
        } else {
            v = v1;
        }
        mpz_set(b.shares[i], v->precomputedShares[bitIndex].b);
        mpz_set(t.shares[i], v->precomputedShares[bitIndex].t);
        secondParty = true;
    }

    verifySub(a.shares, b.shares, y.shares, skippedParty); // y = a-b

    verifyMulScalar(t.shares, t.shares, N, skippedParty); // t*N

    verifySub(y.shares, t.shares, y.shares, skippedParty); // y = a - b - t*N

    verifyZeroAssertion(y, "squareModNExtractLSB", v1, v2);

    verifyLowerThan(b, N, v1, v2, ole1, ole2, circuitDepth);

    verifyLowerThan(a, Nsquare, v1, v2, ole1, ole2, circuitDepth);

    verifyLessThanHalfOrder(t, v1, v2, ole1, ole2, circuitDepth);

    LayerOutput lo;
    lo.lsb = verifyComputeLSB(b, v1, v2, ole1, ole2, circuitDepth);
    lo.next = b;

    return lo;
}

void testBooleanizeLSB() {
    for (int i = 0; i < 100; i++) {
        for (int lsb = 0; lsb <= 1; lsb++) {
            mpz_t lsbN;
            mpz_init_set_si(lsbN, lsb);

            SS ss = secretShareNoTape(lsbN);

            RandomTape *tapes[3];
            for (int i = 0; i < 3; i++) {
                tapes[i] = newRandomTape();
            }

            Views views = createViews();

            char out = booleanizeLSB(ss, tapes, views);
            int b0 = out & 1;
            int b1 = (out & 2) > 1;
            int b2 = (out & 4) >> 2;

            int b[3] = {b0, b1, b2};

            if (b0 ^ b1 ^ b2 != lsb) {
                fprintf(stderr, "result is %d but expected 1\n");
                fflush(stderr);
                exit(2);
            }

            int skippedParty = rand() % 3;


            OpenedView openedViews[2];
            int index = 0;
            for (int j=0; j <= 2; j++) {
                // Reveal only parties that are not the skipped party
                if (j == skippedParty) {
                    continue;
                }

                openedViews[index].party = j;
                openedViews[index].tape = tapeFromSeed(tapes[j]->seed);
                openedViews[index].conversionShares = views.conversions[j];
                openedViews[index].reconstructions = views.reconstructions[j];
                openedViews[index].reconstructions->reconstructionIndex = 0;
                openedViews[index].oleOutputs = views.fromOLE[j];
                index++;
            }

            // clear skipped party
            mpz_clear(ss.shares[skippedParty]);


            ComputedToOLE toOLE1;
            ComputedToOLE toOLE2;
            int depth = 0;
            int convIndex = 0;
            out = verifyBooleanizeLSB(ss, &openedViews[0], &openedViews[1], &toOLE1, &toOLE2, &depth, &convIndex);
            int v0 = out & 1;
            int v1 = (out & 2) > 1;
            int v2 = (out & 4) >> 2;

            int v[3] = {v0, v1, v2};
            for (int j=0; j<3; j++) {
                if (j == skippedParty) {
                    continue;
                }
                if (b[j] != v[j]) {
                    fprintf(stderr, "lsb: %d, b%d: %d, v%d: %d\n", lsb, j, b[j], j, v[j]);
                    exit(2);
                }
            }
        } // for LSB 0 or 1
    } // for i from 0 to 100
}

BooleanShares verifyAlgorithm1(SS z, const char *S, OpenedView *v1, OpenedView *v2, ComputedToOLE *ole1, ComputedToOLE *ole2) {
    int circuitDepth = 0;
    int convIndex = 0;
    BooleanShares result;
    int skippedParty = findSkippedParty(v1, v2);
    for (int i = 0; i < SECRET_LENGTH; i++) {
        int S_i = S[i] - '0';
        LayerOutput currentLayer = verifySquareModNExtractLSB(z, v1, v2, ole1, ole2, &circuitDepth, i);
        SS booleanShares = verifyXORPlaintext(currentLayer.lsb.shares, S_i, v1, v2);
        z = currentLayer.next;
        result.shares[i] = verifyBooleanizeLSB(booleanShares, v1, v2, ole1, ole2, &circuitDepth, &convIndex);
    }
    return result;
}
