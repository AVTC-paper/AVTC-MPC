#include "avtc.h"

// WARNING:
// This code is for educational purposes only.
// Do not use this in any commercial or production environment.

RandomTape *proverTape;

int serializeFieldElement(mpz_t n, void *buff) {
    mpz_export(buff, NULL, 1, 1, 1, 0, n);
    return mpz_sizeinbase(n, 2) / 8;
}

void clearSS(SS s) {
    mpz_clear(s.shares[0]);
    mpz_clear(s.shares[1]);
    mpz_clear(s.shares[2]);
}

SS zero() {
    mpz_t zero;
    mpz_init_set_si(zero, 0);

    SS zeroShares;
    mpz_init_set(zeroShares.shares[0], zero);
    mpz_init_set(zeroShares.shares[1], zero);
    mpz_init_set(zeroShares.shares[2], zero);

    return zeroShares;
}

void randFieldElement(mpz_t out, gmp_randstate_t rand) {
    mpz_init(out);
    mpz_urandomm(out, rand, field);
}

void SHA256Digest(const void *input, const int inputLength, unsigned char *outputBuffer) {
    SHA256_CTX sha256;
    SHA256_Init(&sha256);
    SHA256_Update(&sha256, input, inputLength);
    SHA256_Final(outputBuffer, &sha256);
}

RandomTape *tapeFromSeed(const unsigned char *rawSeed) {
    RandomTape *tape = malloc(sizeof(RandomTape));
    tape->seed = rawSeed;
    tape->invocationNum = 0;

    mpz_t seed;
    mpz_init(seed);
    mpz_import(seed, RANDOM_TAPE_SEED_SIZE, 1, sizeof(char), 0, 0, rawSeed);

    gmp_randinit_mt(tape->r);
    gmp_randseed(tape->r, seed);

    return tape;
}

RandomTape *newRandomTape() {
    unsigned char *rawSeed = malloc(sizeof(char) * RANDOM_TAPE_SEED_SIZE);
    int success = RAND_bytes(rawSeed, sizeof(char) * RANDOM_TAPE_SEED_SIZE);
    if (!success) {
        fprintf(stderr, "failed getting randomness");
        exit(2);
    }
    return tapeFromSeed(rawSeed);
}

int nextRandomTapeEntry(mpz_t out, RandomTape *tape, mpz_t upperBound) {
    mpz_init(out);
    mpz_urandomm(out, tape->r, upperBound);
    int num = tape->invocationNum;
    tape->invocationNum++;
    return num;
}

void randBlindingFactor(RandomTape *tape, mpz_t out, mpz_t blindingFactorUpperbound) {
    mpz_init(out);
    nextRandomTapeEntry(out, tape, blindingFactorUpperbound);
}

void randEvenBlindingFactor(RandomTape *tape, mpz_t out, mpz_t blindingFactorUpperbound) {
    randBlindingFactor(tape, out, blindingFactorUpperbound);
    mpz_clrbit(out, 0);
}

SS secretShare(mpz_t s, RandomTape *tape) {
    mpz_t x1;
    mpz_t x2;
    mpz_t x3;
    mpz_init(x1);
    mpz_init(x2);
    mpz_init(x3);

    randFieldElement(x1, tape->r);
    randFieldElement(x2, tape->r);
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


    mpz_clear(x1);
    mpz_clear(x2);
    mpz_clear(x3);
    mpz_clear(sum);

    return ss;
}

SS secretShareNoTape(mpz_t w) {
    mpz_t x1;
    mpz_t x2;
    mpz_t x3;
    mpz_init(x1);
    mpz_init(x2);
    mpz_init(x3);

    randFieldElement(x1, globalRand);
    randFieldElement(x2, globalRand);
    // (x1+x2+x3)%field = w
    mpz_t sum;
    mpz_init(sum);
    mpz_add(sum, x1, x2);
    mpz_sub(sum, w, sum);     // sum = w - (x1+x2)
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

void mpc_Add(mpz_t x[3], mpz_t y[3], mpz_t z[3]) {
    mpz_add(z[0], x[0], y[0]);
    mpz_mod(z[0], z[0], field);

    mpz_add(z[1], x[1], y[1]);
    mpz_mod(z[1], z[1], field);

    mpz_add(z[2], x[2], y[2]);
    mpz_mod(z[2], z[2], field);
}

void mpc_AddScalar(mpz_t in[3], mpz_t out[3], mpz_t n) {
    mpz_add(out[0], in[0], n);
    mpz_mod(out[0], out[0], field);
    mpz_set(out[1], in[1]);
    mpz_set(out[2], in[2]);
}

SS neg(SS x) {
    SS res;
    mpz_init(res.shares[0]);
    mpz_init(res.shares[1]);
    mpz_init(res.shares[2]);
    mpz_sub(res.shares[0], field, x.shares[0]);
    mpz_sub(res.shares[1], field, x.shares[1]);
    mpz_sub(res.shares[2], field, x.shares[2]);
    return res;
}

void mpc_Sub(mpz_t x[3], mpz_t y[3], mpz_t z[3]) {
    mpz_t negY[3];
    mpz_init(negY[0]);
    mpz_init(negY[1]);
    mpz_init(negY[2]);
    mpz_sub(negY[0], field, y[0]);
    mpz_sub(negY[1], field, y[1]);
    mpz_sub(negY[2], field, y[2]);

    mpc_Add(x, negY, z);

}

void mpc_mulScalar(mpz_t in[3], mpz_t out[3], mpz_t n) {
    for (int i = 0; i < 3; i++) {
        mpz_mul(out[i], in[i], n);
        mpz_mod(out[i], out[i], field);
    }
}

void initRand(gmp_randstate_t rand) {
    void *randomBuff = malloc(32);
    if (!randomBuff) {
        fprintf(stderr, "malloc failed");
        exit(2);
    }
    if (getrandom(randomBuff, 32, GRND_RANDOM) == -1) {
        fprintf(stderr, "not enough randomness");
        exit(2);
    }

    mpz_t randomNumberOf256it;
    mpz_init(randomNumberOf256it);
    mpz_import(randomNumberOf256it, 32, 1, sizeof(randomBuff[0]), 0, 0, randomBuff);

    gmp_randinit_mt(rand);
    gmp_randseed(rand, randomNumberOf256it);
}

void randGroupElement(mpz_t out) {
    gmp_randstate_t r;
    initRand(r);

    mpz_init(out);
    mpz_urandomm(out, r, N);
}

int bitAt(int n, int i) {
    n = n >> i;
    return n & 1;
}

int toIndex(int i, int j) {
    if (i == 2) {
        return j;
    }

    if (i == 0) {
        return j - 1;
    }

    if (i == 1) {
        if (j == 0) {
            return 0;
        }

        if (j == 2) {
            return 1;
        }
    }

    fprintf(stderr, "invalid party number %d, expected one of 0,1,2", i);
    fflush(stderr);
    exit(2);
}

int fromIndex(int i, int j) {
    if (i == 2) {
        return 1;
    }

    if (i == 0) {
        return 0;
    }

    if (i == 1) {
        if (j == 0) {
            return 0;
        }

        if (j == 2) {
            return 1;
        }
    }

    fprintf(stderr, "invalid party number %d, expected one of 0,1,2", i);
    fflush(stderr);
    exit(2);
}

int partyIndex(int i, int j) {
    if (i == 0) {
        if (j == 0) {
            return 1;
        } else {
            return 2;
        }
    }

    if (i == 1) {
        if (j == 0) {
            return 0;
        } else {
            return 2;
        }
    }

    if (i == 2) {
        if (j == 0) {
            return 0;
        } else {
            return 1;
        }
    }

    fprintf(stderr, "invalid party number %d, expected one of 0,1,2", i);
    fflush(stderr);
    exit(2);
}

// out = x - y
void subModField(mpz_t out, mpz_t x, mpz_t y) {
    mpz_sub(out, x, y);
    while (mpz_sgn(out) == -1) {
        mpz_add(out, out, field);
    }
}

// 50 micro-seconds
void mpc_Mul(mpz_t x[3], mpz_t y[3], mpz_t z[3], RandomTape *tapes[3], Views views) {
    int depth = *(views.depth);
    for (int i = 0; i <= 2; i++) { // i is the current party
        nextRandomTapeEntry(views.toOLE[i]->msgs[depth].r0, tapes[i], field);
        nextRandomTapeEntry(views.toOLE[i]->msgs[depth].r1, tapes[i], field);

        for (int j = 0; j <= 1; j++) { // j is either 0 or 1 to indicate the two other parties
            mpz_t XiYj;
            mpz_init(XiYj);

            int jIndex = partyIndex(i, j);
            int from = fromIndex(i, jIndex);

            if (j == 0) {
                //gmp_printf(">>>> party: %d, j: %d msgs[%d].y0 = %Zd\n",i, j, depth, y[i]);
                mpz_init_set(views.toOLE[i]->msgs[depth].y0, y[i]);
            } else {
                //gmp_printf(">>>> party: %d, j: %d msgs[%d].y1 = %Zd\n",i, j, depth, y[i]);
                mpz_init_set(views.toOLE[i]->msgs[depth].y1, y[i]);
            }

            mpz_mul(XiYj, x[i], y[jIndex]);

            if (j == 0) {
                //gmp_printf(">>>> party: %d, j: %d msgs[%d].x0 = %Zd\n",i, j, depth, x[i]);
                mpz_init_set(views.toOLE[i]->msgs[depth].x0, x[i]);
                //mpz_init(views.toOLE[i]->msgs[depth].x1);
            } else {
               // gmp_printf(">>>> party: %d, j: %d msgs[%d].x1 = %Zd\n",i, j, depth, x[i]);
                mpz_init_set(views.toOLE[i]->msgs[depth].x1, x[i]);
                //mpz_init(views.toOLE[i]->msgs[depth].x0);
           }

            mpz_t Sij;
            mpz_init(Sij);
            if (j == 0) {
                mpz_add(Sij, XiYj, views.toOLE[i]->msgs[depth].r0);
            } else {
                mpz_add(Sij, XiYj, views.toOLE[i]->msgs[depth].r1);
            }
            mpz_mod(Sij, Sij, field);

            mpz_clear(XiYj);

            if (from == 0) {
                mpz_set(views.fromOLE[jIndex]->incMsgs0[depth], Sij);
            } else {
                mpz_set(views.fromOLE[jIndex]->incMsgs1[depth], Sij);
            }

            mpz_clear(Sij);
        }
    }

    for (int i = 0; i <= 2; i++) {
        mpz_t XiYi;
        mpz_init(XiYi);
        mpz_mul(XiYi, x[i], y[i]);


        mpz_t Ci;
        mpz_init(Ci);

        mpz_add(Ci, XiYi, views.fromOLE[i]->incMsgs0[depth]);
        mpz_add(Ci, Ci, views.fromOLE[i]->incMsgs1[depth]);
        mpz_mod(Ci, Ci, field);

        subModField(Ci, Ci, views.toOLE[i]->msgs[depth].r0);
        subModField(Ci, Ci, views.toOLE[i]->msgs[depth].r1);

        mpz_set(z[i], Ci);

        mpz_clear(XiYi);
        mpz_clear(Ci);

    }

    if ((*views.depth) == ARITH_MAXDEPTH) {
        fprintf(stderr, "reached max depth, aborting");
        exit(2);
    }

    (*views.depth)++;
}

void initField(void) {
    mpz_init(N);
    mpz_t p;
    mpz_t q;
    mpz_init_set_str(p,
                     "167200376276670642016090512684478075986881043511698258649784440165760531614287751164699732488362198044880237364322370046262248848122701174814756282450532401114183644058420976505502388938601057953905377039740769159837596703470655396816499117519462002074550864047837670311350505073407902935899043479344607595889",
                     10);
    mpz_init_set_str(q,
                     "164581833258898471926047213996711707410139611089624872054206990511537353292876008515121465025240010559853537311638051224316187598264005969799826997295111863901433482543033345672346474526338271863541080475645207694300844235381330374607894855688195616904312060080619907928573074377328119925104754166127982720343",
                     10);
    mpz_mul(N, p, q);
    mpz_mul(Nsquare, N, N);
    mpz_init_set_str(field,
                     "1388230939660854883604537775460115325908969642113808203392122833418608787083976885379292433466725198699739819530446473259960376244446879928087002313292339436598883666515497470356709807215123075367363809584922311420487946848086086294843211415786155612579694258190011298992370027553809649365306209731303237155248187850837214184428302701607764455594871979520937915578956295739808786107224256617701491664745682360489947465422743884606551322139486204517350502996545318082135837447865877182823768742541783545351126326169612408570232810242771524314659891155520769161512177772827674367114815610776650253905035122029561503382630808547589247953666980258740405923223575872886504640445421053757720377937620518373631304054308114772611991472952154112969153357563892340542722161024486536771300508924436079471707550816433780648819578144183252598347990119108523128170828449329527611105118838531013214893192232833667197602496412251835786759535794001110388310810942772219363783125417353030175293444167602300061253288156651122889081497487724958221700669028632746130132659646386854884958861174269403103954638030725972676692723736647396306571922665343860278282285264575747229952057585921131969839919880482416587768088726685064072221203416591623577499739483915460972457002980246502712369218037",
                     10);
    fieldSize = mpz_sizeinbase(field, 2);
    mpz_t two;
    mpz_init_set_si(two, 2);
    mpz_init(blindingFactorUpperbound);
    mpz_powm_ui(blindingFactorUpperbound, two, mpz_sizeinbase(N, 2) + 100, field);
    mpz_init(plaintextShareUpperBound);
    mpz_t nSquare;
    mpz_init(nSquare);
    mpz_mul(nSquare, N, N);
    mpz_mul_si(plaintextShareUpperBound, blindingFactorUpperbound, 2 * 3); // two blinding factors, three parties
    mpz_add(plaintextShareUpperBound, blindingFactorUpperbound, nSquare);  // N^2

    printf("Field is %d bits, blinding factor is %d bits, N is %d bits, plaintext upper bound is %d bits, running %d iterations \n",
           fieldSize, fieldSize - 100, mpz_sizeinbase(N, 2), mpz_sizeinbase(plaintextShareUpperBound, 2), ITERATIONS);
    fflush(stdout);
}


Views createViews() {
    Views views;
    views.depth = (malloc(sizeof(int)));
    *views.depth = 0;
    views.conversionIndex = malloc(sizeof(int));
    *views.conversionIndex = 0;
    for (int i = 0; i < 3; i++) {
        views.reconstructions[i] = malloc(sizeof(ReconstructionShares));
        views.reconstructions[i]->reconstructionIndex = 0;
        views.fromOLE[i] = malloc(sizeof(OLEOutputs));
        views.toOLE[i] = malloc(sizeof(OLEInputs));
        views.conversions[i] = malloc(sizeof(ConversionShares));
        views.extractionShares[i] = malloc(sizeof(AllLSBExtractionShares));
        views.extractionShares[i]->distributionDepth = 0;
        views.extractionShares[i]->lsbExtractionDepth = 0;
        for (int j = 0; j < 4609; j++) {
            mpz_init(views.extractionShares[i]->sharesFromParty0[j]);
            mpz_init(views.extractionShares[i]->sharesFromParty1[j]);
        }
        for (int j = 0; j < ARITH_MAXDEPTH; j++) {
            mpz_init(views.extractionShares[i]->lsbExtractionShares[j].blindingShare0);
            mpz_init(views.extractionShares[i]->lsbExtractionShares[j].blindingShare1);
            mpz_init(views.extractionShares[i]->lsbExtractionShares[j].blindedShare);
        }
    }

    for (int party = 0; party < 3; party++) {
        for (int i = 0; i < ARITH_MAXDEPTH; i++) {
            mpz_init(views.fromOLE[party]->incMsgs0[i]);
            mpz_init(views.fromOLE[party]->incMsgs1[i]);
        }
    }

    return views;
}

void findField() {
    mpz_t n;
    mpz_init(n);
    mpz_set_str(n,
                "27518144449192091323785136890288081002111163008131579389868979123058526612961814169576479652690566315688165545261084691759244102073192775762602945314299857048479498599462787617752344346867312145528536437275377969294777816873672780487846002110944989481483427204154700053455661402700812086819281516246557605984163654710569034300736631406420995589670835608735412121315932554722222442747005055168229165207085101172903012266786259853145875210694315201682190470477846286400598367501391121175336904167177167826762738669787017253904152441033114569221450919087114704937595412378227707650670216958267520527775834443735443469927",
                10);

    mpz_t two;
    mpz_init_set_si(two, 2);

    char bigNumber[5001];
    bigNumber[0] = '1';
    for (int i = 1; i < 5000; i++) {
        bigNumber[i] = '1';
    }
    bigNumber[5000] = '\0';
    mpz_t bigNum;
    mpz_init_set_str(bigNum, bigNumber, 2);

    int targetBitSize = mpz_sizeinbase(n, 2) * 2 + 120;

    mpz_init(field);
    mpz_powm_ui(field, two, targetBitSize, bigNum);

    mpz_t one;
    mpz_init_set_si(one, 1);

    mpz_add(field, field, one); // p = 2^{targetBitSize} + 1

    while (1) {
        mpz_add(field, field, two); // p += 2
        int probably = mpz_probab_prime_p(field, 50);
        if (probably == 0) {
            continue;
        }

        mpz_out_str(stdout, 10, field);
        printf("\nfield size is %d bit", mpz_sizeinbase(field, 2));
        break;
    }
}


void reconstructArithStandalone(SS ss, mpz_t x) {
    mpz_t sum;
    mpz_init(sum);
    mpz_add(sum, ss.shares[0], ss.shares[1]);
    mpz_add(sum, sum, ss.shares[2]);
    mpz_mod(sum, sum, field);
    mpz_set(x, sum);
}

void reconstructArith(Views views, SS ss, mpz_t x) {
    for (int i = 0; i < 3; i++) {
        int partyStorageIndex = (i + 1) % 3;

        // Simulate party i broadcasting its share and saving it in party (i+1)%3
        int index = views.reconstructions[partyStorageIndex]->reconstructionIndex;
        mpz_init(views.reconstructions[partyStorageIndex]->shares[index]);
        mpz_set(views.reconstructions[partyStorageIndex]->shares[index], ss.shares[i]);

        views.reconstructions[partyStorageIndex]->reconstructionIndex++;
    }

    mpz_t sums[3];
    // Each party collects all extractionShares from other partis
    for (int i = 0; i < 3; i++) {
        mpz_init(sums[i]);
        mpz_add(sums[i], ss.shares[0], ss.shares[1]);
        mpz_add(sums[i], sums[i], ss.shares[2]);
        mpz_mod(sums[i], sums[i], field);

    }

    // Sanity check - compare reconstruction results across parties
    for (int i = 0; i < 2; i++) {
        int eq = mpz_cmp(sums[i], sums[i + 1]);
        if (eq != 0) {
            fprintf(stderr, "reconstruction yielded different results across parties %d and %d", i, i + 1);
            fflush(stderr);
            exit(2);
        }
    }

    mpz_set(x, sums[0]); // Without loss of generality, use the input of the first party
}

void assertZero(Views views, SS ss, const char *caller) {
    mpz_t n;
    mpz_init(n);
    reconstructArith(views, ss, n);

    if (mpz_sgn(n) != 0) {
        fprintf(stderr, "reconstruction of '%s' has %d bits ( value=%d )\n", caller, mpz_sizeinbase(n, 2),
                mpz_get_si(n));
        fflush(stderr);
        exit(2);
    }

    mpz_clear(n);
}

void assertLowerThan(Views views, SS v, mpz_t u) {
    mpz_t vplain;
    mpz_init(vplain);
    reconstructArith(views, v, vplain);
    int cmp = mpz_cmp(vplain, u);
    if (cmp != -1) {
        fprintf(stderr, "v >= u");
        exit(2);
    }

    mpz_clear(vplain);
}

SS secretShareBit(int x, RandomTape *tape) {
    if (x > 1 || x < 0) {
        fprintf(stderr, "x is %d but should be 0 or 1\n", x);
        fflush(stderr);
        exit(2);
    }
    mpz_t n;
    mpz_init_set_si(n, x);
    return secretShare(n, tape);
}


SS mpc_xor(SS x, SS y, RandomTape *tapes[3], Views views) {
    SS xplusy;
    mpz_init(xplusy.shares[0]);
    mpz_init(xplusy.shares[1]);
    mpz_init(xplusy.shares[2]);

    // (x xor y) = x+y - 2*x*y

    mpc_Add(x.shares, y.shares, xplusy.shares);

    mpz_t xy[3];
    mpz_init(xy[0]);
    mpz_init(xy[1]);
    mpz_init(xy[2]);

    mpc_Mul(x.shares, y.shares, xy, tapes, views);

    SS twoXY;
    mpz_init(twoXY.shares[0]);
    mpz_init(twoXY.shares[1]);
    mpz_init(twoXY.shares[2]);
    mpc_Add(xy, xy, twoXY.shares); // 2xy

    SS twoXYNeg = neg(twoXY);

    SS z;
    mpz_init(z.shares[0]);
    mpz_init(z.shares[1]);
    mpz_init(z.shares[2]);
    mpc_Add(twoXYNeg.shares, xplusy.shares, z.shares);

    clearSS(xplusy);

    return z;
}

SS mpc_XORPlaintext(mpz_t x[3], int y) {
    SS zs = zero();

    // (x xor y) = x+y - 2*x*y
    // if y == 1:  (x xor y) = x+1 - 2*x
    // if y == 0:  (x xor y) = x
    if (y == 0) {
        for (int i = 0; i < 3; i++) {
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

    mpc_AddScalar(x, zs.shares, one);

    SS twoX;
    mpz_init(twoX.shares[0]);
    mpz_init(twoX.shares[1]);
    mpz_init(twoX.shares[2]);
    mpc_Add(x, x, twoX.shares);

    SS minusTwoX = neg(twoX);
    mpc_Add(zs.shares, minusTwoX.shares, zs.shares);
    return zs;
}

void testXorPlaintext() {
    mpz_t n;
    mpz_init_set_si(n, 1);

    SS ns = secretShareNoTape(n);

    SS result = mpc_XORPlaintext(ns.shares, 1);

    mpz_t out;
    mpz_init(out);
    reconstructArithStandalone(result, out);

    int res = mpz_get_si(out);
    printf("%d\n", res);
}

void testAdd() {
    mpz_t x;
    mpz_init_set_si(x, 100);

    mpz_t y;
    mpz_init_set_si(y, 200);

    SS xs = secretShareNoTape(x);
    SS ys = secretShareNoTape(y);

    SS zs;
    mpz_init(zs.shares[0]);
    mpz_init(zs.shares[1]);
    mpz_init(zs.shares[2]);

    mpc_Add(xs.shares, ys.shares, zs.shares);

    mpz_t z;
    mpz_init(z);
    reconstructArithStandalone(zs, z);
    mpz_out_str(stdout, 10, z);
}

void precomputeSquareModN(PrecomputedShares *result, mpz_t in) {
    mpz_t zSquare;
    mpz_init(zSquare);

    mpz_t t;
    mpz_init(t);

    mpz_t z;
    mpz_init(z);
    mpz_set(z, in);


    mpz_t one;
    mpz_init_set_ui(one, 1);

    mpz_t lsb;
    mpz_init(lsb);
    mpz_and(lsb, z, one);

    char str[2];
    mpz_get_str(str, 10, lsb);
    result->lsb = str[0] - '0';

    mpz_mul(zSquare, z, z); // zSquare = z^2
    mpz_mod(z, zSquare, N); // z = z^2 % N
    mpz_sub(t, zSquare, z); // t = z^2 - (z^2 % N)
    mpz_div(t, t, N);       // t = (z^2 - (z^2 % N)) / N

    mpz_init(result->t);
    mpz_init(result->b);

    mpz_set(result->t, t);
    mpz_set(result->b, z);

    mpz_clear(zSquare);
    mpz_clear(t);
    mpz_clear(z);
    mpz_clear(lsb);
    mpz_clear(one);
}

ModReducePrecomp computeCircuitPlaintext(mpz_t w) {
    ModReducePrecomp result;
    PrecomputedShares *elements = malloc(sizeof(PrecomputedShares) * (SECRET_LENGTH + 1));
    result.elements = elements;
    result.length = SECRET_LENGTH + 1;

    // The circuit
    mpz_t z;
    mpz_init(z);

    precomputeSquareModN(&result.elements[0], w);
    mpz_set(z, result.elements[0].b);

    result.elements[0].b_shares = secretShareNoTape(result.elements[0].b);
    result.elements[0].t_shares = secretShareNoTape(result.elements[0].t);

    for (int i = 1; i <= SECRET_LENGTH; i++) {
        precomputeSquareModN(&result.elements[i], z);
        result.elements[i].b_shares = secretShareNoTape(result.elements[i].b);
        result.elements[i].t_shares = secretShareNoTape(result.elements[i].t);

        mpz_set(z, result.elements[i].b);
    }
    return result;

    mpz_clear(z);
}

void splitA(mpz_t a, mpz_t a1, mpz_t a2, mpz_t a3) {
    mpz_t three;
    mpz_init_set_si(three, 3);

    mpz_div(a1, a, three);
    mpz_div(a2, a, three);

    mpz_t a1plusa2;
    mpz_init(a1plusa2);
    mpz_add(a1plusa2, a1, a2);
    mpz_sub(a3, a, a1plusa2);

    mpz_clear(three);
}

int extractLSBByReconstruction(Views views, SS v, mpz_t vPlaintext, mpz_t lsb) {
    mpz_init(vPlaintext);
    reconstructArithStandalone(v, vPlaintext);

    int plaintextLSB = mpz_get_si(vPlaintext) % 2;

    if (plaintextLSB < 0 || plaintextLSB > 1) {
        fprintf(stderr, "plaintextLSB is %d but should be a bit", plaintextLSB);
        fflush(stderr);
        exit(2);
    }
    mpz_init_set_si(lsb, plaintextLSB);

    return plaintextLSB;
}

SS computeLSB(SS v, RandomTape *tapes[3], Views views) {
    // Offline phase:
    mpz_t vPlaintext;
    mpz_t lsb;
    int plaintextLSB = extractLSBByReconstruction(views, v, vPlaintext, lsb);

    mpz_t a;
    mpz_init(a);

    mpz_sub(a, vPlaintext, lsb);

    mpz_t a0;
    mpz_init(a0);

    mpz_t a1;
    mpz_init(a1);

    mpz_t a2;
    mpz_init(a2);

    splitA(a, a0, a1, a2);

    int b0 = 0;
    int b1 = 0;
    int b2 = plaintextLSB;

    if (plaintextLSB > 1 || plaintextLSB < 0) {
        fprintf(stderr, "plaintextLSB is %d\n", plaintextLSB);
        exit(2);
    }

    int e = rand();

    int e01 = bitAt(e, 0);
    int e02 = bitAt(e, 1);
    int e10 = bitAt(e, 2);
    int e12 = bitAt(e, 3);
    int e20 = bitAt(e, 4);
    int e21 = e01 ^e02 ^e10 ^e12 ^e20; // Make sure they all sum to 0

    if (e01 ^ e02 ^ e10 ^ e12 ^ e20 ^ e21 != 0) {
        fprintf(stderr, "eij do not sum to 0\n");
        exit(2);
    }


    mpz_t r01;
    mpz_t r02;
    mpz_t r10;
    mpz_t r12;
    mpz_t r20;
    mpz_t r21;


    randEvenBlindingFactor(proverTape, r01, blindingFactorUpperbound);
    randEvenBlindingFactor(proverTape, r02, blindingFactorUpperbound);
    randEvenBlindingFactor(proverTape, r10, blindingFactorUpperbound);
    randEvenBlindingFactor(proverTape, r12, blindingFactorUpperbound);
    randEvenBlindingFactor(proverTape, r20, blindingFactorUpperbound);
    randEvenBlindingFactor(proverTape, r21, blindingFactorUpperbound);


    mpz_t blindedPlaintext[3];


    // Party 0 input:
    mpz_init(blindedPlaintext[0]);
    mpz_add(blindedPlaintext[0], a0, r01);
    mpz_add(blindedPlaintext[0], blindedPlaintext[0], r02);
    int b0blinded = e01 ^e02 ^b0;

    // Party 1 input:
    mpz_init(blindedPlaintext[1]);
    mpz_add(blindedPlaintext[1], a1, r10);
    mpz_add(blindedPlaintext[1], blindedPlaintext[1], r12);
    int b1blinded = e10 ^e12 ^b1;

    // Party 2 input:
    mpz_init(blindedPlaintext[2]);
    mpz_add(blindedPlaintext[2], a2, r20);
    mpz_add(blindedPlaintext[2], blindedPlaintext[2], r21);
    int b2blinded = e21 ^e20 ^b2;

    // Online phase:
    // Initialize input for each party
    mpz_t blindingFactors[3][2];

    mpz_init(blindingFactors[0][0]);
    mpz_init(blindingFactors[0][1]);
    mpz_init(blindingFactors[1][0]);
    mpz_init(blindingFactors[1][1]);
    mpz_init(blindingFactors[2][0]);
    mpz_init(blindingFactors[2][1]);

    mpz_set(blindingFactors[0][0], r01);
    mpz_set(blindingFactors[0][1], r02);
    mpz_set(blindingFactors[1][0], r10);
    mpz_set(blindingFactors[1][1], r12);
    mpz_set(blindingFactors[2][0], r20);
    mpz_set(blindingFactors[2][1], r21);

    int blindingBits[3][2] = {
            {e01, e02},
            {e10, e12},
            {e20, e21},
    };
    int blindedBits[3] = {b0blinded, b1blinded, b2blinded};


    for (int i = 0; i < 3; i++) {
        if (blindingBits[i][0] > 1 || blindingBits[i][0] < 0) {
            fprintf(stderr, "blinding bit [%d][0] is %d", i, blindingBits[i][0]);
            fflush(stderr);
            exit(2);
        }
        if (blindingBits[i][1] > 1 || blindingBits[i][1] < 0) {
            fprintf(stderr, "blinding bit [%d][1] is %d", i, blindingBits[i][1]);
            fflush(stderr);
            exit(2);
        }
    }

    // Input validation:
    // The verifier should check that blinded bits are indeed bits and blinded values and blinding factors
    // are within appropriate range. We do not check it here since we're the prover.

    // Secret sharing
    SS blindedPlaintextShares[3];
    SS blindingFactorsShares1[3];
    SS blindingFactorsShares2[3];
    SS blindedBitsShares[3];
    SS blindingBitsShares1[3];
    SS blindingBitsShares2[3];

    // Set shares given to parties by the prover
    for (int i = 0; i < 3; i++) {
        mpz_set(views.extractionShares[i]->lsbExtractionShares[views.extractionShares[i]->lsbExtractionDepth].blindedShare,
                blindedPlaintext[i]);
        mpz_set(views.extractionShares[i]->lsbExtractionShares[views.extractionShares[i]->lsbExtractionDepth].blindingShare0,
                blindingFactors[i][0]);
        mpz_set(views.extractionShares[i]->lsbExtractionShares[views.extractionShares[i]->lsbExtractionDepth].blindingShare1,
                blindingFactors[i][1]);
        views.extractionShares[i]->lsbExtractionShares[views.extractionShares[i]->lsbExtractionDepth].blindedBit = blindedBits[i];
        views.extractionShares[i]->lsbExtractionShares[views.extractionShares[i]->lsbExtractionDepth].blindingBit0 = blindingBits[i][0];
        views.extractionShares[i]->lsbExtractionShares[views.extractionShares[i]->lsbExtractionDepth].blindingBit1 = blindingBits[i][1];
        views.extractionShares[i]->lsbExtractionDepth++;
    }

    // Parties secret share their local input
    for (int i = 0; i < 3; i++) {
        blindedPlaintextShares[i] = secretShare(blindedPlaintext[i], tapes[i]);
        blindingFactorsShares1[i] = secretShare(blindingFactors[i][0], tapes[i]);
        blindingFactorsShares2[i] = secretShare(blindingFactors[i][1], tapes[i]);
        blindedBitsShares[i] = secretShareBit(blindedBits[i], tapes[i]);
        blindingBitsShares1[i] = secretShareBit(blindingBits[i][0], tapes[i]);
        blindingBitsShares2[i] = secretShareBit(blindingBits[i][1], tapes[i]);
    }

    int fromArr[3][2] = {
            {1, 2},
            {0, 2},
            {0, 1},
    };

    // Simulate secret share distribution among parties
    for (int i = 0; i < 3; i++) {
        int from0 = fromArr[i][0];
        int from1 = fromArr[i][1];
        mpz_set(views.extractionShares[i]->sharesFromParty0[views.extractionShares[i]->distributionDepth],
                blindedPlaintextShares[from0].shares[i]);
        mpz_set(views.extractionShares[i]->sharesFromParty1[views.extractionShares[i]->distributionDepth],
                blindedPlaintextShares[from1].shares[i]);
        views.extractionShares[i]->distributionDepth++;
        mpz_set(views.extractionShares[i]->sharesFromParty0[views.extractionShares[i]->distributionDepth],
                blindingFactorsShares1[from0].shares[i]);
        mpz_set(views.extractionShares[i]->sharesFromParty1[views.extractionShares[i]->distributionDepth],
                blindingFactorsShares1[from1].shares[i]);
        views.extractionShares[i]->distributionDepth++;
        mpz_set(views.extractionShares[i]->sharesFromParty0[views.extractionShares[i]->distributionDepth],
                blindingFactorsShares2[from0].shares[i]);
        mpz_set(views.extractionShares[i]->sharesFromParty1[views.extractionShares[i]->distributionDepth],
                blindingFactorsShares2[from1].shares[i]);
        views.extractionShares[i]->distributionDepth++;
        mpz_set(views.extractionShares[i]->sharesFromParty0[views.extractionShares[i]->distributionDepth],
                blindedBitsShares[from0].shares[i]);
        mpz_set(views.extractionShares[i]->sharesFromParty1[views.extractionShares[i]->distributionDepth],
                blindedBitsShares[from1].shares[i]);
        views.extractionShares[i]->distributionDepth++;
        mpz_set(views.extractionShares[i]->sharesFromParty0[views.extractionShares[i]->distributionDepth],
                blindingBitsShares1[from0].shares[i]);
        mpz_set(views.extractionShares[i]->sharesFromParty1[views.extractionShares[i]->distributionDepth],
                blindingBitsShares1[from1].shares[i]);
        views.extractionShares[i]->distributionDepth++;
        mpz_set(views.extractionShares[i]->sharesFromParty0[views.extractionShares[i]->distributionDepth],
                blindingBitsShares2[from0].shares[i]);
        mpz_set(views.extractionShares[i]->sharesFromParty1[views.extractionShares[i]->distributionDepth],
                blindingBitsShares2[from1].shares[i]);
        views.extractionShares[i]->distributionDepth++;
    }

    SS v0 = zero();
    for (int i = 0; i < 3; i++) {
        v0 = mpc_xor(v0, blindedBitsShares[i], tapes, views);
        v0 = mpc_xor(v0, blindingBitsShares1[i], tapes, views);
        v0 = mpc_xor(v0, blindingBitsShares2[i], tapes, views);
    }

    SS sumBlinded = zero();
    for (int i = 0; i < 3; i++) {
        mpc_Add(sumBlinded.shares, blindedPlaintextShares[i].shares, sumBlinded.shares);
    }

    SS sumBlindingFactors = zero();
    for (int i = 0; i < 3; i++) {
        mpc_Add(sumBlindingFactors.shares, blindingFactorsShares1[i].shares, sumBlindingFactors.shares);
        mpc_Add(sumBlindingFactors.shares, blindingFactorsShares2[i].shares, sumBlindingFactors.shares);
    }

    SS shouldBeZero = zero();
    mpc_Sub(v.shares, sumBlinded.shares, shouldBeZero.shares);
    mpc_Add(shouldBeZero.shares, sumBlindingFactors.shares, shouldBeZero.shares);
    mpc_Sub(shouldBeZero.shares, v0.shares, shouldBeZero.shares);


    assertZero(views, shouldBeZero, "computeLSB");


   mpz_clears(vPlaintext, lsb, a, a0, a1, a2, r01, r02, r10, r12, r20, r21, blindedPlaintext[0], blindedPlaintext[1], blindedPlaintext[2],
             blindingFactors[0][0], blindingFactors[0][1], blindingFactors[1][0], blindingFactors[1][1], blindingFactors[2][0], blindingFactors[2][1], NULL);

    clearSS(shouldBeZero);
    clearSS(sumBlindingFactors);
    clearSS(sumBlinded);

    clearSS(blindedPlaintextShares[0]);
    clearSS(blindedPlaintextShares[1]);
    clearSS(blindedPlaintextShares[2]);
    clearSS(blindingFactorsShares1[0]);
    clearSS(blindingFactorsShares1[1]);
    clearSS(blindingFactorsShares1[2]);
    clearSS(blindingFactorsShares2[0]);
    clearSS(blindingFactorsShares2[1]);
    clearSS(blindingFactorsShares2[2]);
    clearSS(blindedBitsShares[0]);
    clearSS(blindedBitsShares[1]);
    clearSS(blindedBitsShares[2]);
    clearSS(blindingBitsShares1[0]);
    clearSS(blindingBitsShares1[1]);
    clearSS(blindingBitsShares1[2]);
    clearSS(blindingBitsShares2[0]);
    clearSS(blindingBitsShares2[1]);
    clearSS(blindingBitsShares2[2]);

    return v0;
}

void lessThanHalfOrder(SS v, RandomTape *tapes[3], Views views) {
    SS w;
    mpz_init(w.shares[0]);
    mpz_init(w.shares[1]);
    mpz_init(w.shares[2]);
    mpc_Add(v.shares, v.shares, w.shares);
    SS shouldBeZero = computeLSB(w, tapes, views);
    assertZero(views, shouldBeZero, "lessThanHalfOrder");

    clearSS(w);
    clearSS(shouldBeZero);
}

void lowerThan(SS v, mpz_t u, RandomTape *tapes[3], Views views) {
    lessThanHalfOrder(v, tapes, views);
    SS w = neg(v);
    SS x = zero();
    mpc_AddScalar(w.shares, x.shares, u);
    lessThanHalfOrder(x, tapes, views);

    clearSS(w);
    clearSS(x);
}

LayerOutput squareModNExtractLSB(SS x_i, SS b, SS t, RandomTape *tapes[3], Views views) {
    SS a;
    mpz_init(a.shares[0]);
    mpz_init(a.shares[1]);
    mpz_init(a.shares[2]);

    mpc_Mul(x_i.shares, x_i.shares, a.shares, tapes, views); // a = x^2

    SS y;
    mpz_init(y.shares[0]);
    mpz_init(y.shares[1]);
    mpz_init(y.shares[2]);
    mpc_Sub(a.shares, b.shares, y.shares); // y = a-b

    mpc_mulScalar(t.shares, t.shares, N); // t*N

    mpc_Sub(y.shares, t.shares, y.shares); // y = a - b - t*N
    assertZero(views, y, "squareModNExtractLSB");

    lowerThan(b, N, tapes, views);
    lowerThan(a, Nsquare, tapes, views);
    lessThanHalfOrder(t, tapes, views);

    LayerOutput lo;
    lo.lsb = computeLSB(b, tapes, views);
    lo.next = b;


    clearSS(a);
    clearSS(y);

    return lo;
}

char booleanizeLSB(SS in, RandomTape *tapes[3], Views views) {
    int conversionIndex = *views.conversionIndex;
    (*views.conversionIndex)++;

    mpz_t ourShare[3];
    char ourBit[3];

    for (int i = 0; i < 3; i++) {
        // Create blinding factors
        mpz_t r_i;
        randBlindingFactor(tapes[i], r_i, blindingFactorUpperbound);

        int lsb = mpz_tstbit(r_i, 0);

        // Secret share the blinding factors in both domains;
        SS rs = secretShare(r_i, tapes[i]);

        mpz_init_set(ourShare[i], rs.shares[i]);

        mpz_t tmp;
        randBlindingFactor(tapes[i], tmp, blindingFactorUpperbound);
        //gmp_printf("randBlindingFactor[%d]=%Zd\n",i, tmp);
        int b0 = mpz_tstbit(tmp, 0); // random bit for party 0
        int b1 = mpz_tstbit(tmp, 1); // random bit for party 1
        ourBit[i] = b0 ^ b1 ^ lsb; // We get the XOR of the rest of the bits

        // Send the blinding factors to the other parties;
        for (int j = 0; j <= 1; j++) {
            int jIndex = partyIndex(i, j);
            int from = fromIndex(i, jIndex);
            if (from == 0) {
                mpz_init_set(views.conversions[jIndex]->shares[conversionIndex].r_0, rs.shares[jIndex]);
                views.conversions[jIndex]->shares[conversionIndex].b_0 = (j == 0 ? b0 : b1);
            } else {
                mpz_init_set(views.conversions[jIndex]->shares[conversionIndex].r_1, rs.shares[jIndex]);
                views.conversions[jIndex]->shares[conversionIndex].b_1 = (j == 0 ? b0 : b1);
            }
        }
    }

    // XOR all boolean domain blinding factors and create a sharing represented by 3 bits
    char lsbP0 = ourBit[0] ^views.conversions[0]->shares[conversionIndex].b_0 ^
                 views.conversions[0]->shares[conversionIndex].b_1;
    char lsbP1 = ourBit[1] ^views.conversions[1]->shares[conversionIndex].b_0 ^
                 views.conversions[1]->shares[conversionIndex].b_1;
    char lsbP2 = ourBit[2] ^views.conversions[2]->shares[conversionIndex].b_0 ^
                 views.conversions[2]->shares[conversionIndex].b_1;


    // Sum all arithmetic blinding factors and create a sharing 'r'
    SS r;
    for (int i = 0; i < 3; i++) {
        mpz_init(r.shares[i]);
        mpz_add(r.shares[i], ourShare[i], views.conversions[i]->shares[conversionIndex].r_0);
        mpz_mod(r.shares[i], r.shares[i], field);
        mpz_add(r.shares[i], r.shares[i], views.conversions[i]->shares[conversionIndex].r_1);
        mpz_mod(r.shares[i], r.shares[i], field);
    }

    SS x;
    mpz_init(x.shares[0]);
    mpz_init(x.shares[1]);
    mpz_init(x.shares[2]);

    // Blind arithmetic secret with blinding factor
    mpc_Add(r.shares, in.shares, x.shares);

    // Reconstruct and reveal blinded secret
    mpz_t xPlaintext;
    mpz_init(xPlaintext);

    reconstructArith(views, x, xPlaintext);

    int lsb = mpz_tstbit(xPlaintext, 0);
    // Without loss of generality, add LSB to party 0's share
    lsbP0 = lsbP0 ^ lsb;

    // Encode output by putting party i's sharing in bit i.
    char out = (lsbP2 << 2) | (lsbP1 << 1) | (lsbP0);

    mpz_clear(r.shares[0]);
    mpz_clear(r.shares[1]);
    mpz_clear(r.shares[2]);

    mpz_clear(xPlaintext);

    mpz_clear(x.shares[0]);
    mpz_clear(x.shares[1]);
    mpz_clear(x.shares[2]);


    mpz_clear(ourShare[0]);
    mpz_clear(ourShare[1]);
    mpz_clear(ourShare[2]);

    return out;
}

BooleanShares algorithm1(SS z, Views views, RandomTape *tapes[3], ModReducePrecomp modReducePrecomp, const char *S) {
    BooleanShares result;
    for (int i = 0; i < SECRET_LENGTH; i++) {
        int S_i = S[i] - '0';
        SS b = modReducePrecomp.elements[i].b_shares;
        SS t = modReducePrecomp.elements[i].t_shares;
        LayerOutput currentLayer = squareModNExtractLSB(z, b, t, tapes, views);
        SS booleanShares = mpc_XORPlaintext(currentLayer.lsb.shares, S_i);
        z = currentLayer.next;
        result.shares[i] = booleanizeLSB(booleanShares, tapes, views);
    }

    return result;
}

Algorithm1Output runAlgorithm1(char *S, mpz_t w, SS z) {
    Algorithm1Output result;
    omp_set_num_threads(ITERATIONS);


#pragma omp parallel for default(none) shared(result, w, z, S)
    for (int i = 0; i < ITERATIONS; i++) {
        ModReducePrecomp modReducePrecomp = computeCircuitPlaintext(w);
        for (int j = 0; j < SECRET_LENGTH + 1; j++) {
            for (int party = 0; party < 3; party++) {
                mpz_init_set(result.simulations[i].committedViews[party].precomputedShares[j].b,
                             modReducePrecomp.elements[j].b_shares.shares[party]);
                mpz_init_set(result.simulations[i].committedViews[party].precomputedShares[j].t,
                             modReducePrecomp.elements[j].t_shares.shares[party]);
            }
        }
        Views views = createViews();
        RandomTape *tapes[3] = {newRandomTape(), newRandomTape(), newRandomTape()};
        result.booleanShares[i] = algorithm1(z, views, tapes, modReducePrecomp, S);

        for (int party = 0; party < 3; party++) {
            result.simulations[i].committedViews[party].randomnessSeed = tapes[party]->seed;
            result.simulations[i].committedViews[party].oleOutputs = views.fromOLE[party];
            result.simulations[i].committedViews[party].oleInputs = views.toOLE[party];
            result.simulations[i].committedViews[party].conversionShares = views.conversions[party];
            result.simulations[i].committedViews[party].reconstructions = views.reconstructions[party];
            result.simulations[i].committedViews[party].lsbExtractionShares = views.extractionShares[party];
        }
    }
    return result;

}

OpenedViewTuple *openView(Algorithm1Output *out, int i, int nonRevealedParty) {
    OpenedViewTuple *viewTuple;
    viewTuple = malloc(sizeof(OpenedViewTuple));
    int n1 = (nonRevealedParty + 1) % 3;
    int n2 = (nonRevealedParty + 2) % 3;
    int p0 = MIN(n1, n2);
    int p1 = MAX(n1, n2);
    int *depth = malloc(sizeof(int));
    *(depth) = 0;
    viewTuple->depth = depth;
    viewTuple->views[0].party = p0;
    viewTuple->views[1].party = p1;
    viewTuple->views[0].oleOutputs = out->simulations[i].committedViews[p0].oleOutputs;
    viewTuple->views[1].oleOutputs = out->simulations[i].committedViews[p1].oleOutputs;
    viewTuple->views[0].conversionShares = out->simulations[i].committedViews[p0].conversionShares;
    viewTuple->views[1].conversionShares = out->simulations[i].committedViews[p1].conversionShares;
    viewTuple->views[0].lsbExtractionShares = out->simulations[i].committedViews[p0].lsbExtractionShares;
    viewTuple->views[1].lsbExtractionShares = out->simulations[i].committedViews[p1].lsbExtractionShares;
    viewTuple->views[0].lsbExtractionShares->lsbExtractionDepth = 0;
    viewTuple->views[1].lsbExtractionShares->lsbExtractionDepth = 0;
    viewTuple->views[0].lsbExtractionShares->distributionDepth = 0;
    viewTuple->views[1].lsbExtractionShares->distributionDepth = 0;
    viewTuple->views[0].reconstructions = out->simulations[i].committedViews[p0].reconstructions;
    viewTuple->views[1].reconstructions = out->simulations[i].committedViews[p1].reconstructions;
    viewTuple->views[0].reconstructions->reconstructionIndex = 0;
    viewTuple->views[1].reconstructions->reconstructionIndex = 0;
    viewTuple->views[0].tape = tapeFromSeed(out->simulations[i].committedViews[p0].randomnessSeed);
    viewTuple->views[1].tape = tapeFromSeed(out->simulations[i].committedViews[p1].randomnessSeed);
    for (int j = 0; j < SECRET_LENGTH + 1; j++) {
        mpz_init(viewTuple->views[0].precomputedShares[j].b);
        mpz_set(viewTuple->views[0].precomputedShares[j].b,
                out->simulations[i].committedViews[p0].precomputedShares[j].b);
        mpz_init(viewTuple->views[1].precomputedShares[j].b);
        mpz_set(viewTuple->views[1].precomputedShares[j].b,
                out->simulations[i].committedViews[p1].precomputedShares[j].b);
        mpz_init(viewTuple->views[0].precomputedShares[j].t);
        mpz_set(viewTuple->views[0].precomputedShares[j].t,
                out->simulations[i].committedViews[p0].precomputedShares[j].t);
        mpz_init(viewTuple->views[1].precomputedShares[j].t);
        mpz_set(viewTuple->views[1].precomputedShares[j].t,
                out->simulations[i].committedViews[p1].precomputedShares[j].t);
    }


    return viewTuple;
}

void openViews(OpenedViewTuple *openedViews[ITERATIONS], int *es, Algorithm1Output *out) {
    for (int i = 0; i < ITERATIONS; i++) {
        openedViews[i] = openView(out, i, es[i]);
    }
}

void testTapes() {
    RandomTape *t1 = newRandomTape();

    mpz_t a;
    nextRandomTapeEntry(a, t1, field);

    RandomTape *t2 = tapeFromSeed(t1->seed);

    mpz_t b;
    nextRandomTapeEntry(b, t2, field);

    if (mpz_cmp(a, b) != 0) {
        fprintf(stderr, "random elements do not match!\n");
        exit(2);
    }
}

int countProofSize(OpenedViewTuple *ov[ITERATIONS]) {
    long size = 0;
    for (int iteration = 0; iteration < ITERATIONS; iteration++) {
        for (int party = 0; party < 2; party++) {
            for (int i = 0; i < SECRET_LENGTH + 1; i++) {
                size += mpz_sizeinbase(ov[iteration]->views[party].precomputedShares[i].b, 2);
                size += mpz_sizeinbase(ov[iteration]->views[party].precomputedShares[i].t, 2);
            }
            for (int i = 0; i < 4609; i++) {
                // We only need 1 of the shares, as we can compute the other share ourselves
                size += mpz_sizeinbase(ov[iteration]->views[party].lsbExtractionShares->sharesFromParty0[i], 2);
            }
            for (int i = 0; i < ARITH_MAXDEPTH; i++) {
                size += mpz_sizeinbase(
                        ov[iteration]->views[party].lsbExtractionShares->lsbExtractionShares[i].blindingShare0, 2);
                size += mpz_sizeinbase(
                        ov[iteration]->views[party].lsbExtractionShares->lsbExtractionShares[i].blindingShare1, 2);
                size += mpz_sizeinbase(
                        ov[iteration]->views[party].lsbExtractionShares->lsbExtractionShares[i].blindedShare, 2);
            }
            for (int i = 0; i < SECRET_LENGTH; i++) {
                // We only need 1 of the shares, as we can compute the other share ourselves
                size += mpz_sizeinbase(ov[iteration]->views[party].conversionShares->shares[i].r_0, 2);
            }
            for (int i = 0; i < 1664; i++) {
                size += mpz_sizeinbase(ov[iteration]->views[party].reconstructions->shares[i], 2);
            }
            for (int i = 0; i < ARITH_MAXDEPTH; i++) {
                size += mpz_sizeinbase(ov[iteration]->views[party].oleOutputs->incMsgs0[i], 2);
                size += mpz_sizeinbase(ov[iteration]->views[party].oleOutputs->incMsgs1[i], 2);
            }
        }
    }


    size /= 8;
    return size;
}


void convertBoolShares(BooleanShares inputToSHA256, unsigned char inputs[2][SECRET_LENGTH / 8], int skippedParty) {
    for (int j = 0; j < SECRET_LENGTH / 8; j++) {
        char c0 = 0;
        char c1 = 0;
        for (int k = 0; k < 8; k++) {
            char currentPartiesShares = inputToSHA256.shares[j * 8 + k];
            bool secondParty = false;
            char bitP0;
            char bitP1;
            for (int p = 0; p < 3; p++) {
                if (p == skippedParty) {
                    continue;
                }
                if (secondParty) {
                    bitP1 = (currentPartiesShares >> p) & 1;
                } else {
                    bitP0 = (currentPartiesShares >> p) & 1;
                }
                secondParty = true;
            }
            c0 |= (bitP0 << k);
            c1 |= (bitP1 << k);
        }
        inputs[0][j] = c0;
        inputs[1][j] = c1;
    }
}

void
bitSharesTo8bitWords(BooleanShares booleanShares[ITERATIONS], unsigned char shares[ITERATIONS][3][SECRET_LENGTH / 8]) {
    for (int iteration = 0; iteration < ITERATIONS; iteration++) {
        for (int j = 0; j < SECRET_LENGTH / 8; j++) {
            shares[iteration][0][j] = 0;
            shares[iteration][1][j] = 0;
            shares[iteration][2][j] = 0;
            int offset = j * 8;
            for (int k = 0; k < 8; k++) {
                char sharesForCurrentBit = booleanShares[iteration].shares[offset + k];
                shares[iteration][0][j] |= ((sharesForCurrentBit & 1) << k);
                shares[iteration][1][j] |= (((sharesForCurrentBit & 2) >> 1) << k);
                shares[iteration][2][j] |= (((sharesForCurrentBit & 4) >> 2) << k);
            }
        }
    }
}


void compare8bitWords(int iteration, int skippedParty, unsigned char sha256Input[2][SECRET_LENGTH / 8],  unsigned char eightBitWords[ITERATIONS][3][SECRET_LENGTH / 8]) {
    int index = 0;
    for (int p = 0; p < 3; p++) {
        if (skippedParty == p) {
            continue;
        }

        for (int k = 0; k < SECRET_LENGTH / 8; k++) {
            bool eq = sha256Input[index][k] == eightBitWords[iteration][p][k];
            if (!eq) {
                printf("not equal 2 (%d) is %d vs %d ", k, sha256Input[index][k], eightBitWords[iteration][p][k]);
                exit(2);
            }
        }

        index++;
    }
}

void commitToOLEInputs(Algorithm1Output *out, int iteration, unsigned char *c1, unsigned char *c2, unsigned char *c3) {
    SHA256_CTX ctx1;
    SHA256_Init(&ctx1);

    SHA256_CTX ctx2;
    SHA256_Init(&ctx2);

    SHA256_CTX ctx3;
    SHA256_Init(&ctx3);

    int size = 528 * ARITH_MAXDEPTH * 6;

    unsigned char *buff = calloc(size, 1);

    for (int party = 0; party <= 2; party++) {
        int pos = 0;

        for (int i = 0; i < ARITH_MAXDEPTH; i++) {
            pos += serializeFieldElement(out->simulations[iteration].committedViews[party].oleInputs->msgs[i].r0, buff + pos);
            pos += serializeFieldElement(out->simulations[iteration].committedViews[party].oleInputs->msgs[i].r1, buff + pos);
            pos += serializeFieldElement(out->simulations[iteration].committedViews[party].oleInputs->msgs[i].x0, buff + pos);
            pos += serializeFieldElement(out->simulations[iteration].committedViews[party].oleInputs->msgs[i].x1, buff + pos);
            pos += serializeFieldElement(out->simulations[iteration].committedViews[party].oleInputs->msgs[i].y0, buff + pos);
            pos += serializeFieldElement(out->simulations[iteration].committedViews[party].oleInputs->msgs[i].y1, buff + pos);
        }

        if (party == 0) {
            SHA256_Update(&ctx1, buff, pos);
            SHA256_Final(c1, &ctx1);
        }

        if (party == 1) {
            SHA256_Update(&ctx2, buff, pos);
            SHA256_Final(c2, &ctx2);
        }

        if (party == 2) {
            SHA256_Update(&ctx3, buff, pos);
            SHA256_Final(c3, &ctx3);
        }

        memset(buff, 0, size);
    }

    free(buff);
}

void commitToViews(Algorithm1Output *out, int iteration, unsigned char *c1, unsigned char *c2, unsigned char *c3) {
    SHA256_CTX ctx1;
    SHA256_Init(&ctx1);

    SHA256_CTX ctx2;
    SHA256_Init(&ctx2);

    SHA256_CTX ctx3;
    SHA256_Init(&ctx3);

    int size = 528 * ((SECRET_LENGTH + 1) * 2 + 4609 * 2 + ARITH_MAXDEPTH * 3 + SECRET_LENGTH * 2 + 1664 + ARITH_MAXDEPTH * 2 );

    unsigned char *buff = calloc(size, 1);

    for (int party = 0; party <= 2; party++) {
        int pos = 0;

        for (int i = 0; i < SECRET_LENGTH + 1; i++) {
            pos += serializeFieldElement(out->simulations[iteration].committedViews[party].precomputedShares[i].b, buff + pos);
            pos += serializeFieldElement(out->simulations[iteration].committedViews[party].precomputedShares[i].t, buff + pos);
        }
        for (int i = 0; i < 4609; i++) {
            pos += serializeFieldElement(out->simulations[iteration].committedViews[party].lsbExtractionShares->sharesFromParty0[i], buff + pos);
            pos += serializeFieldElement(out->simulations[iteration].committedViews[party].lsbExtractionShares->sharesFromParty1[i], buff + pos);
        }
        for (int i = 0; i < ARITH_MAXDEPTH; i++) {
            pos += serializeFieldElement(out->simulations[iteration].committedViews[party].lsbExtractionShares->lsbExtractionShares[i].blindingShare0, buff + pos);
            pos += serializeFieldElement(out->simulations[iteration].committedViews[party].lsbExtractionShares->lsbExtractionShares[i].blindingShare1, buff + pos);
            pos += serializeFieldElement(out->simulations[iteration].committedViews[party].lsbExtractionShares->lsbExtractionShares[i].blindedShare, buff + pos);
        }
        for (int i = 0; i < SECRET_LENGTH; i++) {
            pos += serializeFieldElement(out->simulations[iteration].committedViews[party].conversionShares->shares[i].r_0, buff + pos);
            pos += serializeFieldElement(out->simulations[iteration].committedViews[party].conversionShares->shares[i].r_1, buff + pos);
        }
        for (int i = 0; i < 1664; i++) {
            pos += serializeFieldElement(out->simulations[iteration].committedViews[party].reconstructions->shares[i], buff + pos);
        }
        for (int i = 0; i < ARITH_MAXDEPTH; i++) {
            pos += serializeFieldElement(out->simulations[iteration].committedViews[party].oleOutputs->incMsgs0[i], buff + pos);
            pos += serializeFieldElement(out->simulations[iteration].committedViews[party].oleOutputs->incMsgs1[i], buff + pos);
        }

        if (party == 0) {
            SHA256_Update(&ctx1, buff, pos);
            SHA256_Final(c1, &ctx1);
        }

        if (party == 1) {
            SHA256_Update(&ctx2, buff, pos);
            SHA256_Final(c2, &ctx2);
        }

        if (party == 2) {
            SHA256_Update(&ctx3, buff, pos);
            SHA256_Final(c3, &ctx3);
        }

        memset(buff, 0, size);
    }

    free(buff);
}

void verifyCommittedOLEInputs(ComputedToOLE ole1, ComputedToOLE ole2, unsigned char *c1, unsigned char *c2) {
    SHA256_CTX ctx1;
    SHA256_Init(&ctx1);

    SHA256_CTX ctx2;
    SHA256_Init(&ctx2);

    unsigned char expected1[SHA256_DIGEST_LENGTH];
    unsigned char expected2[SHA256_DIGEST_LENGTH];

    int size = 528 * ARITH_MAXDEPTH * 6;

    unsigned char *buff = calloc(size, 1);

    for (int party = 0; party <= 1; party++) {
        int pos = 0;

        if (party == 1) {
            for (int i = 0; i < ARITH_MAXDEPTH; i++) {
                pos += serializeFieldElement(ole2.msgs[i].r0, buff + pos);
                pos += serializeFieldElement(ole2.msgs[i].r1, buff + pos);
                pos += serializeFieldElement(ole2.msgs[i].x0, buff + pos);
                pos += serializeFieldElement(ole2.msgs[i].x1, buff + pos);
                pos += serializeFieldElement(ole2.msgs[i].y0, buff + pos);
                pos += serializeFieldElement(ole2.msgs[i].y1, buff + pos);
            }
        } else {
            for (int i = 0; i < ARITH_MAXDEPTH; i++) {
                pos += serializeFieldElement(ole1.msgs[i].r0, buff + pos);
                pos += serializeFieldElement(ole1.msgs[i].r1, buff + pos);
                pos += serializeFieldElement(ole1.msgs[i].x0, buff + pos);
                pos += serializeFieldElement(ole1.msgs[i].x1, buff + pos);
                pos += serializeFieldElement(ole1.msgs[i].y0, buff + pos);
                pos += serializeFieldElement(ole1.msgs[i].y1, buff + pos);
            }
        }

        if (party == 0) {
            SHA256_Update(&ctx1, buff, pos);
            SHA256_Final(expected1, &ctx1);
            memset(buff, 0, size);
        } else {
            SHA256_Update(&ctx2, buff, pos);
            SHA256_Final(expected2, &ctx2);
        }
    }

    if (memcmp(expected1, c1, SHA256_DIGEST_LENGTH) != 0 || memcmp(expected2, c2, SHA256_DIGEST_LENGTH) != 0) {
        fprintf(stderr, "commitment to OLEs do not match");
        exit(2);
    }

    free(buff);
}

void verifyCommittedViewPair(OpenedViewTuple *ov, unsigned char *c1, unsigned char *c2) {
    SHA256_CTX ctx1;
    SHA256_Init(&ctx1);

    SHA256_CTX ctx2;
    SHA256_Init(&ctx2);

    unsigned char expected1[SHA256_DIGEST_LENGTH];
    unsigned char expected2[SHA256_DIGEST_LENGTH];

    int size = 528 * ((SECRET_LENGTH + 1) * 2 + 4609 * 2 + ARITH_MAXDEPTH * 3 + SECRET_LENGTH * 2 + 1664 + ARITH_MAXDEPTH * 2 );

    unsigned char *buff = calloc(size, 1);

    for (int party = 0; party <= 1; party++) {
        int pos = 0;

        for (int i = 0; i < SECRET_LENGTH + 1; i++) {
            pos += serializeFieldElement(ov->views[party].precomputedShares[i].b, buff + pos);
            pos += serializeFieldElement(ov->views[party].precomputedShares[i].t, buff + pos);
        }
        for (int i = 0; i < 4609; i++) {
            pos += serializeFieldElement(ov->views[party].lsbExtractionShares->sharesFromParty0[i], buff + pos);
            pos += serializeFieldElement(ov->views[party].lsbExtractionShares->sharesFromParty1[i], buff + pos);
        }
        for (int i = 0; i < ARITH_MAXDEPTH; i++) {
            pos += serializeFieldElement(ov->views[party].lsbExtractionShares->lsbExtractionShares[i].blindingShare0, buff + pos);
            pos += serializeFieldElement(ov->views[party].lsbExtractionShares->lsbExtractionShares[i].blindingShare1, buff + pos);
            pos += serializeFieldElement(ov->views[party].lsbExtractionShares->lsbExtractionShares[i].blindedShare, buff + pos);
        }
        for (int i = 0; i < SECRET_LENGTH; i++) {
            pos += serializeFieldElement(ov->views[party].conversionShares->shares[i].r_0, buff + pos);
            pos += serializeFieldElement(ov->views[party].conversionShares->shares[i].r_1, buff + pos);
        }
        for (int i = 0; i < 1664; i++) {
            pos += serializeFieldElement(ov->views[party].reconstructions->shares[i], buff + pos);
        }
        for (int i = 0; i < ARITH_MAXDEPTH; i++) {
            pos += serializeFieldElement(ov->views[party].oleOutputs->incMsgs0[i], buff + pos);
            pos += serializeFieldElement(ov->views[party].oleOutputs->incMsgs1[i], buff + pos);
        }

        if (party == 0) {
            SHA256_Update(&ctx1, buff, pos);
            SHA256_Final(expected1, &ctx1);
            memset(buff, 0, size);
        } else {
            SHA256_Update(&ctx2, buff, pos);
            SHA256_Final(expected2, &ctx2);
        }
    }

    if (memcmp(expected1, c1, SHA256_DIGEST_LENGTH) != 0 || memcmp(expected2, c2, SHA256_DIGEST_LENGTH) != 0) {
        fprintf(stderr, "commitment to views do not match");
    }

    free(buff);
}

int main(void) {
    initRand(globalRand);
    initField();
    proverTape = newRandomTape();



    // ********************  Prover PART ************************
    char *S = malloc(sizeof(char) * SECRET_LENGTH);
    for (int i = 0; i < SECRET_LENGTH; i++) {
        // The string is all zeroes because for evaluation, we use a secret that is all zero.
        // If we want to commit to a real secret, we need to replace this with the secret,
        // and then the XOR of the BBS PRG tail with this will be the commitment string.
        S[i] = '0';
    }

    mpz_t w;
    randGroupElement(w);

    SS input = secretShareNoTape(w);

    struct timeval start;
    gettimeofday(&start, NULL);
    Algorithm1Output out = runAlgorithm1(S, w, input);
    unsigned char sha256ProofShares[ITERATIONS][3][SECRET_LENGTH / 8];
    bitSharesTo8bitWords(out.booleanShares, sha256ProofShares);
    sha256TotalViews sha256Views = runSHA256Proof(sha256ProofShares);
    struct timeval end;
    gettimeofday(&end, NULL);
    time_t elapsedSeconds = end.tv_sec - start.tv_sec;
    printf("Proving elapsed: %d seconds\n", elapsedSeconds);

    // Commit to views
    unsigned char commitmentsOnViews[ITERATIONS][3][SHA256_DIGEST_LENGTH];
    // Commit to OLE outputs
    unsigned char commitmentsOnOLEs[ITERATIONS][3][SHA256_DIGEST_LENGTH];
#pragma omp parallel for default(none) shared(out, commitmentsOnViews, commitmentsOnOLEs)
    for (int i = 0; i < ITERATIONS; i++) {
        commitToViews(&out, i, commitmentsOnViews[i][0], commitmentsOnViews[i][1], commitmentsOnViews[i][2]);
        commitToOLEInputs(&out, i, commitmentsOnOLEs[i][0], commitmentsOnOLEs[i][1], commitmentsOnOLEs[i][2]);
    }


    // ********************  VERIFIER PART ************************
    int es[ITERATIONS];
    // Get random challenges from verifier.
    // We simulate an interactive verifier that asks to open random views.
    // If we wanted to use random oracle style NIZK then we would've needed to increase the rounds,
    // because our round count is only good enough for interactive zero knowledge.
    simulateVerifierViewSelection(es);
    Z *zs = packSHA256A(es, sha256Views);


    // Open the views according to what the verifier requested
    OpenedViewTuple *ov[ITERATIONS];
    openViews(ov, es, &out);
    int proofSize = countProofSize(ov);
    printf("Total proof size ~ %d MB\n", proofSize  / (1024 * 1024));


#pragma omp parallel for default(none) shared(es, ov, commitmentsOnViews)
    for (int i = 0; i < ITERATIONS; i++) {
        int n0 = (es[i] + 1) % 3;
        int n1 = (es[i] + 2) % 3;
        int p0 = MIN(n0, n1);
        int p1 = MAX(n0, n1);
        verifyCommittedViewPair(ov[i], commitmentsOnViews[i][p0], commitmentsOnViews[i][p1]);
    }

    // Run verification
    gettimeofday(&start, NULL);


#pragma omp parallel for default(none) shared(sha256ProofShares, out, input, S, ov, sha256Views, es, zs, commitmentsOnOLEs)
    for (int i = 0; i < ITERATIONS; i++) {
        int skippedParty = es[i];
        int n0 = (es[i] + 1) % 3;
        int n1 = (es[i] + 2) % 3;
        int p0 = MIN(n0, n1);
        int p1 = MAX(n0, n1);
        ComputedToOLE toOLE0;
        ComputedToOLE toOLE2;
        BooleanShares boolShares = verifyAlgorithm1(input, S, &(ov[i]->views[0]), &(ov[i]->views[1]), &toOLE0, &toOLE2);
        verifyCommittedOLEInputs(toOLE0, toOLE2, commitmentsOnOLEs[i][p0], commitmentsOnOLEs[i][p1]);
        // Sanity test - ensure bool shares from proof match bool shares from verification.
        for (int j = 0; j < SECRET_LENGTH; j++) {
            for (int p = 0; p < 3; p++) {
                if (skippedParty == p) {
                    continue;
                }
                char verificationProofShare = (boolShares.shares[j] >> p) & 1;
                char proofShare = (out.booleanShares[i].shares[j] >> p) & 1;
                if (verificationProofShare != proofShare) {
                    printf("proof share is %d but verification share is %d\n", proofShare, verificationProofShare);
                    exit(2);
                }
            }
        }

        // Convert AVTC boolean shares to SHA256 circuit representation
        unsigned char sha256VerificationShares[2][SECRET_LENGTH / 8];
        convertBoolShares(boolShares, sha256VerificationShares, skippedParty);
        compare8bitWords(i, skippedParty, sha256VerificationShares, sha256ProofShares);
        // Run the SHA256 verifier, pass it the boolean shares from the AVTC circuit
        // and override the input of the parties to force the SHA256 act on the AVTC circuit's output
        // Overwrite inputs with the output of the AVTC circuit
        int firstParty = 0;
        if (skippedParty == 1) {
            firstParty = 1;
        }
        for (int j = 0; j < 16; j++) {
            zs[i].ve.x[j] = sha256VerificationShares[firstParty][j];
            zs[i].ve1.x[j] = sha256VerificationShares[1 - firstParty][j];
        }
        int success = verifySHA256(sha256Views.as[i], (es[i] + 1) % 3, zs[i], sha256VerificationShares);
        if (success != 0) {
            exit(2);
        }
    }

    gettimeofday(&end, NULL);

    elapsedSeconds = end.tv_sec - start.tv_sec;
    printf("Verification elapsed: %d seconds\n", elapsedSeconds);


    return 0;
}
