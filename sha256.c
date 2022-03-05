/*
 ============================================================================
 Name        : MPC_SHA256.c
 Author      : Sobuno
 Version     : 0.1
 Description : MPC SHA256 for one block only
 ============================================================================

 Modified to fit runSHA256Proof
 */


#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "shared.h"
#include "omp.h"
#include "avtc.h"


#define CH(e, f, g) ((e & f) ^ ((~e) & g))

void CMT(unsigned char k[16], V v, unsigned char r[4], unsigned char hash[SHA256_DIGEST_LENGTH]) {
    SHA256_CTX ctx;
    SHA256_Init(&ctx);
    SHA256_Update(&ctx, k, 16);
    SHA256_Update(&ctx, &v, sizeof(v));
    SHA256_Update(&ctx, r, 4);
    SHA256_Final(hash, &ctx);
}

int totalRandom = 0;
int totalSha = 0;
int totalSS = 0;
int totalHash = 0;


uint32_t rand32() {
    uint32_t x;
    x = rand() & 0xff;
    x |= (rand() & 0xff) << 8;
    x |= (rand() & 0xff) << 16;
    x |= (rand() & 0xff) << 24;

    return x;
}

void printbits(uint32_t n) {
    if (n) {
        printbits(n >> 1);
        printf("%d", n & 1);
    }

}


void mpc_XOR(uint32_t x[3], uint32_t y[3], uint32_t z[3]) {
    z[0] = x[0] ^ y[0];
    z[1] = x[1] ^ y[1];
    z[2] = x[2] ^ y[2];
}


void mpc_AND(uint32_t x[3], uint32_t y[3], uint32_t z[3], unsigned char *randomness[3], int *randCount, View views[3],
             int *countY) {
    uint32_t r[3] = {getRandom32(randomness[0], *randCount), getRandom32(randomness[1], *randCount),
                     getRandom32(randomness[2], *randCount)};
    *randCount += 4;
    uint32_t t[3] = {0};

    t[0] = (x[0] & y[1]) ^ (x[1] & y[0]) ^ (x[0] & y[0]) ^ r[0] ^ r[1];
    t[1] = (x[1] & y[2]) ^ (x[2] & y[1]) ^ (x[1] & y[1]) ^ r[1] ^ r[2];
    t[2] = (x[2] & y[0]) ^ (x[0] & y[2]) ^ (x[2] & y[2]) ^ r[2] ^ r[0];
    z[0] = t[0];
    z[1] = t[1];
    z[2] = t[2];
    views[0].y[*countY] = z[0];
    views[1].y[*countY] = z[1];
    views[2].y[*countY] = z[2];
    (*countY)++;
}


void mpc_NEGATE(uint32_t x[3], uint32_t z[3]) {
    z[0] = ~x[0];
    z[1] = ~x[1];
    z[2] = ~x[2];
}


void mpc_ADD(uint32_t x[3], uint32_t y[3], uint32_t z[3], unsigned char *randomness[3], int *randCount, View views[3],
             int *countY) {
    uint32_t c[3] = {0};
    uint32_t r[3] = {getRandom32(randomness[0], *randCount), getRandom32(randomness[1], *randCount),
                     getRandom32(randomness[2], *randCount)};
    *randCount += 4;

    uint8_t a[3], b[3];

    uint8_t t;

    for (int i = 0; i < 31; i++) {
        a[0] = GETBIT(x[0] ^ c[0], i);
        a[1] = GETBIT(x[1] ^ c[1], i);
        a[2] = GETBIT(x[2] ^ c[2], i);

        b[0] = GETBIT(y[0] ^ c[0], i);
        b[1] = GETBIT(y[1] ^ c[1], i);
        b[2] = GETBIT(y[2] ^ c[2], i);

        t = (a[0] & b[1]) ^ (a[1] & b[0]) ^ GETBIT(r[1], i);
        SETBIT(c[0], i + 1, t ^ (a[0] & b[0]) ^ GETBIT(c[0], i) ^ GETBIT(r[0], i));

        t = (a[1] & b[2]) ^ (a[2] & b[1]) ^ GETBIT(r[2], i);
        SETBIT(c[1], i + 1, t ^ (a[1] & b[1]) ^ GETBIT(c[1], i) ^ GETBIT(r[1], i));

        t = (a[2] & b[0]) ^ (a[0] & b[2]) ^ GETBIT(r[0], i);
        SETBIT(c[2], i + 1, t ^ (a[2] & b[2]) ^ GETBIT(c[2], i) ^ GETBIT(r[2], i));


    }

    z[0] = x[0] ^ y[0] ^ c[0];
    z[1] = x[1] ^ y[1] ^ c[1];
    z[2] = x[2] ^ y[2] ^ c[2];


    views[0].y[*countY] = c[0];
    views[1].y[*countY] = c[1];
    views[2].y[*countY] = c[2];
    *countY += 1;


}


void mpc_ADDK(uint32_t x[3], uint32_t y, uint32_t z[3], unsigned char *randomness[3], int *randCount, View views[3],
              int *countY) {
    uint32_t c[3] = {0};
    uint32_t r[3] = {getRandom32(randomness[0], *randCount), getRandom32(randomness[1], *randCount),
                     getRandom32(randomness[2], *randCount)};
    *randCount += 4;

    uint8_t a[3], b[3];

    uint8_t t;

    for (int i = 0; i < 31; i++) {
        a[0] = GETBIT(x[0] ^ c[0], i);
        a[1] = GETBIT(x[1] ^ c[1], i);
        a[2] = GETBIT(x[2] ^ c[2], i);

        b[0] = GETBIT(y ^ c[0], i);
        b[1] = GETBIT(y ^ c[1], i);
        b[2] = GETBIT(y ^ c[2], i);

        t = (a[0] & b[1]) ^ (a[1] & b[0]) ^ GETBIT(r[1], i);
        SETBIT(c[0], i + 1, t ^ (a[0] & b[0]) ^ GETBIT(c[0], i) ^ GETBIT(r[0], i));

        t = (a[1] & b[2]) ^ (a[2] & b[1]) ^ GETBIT(r[2], i);
        SETBIT(c[1], i + 1, t ^ (a[1] & b[1]) ^ GETBIT(c[1], i) ^ GETBIT(r[1], i));

        t = (a[2] & b[0]) ^ (a[0] & b[2]) ^ GETBIT(r[0], i);
        SETBIT(c[2], i + 1, t ^ (a[2] & b[2]) ^ GETBIT(c[2], i) ^ GETBIT(r[2], i));


    }

    z[0] = x[0] ^ y ^ c[0];
    z[1] = x[1] ^ y ^ c[1];
    z[2] = x[2] ^ y ^ c[2];


    views[0].y[*countY] = c[0];
    views[1].y[*countY] = c[1];
    views[2].y[*countY] = c[2];
    *countY += 1;

}


int sha256(unsigned char *result, unsigned char *input, int numBits) {
    uint32_t hA[8] = {0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
                      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19};


    if (numBits > 447) {
        printf("Input too long, aborting!");
        return -1;
    }
    int chars = numBits >> 3;
    unsigned char *chunk = calloc(64, 1); //512 bits
    memcpy(chunk, input, chars);
    chunk[chars] = 0x80;
    //Last 8 chars used for storing length of input without padding, in big-endian.
    //Since we only care for one block, we are safe with just using last 9 bits and 0'ing the rest

    //chunk[60] = numBits >> 24;
    //chunk[61] = numBits >> 16;
    chunk[62] = numBits >> 8;
    chunk[63] = numBits;

    uint32_t w[64];
    int i;
    for (i = 0; i < 16; i++) {
        w[i] = (chunk[i * 4] << 24) | (chunk[i * 4 + 1] << 16)
               | (chunk[i * 4 + 2] << 8) | chunk[i * 4 + 3];
    }

    uint32_t s0, s1;
    for (i = 16; i < 64; i++) {
        s0 = RIGHTROTATE(w[i - 15], 7) ^ RIGHTROTATE(w[i - 15], 18)
             ^ (w[i - 15] >> 3);
        s1 = RIGHTROTATE(w[i - 2], 17) ^ RIGHTROTATE(w[i - 2], 19)
             ^ (w[i - 2] >> 10);
        w[i] = w[i - 16] + s0 + w[i - 7] + s1;
    }

    uint32_t a, b, c, d, e, f, g, h, temp1, temp2, maj;
    a = hA[0];
    b = hA[1];
    c = hA[2];
    d = hA[3];
    e = hA[4];
    f = hA[5];
    g = hA[6];
    h = hA[7];

    for (i = 0; i < 64; i++) {
        s1 = RIGHTROTATE(e, 6) ^ RIGHTROTATE(e, 11) ^ RIGHTROTATE(e, 25);

        temp1 = h + s1 + CH(e, f, g) + k[i] + w[i];
        s0 = RIGHTROTATE(a, 2) ^ RIGHTROTATE(a, 13) ^ RIGHTROTATE(a, 22);


        maj = (a & (b ^ c)) ^ (b & c);
        temp2 = s0 + maj;


        h = g;
        g = f;
        f = e;
        e = d + temp1;
        d = c;
        c = b;
        b = a;
        a = temp1 + temp2;

    }

    hA[0] += a;
    hA[1] += b;
    hA[2] += c;
    hA[3] += d;
    hA[4] += e;
    hA[5] += f;
    hA[6] += g;
    hA[7] += h;

    for (i = 0; i < 8; i++) {
        result[i * 4] = (hA[i] >> 24);
        result[i * 4 + 1] = (hA[i] >> 16);
        result[i * 4 + 2] = (hA[i] >> 8);
        result[i * 4 + 3] = hA[i];
    }
    return 0;
}

void mpc_RIGHTROTATE(uint32_t x[], int i, uint32_t z[]) {
    z[0] = RIGHTROTATE(x[0], i);
    z[1] = RIGHTROTATE(x[1], i);
    z[2] = RIGHTROTATE(x[2], i);
}


void mpc_RIGHTSHIFT(uint32_t x[3], int i, uint32_t z[3]) {
    z[0] = x[0] >> i;
    z[1] = x[1] >> i;
    z[2] = x[2] >> i;
}


void mpc_MAJ(uint32_t a[], uint32_t b[3], uint32_t c[3], uint32_t z[3], unsigned char *randomness[3], int *randCount,
             View views[3], int *countY) {
    uint32_t t0[3];
    uint32_t t1[3];

    mpc_XOR(a, b, t0);
    mpc_XOR(a, c, t1);
    mpc_AND(t0, t1, z, randomness, randCount, views, countY);
    mpc_XOR(z, a, z);
}


void mpc_CH(uint32_t e[], uint32_t f[3], uint32_t g[3], uint32_t z[3], unsigned char *randomness[3], int *randCount,
            View views[3], int *countY) {
    uint32_t t0[3];

    //e & (f^g) ^ g
    mpc_XOR(f, g, t0);
    mpc_AND(e, t0, t0, randomness, randCount, views, countY);
    mpc_XOR(t0, g, z);

}


int mpc_sha256(unsigned char *results[3], unsigned char *inputs[3], int numBits, unsigned char *randomness[3],
               View views[3], int *countY) {


    if (numBits > 447) {
        printf("Input too long, aborting!");
        return -1;
    }

    int *randCount = calloc(1, sizeof(int));

    int chars = numBits >> 3;
    unsigned char *chunks[3];
    uint32_t w[64][3];

    for (int i = 0; i < 3; i++) {
        chunks[i] = calloc(64, 1); //512 bits
        memcpy(chunks[i], inputs[i], chars);
        chunks[i][chars] = 0x80;
        //Last 8 chars used for storing length of input without padding, in big-endian.
        //Since we only care for one block, we are safe with just using last 9 bits and 0'ing the rest

        //chunk[60] = numBits >> 24;
        //chunk[61] = numBits >> 16;
        chunks[i][62] = numBits >> 8;
        chunks[i][63] = numBits;
        memcpy(views[i].x, chunks[i], 64);

        for (int j = 0; j < 16; j++) {
            w[j][i] = (chunks[i][j * 4] << 24) | (chunks[i][j * 4 + 1] << 16)
                      | (chunks[i][j * 4 + 2] << 8) | chunks[i][j * 4 + 3];
        }
        free(chunks[i]);
    }


    uint32_t s0[3], s1[3];
    uint32_t t0[3], t1[3];
    for (int j = 16; j < 64; j++) {
        //s0[i] = RIGHTROTATE(w[i][j-15],7) ^ RIGHTROTATE(w[i][j-15],18) ^ (w[i][j-15] >> 3);
        mpc_RIGHTROTATE(w[j - 15], 7, t0);

        mpc_RIGHTROTATE(w[j - 15], 18, t1);
        mpc_XOR(t0, t1, t0);
        mpc_RIGHTSHIFT(w[j - 15], 3, t1);
        mpc_XOR(t0, t1, s0);

        //s1[i] = RIGHTROTATE(w[i][j-2],17) ^ RIGHTROTATE(w[i][j-2],19) ^ (w[i][j-2] >> 10);
        mpc_RIGHTROTATE(w[j - 2], 17, t0);
        mpc_RIGHTROTATE(w[j - 2], 19, t1);

        mpc_XOR(t0, t1, t0);
        mpc_RIGHTSHIFT(w[j - 2], 10, t1);
        mpc_XOR(t0, t1, s1);

        //w[i][j] = w[i][j-16]+s0[i]+w[i][j-7]+s1[i];

        mpc_ADD(w[j - 16], s0, t1, randomness, randCount, views, countY);
        mpc_ADD(w[j - 7], t1, t1, randomness, randCount, views, countY);
        mpc_ADD(t1, s1, w[j], randomness, randCount, views, countY);

    }

    uint32_t a[3] = {hA[0], hA[0], hA[0]};
    uint32_t b[3] = {hA[1], hA[1], hA[1]};
    uint32_t c[3] = {hA[2], hA[2], hA[2]};
    uint32_t d[3] = {hA[3], hA[3], hA[3]};
    uint32_t e[3] = {hA[4], hA[4], hA[4]};
    uint32_t f[3] = {hA[5], hA[5], hA[5]};
    uint32_t g[3] = {hA[6], hA[6], hA[6]};
    uint32_t h[3] = {hA[7], hA[7], hA[7]};
    uint32_t temp1[3], temp2[3], maj[3];
    for (int i = 0; i < 64; i++) {
        //s1 = RIGHTROTATE(e,6) ^ RIGHTROTATE(e,11) ^ RIGHTROTATE(e,25);
        mpc_RIGHTROTATE(e, 6, t0);
        mpc_RIGHTROTATE(e, 11, t1);
        mpc_XOR(t0, t1, t0);

        mpc_RIGHTROTATE(e, 25, t1);
        mpc_XOR(t0, t1, s1);


        //ch = (e & f) ^ ((~e) & g);
        //temp1 = h + s1 + CH(e,f,g) + k[i]+w[i];

        //t0 = h + s1

        mpc_ADD(h, s1, t0, randomness, randCount, views, countY);


        mpc_CH(e, f, g, t1, randomness, randCount, views, countY);

        //t1 = t0 + t1 (h+s1+ch)
        mpc_ADD(t0, t1, t1, randomness, randCount, views, countY);

        mpc_ADDK(t1, k[i], t1, randomness, randCount, views, countY);

        mpc_ADD(t1, w[i], temp1, randomness, randCount, views, countY);

        //s0 = RIGHTROTATE(a,2) ^ RIGHTROTATE(a,13) ^ RIGHTROTATE(a,22);
        mpc_RIGHTROTATE(a, 2, t0);
        mpc_RIGHTROTATE(a, 13, t1);
        mpc_XOR(t0, t1, t0);
        mpc_RIGHTROTATE(a, 22, t1);
        mpc_XOR(t0, t1, s0);


        mpc_MAJ(a, b, c, maj, randomness, randCount, views, countY);

        //temp2 = s0+maj;
        mpc_ADD(s0, maj, temp2, randomness, randCount, views, countY);

        memcpy(h, g, sizeof(uint32_t) * 3);
        memcpy(g, f, sizeof(uint32_t) * 3);
        memcpy(f, e, sizeof(uint32_t) * 3);
        //e = d+temp1;
        mpc_ADD(d, temp1, e, randomness, randCount, views, countY);
        memcpy(d, c, sizeof(uint32_t) * 3);
        memcpy(c, b, sizeof(uint32_t) * 3);
        memcpy(b, a, sizeof(uint32_t) * 3);
        //a = temp1+temp2;

        mpc_ADD(temp1, temp2, a, randomness, randCount, views, countY);
    }

    uint32_t hHa[8][3] = {{hA[0], hA[0], hA[0]},
                          {hA[1], hA[1], hA[1]},
                          {hA[2], hA[2], hA[2]},
                          {hA[3], hA[3], hA[3]},
                          {hA[4], hA[4], hA[4]},
                          {hA[5], hA[5], hA[5]},
                          {hA[6], hA[6], hA[6]},
                          {hA[7], hA[7], hA[7]}};
    mpc_ADD(hHa[0], a, hHa[0], randomness, randCount, views, countY);
    mpc_ADD(hHa[1], b, hHa[1], randomness, randCount, views, countY);
    mpc_ADD(hHa[2], c, hHa[2], randomness, randCount, views, countY);
    mpc_ADD(hHa[3], d, hHa[3], randomness, randCount, views, countY);
    mpc_ADD(hHa[4], e, hHa[4], randomness, randCount, views, countY);
    mpc_ADD(hHa[5], f, hHa[5], randomness, randCount, views, countY);
    mpc_ADD(hHa[6], g, hHa[6], randomness, randCount, views, countY);
    mpc_ADD(hHa[7], h, hHa[7], randomness, randCount, views, countY);

    for (int i = 0; i < 8; i++) {
        mpc_RIGHTSHIFT(hHa[i], 24, t0);
        results[0][i * 4] = t0[0];
        results[1][i * 4] = t0[1];
        results[2][i * 4] = t0[2];
        mpc_RIGHTSHIFT(hHa[i], 16, t0);
        results[0][i * 4 + 1] = t0[0];
        results[1][i * 4 + 1] = t0[1];
        results[2][i * 4 + 1] = t0[2];
        mpc_RIGHTSHIFT(hHa[i], 8, t0);
        results[0][i * 4 + 2] = t0[0];
        results[1][i * 4 + 2] = t0[1];
        results[2][i * 4 + 2] = t0[2];

        results[0][i * 4 + 3] = hHa[i][0];
        results[1][i * 4 + 3] = hHa[i][1];
        results[2][i * 4 + 3] = hHa[i][2];
    }
    free(randCount);

    return 0;
}


int writeToFile(char filename[], void *data, int size, int numItems) {
    FILE *file;

    file = fopen(filename, "wb");
    if (!file) {
        printf("Unable to open file!");
        return 1;
    }
    fwrite(data, size, numItems, file);
    fclose(file);
    return 0;
}


a commit(int numBytes, unsigned char shares[3][numBytes], unsigned char *randomness[3], View views[3]) {
    unsigned char *inputs[3];
    inputs[0] = shares[0];
    inputs[1] = shares[1];
    inputs[2] = shares[2];
    unsigned char *hashes[3];
    hashes[0] = malloc(32);
    hashes[1] = malloc(32);
    hashes[2] = malloc(32);

    int *countY = calloc(1, sizeof(int));
    mpc_sha256(hashes, inputs, numBytes * 8, randomness, views, countY);


    //Explicitly add y to view
    for (int i = 0; i < 8; i++) {
        views[0].y[*countY] = (hashes[0][i * 4] << 24) | (hashes[0][i * 4 + 1] << 16)
                              | (hashes[0][i * 4 + 2] << 8) | hashes[0][i * 4 + 3];

        views[1].y[*countY] = (hashes[1][i * 4] << 24) | (hashes[1][i * 4 + 1] << 16)
                              | (hashes[1][i * 4 + 2] << 8) | hashes[1][i * 4 + 3];
        views[2].y[*countY] = (hashes[2][i * 4] << 24) | (hashes[2][i * 4 + 1] << 16)
                              | (hashes[2][i * 4 + 2] << 8) | hashes[2][i * 4 + 3];
        *countY += 1;
    }
    free(countY);
    free(hashes[0]);
    free(hashes[1]);
    free(hashes[2]);

    uint32_t *result1 = malloc(32);
    output(views[0], result1);
    uint32_t *result2 = malloc(32);
    output(views[1], result2);
    uint32_t *result3 = malloc(32);
    output(views[2], result3);

    a a;
    memcpy(a.yp[0], result1, 32);
    memcpy(a.yp[1], result2, 32);
    memcpy(a.yp[2], result3, 32);

    free(result1);
    free(result2);
    free(result3);

    return a;
}

Z prove(int e, unsigned char keys[3][16], unsigned char rs[3][4], V views[3]) {
    Z z;
    memcpy(z.ke, keys[e], 16);
    memcpy(z.ke1, keys[(e + 1) % 3], 16);
    z.ve = views[e];
    z.ve1 = views[(e + 1) % 3];
    memcpy(z.re, rs[e], 4);
    memcpy(z.re1, rs[(e + 1) % 3], 4);

    return z;
}

Z *packSHA256A(int *es, sha256TotalViews sha256Views) {
    Z *zs = malloc(sizeof(Z) * ITERATIONS);

    for (int i = 0; i < ITERATIONS; i++) {
        zs[i] = prove((es[i] + 1) % 3, sha256Views.keys[i], sha256Views.rs[i], sha256Views.localViews[i]);
    }

    return zs;
}

void out(V v, uint32_t *result) {
    memcpy(result, &v.y[ySize - 8], 32);
}

z zFromZ(Z zz) {
    z z;

    memcpy(z.ve.x, zz.ve.x, sizeof(char) * 64);
    memcpy(z.ve.y, zz.ve.y, sizeof(uint32_t) * ySize);

    memcpy(z.ve1.x, zz.ve1.x, sizeof(char) * 64);
    memcpy(z.ve1.y, zz.ve1.y, sizeof(uint32_t) * ySize);

    memcpy(z.re, zz.re, sizeof(char) * 4);
    memcpy(z.re1, zz.re1, sizeof(char) * 4);

    memcpy(z.ke, zz.ke, sizeof(char) * 16);
    memcpy(z.ke1, zz.ke1, sizeof(char) * 16);

    return z;
}

sha256TotalViews runSHA256Proof(unsigned char shares[ITERATIONS][3][SECRET_LENGTH / 8]) {
    setbuf(stdout, NULL);
    srand((unsigned) time(NULL));
    init_EVP();
    openmp_thread_setup();

    unsigned char garbage[4];
    if (RAND_bytes(garbage, 4) != 1) {
        printf("RAND_bytes failed crypto, aborting\n");
        exit(2);
    }


    sha256TotalViews sha256Views;

    //Generating keys
    if (RAND_bytes(sha256Views.keys, ITERATIONS * 3 * 16) != 1) {
        printf("RAND_bytes failed crypto, aborting\n");
        exit(2);
    }
    if (RAND_bytes(sha256Views.keys, ITERATIONS * 3 * 4) != 1) {
        printf("RAND_bytes failed crypto, aborting\n");
        exit(2);
    }



    //Generating randomness
    unsigned char *randomness[ITERATIONS][3];
#pragma omp parallel for default(none) shared(randomness, sha256Views)
    for (int k = 0; k < ITERATIONS; k++) {
        for (int j = 0; j < 3; j++) {
            randomness[k][j] = malloc(2912 * sizeof(unsigned char));
            getAllRandomness(sha256Views.keys[k][j], randomness[k][j]);
        }
    }

    //Running MPC-SHA2
#pragma omp parallel for default(none) shared(sha256Views, shares, randomness)
    for (int k = 0; k < ITERATIONS; k++) {
        a iterationViews = commit(SECRET_LENGTH / 8, shares[k], randomness[k], sha256Views.localViews[k]);
        memcpy(sha256Views.as[k].yp, &iterationViews.yp, sizeof(int) * 3 * 8);
        for (int j = 0; j < 3; j++) {
            free(randomness[k][j]);
        }
    }

    for (int k = 0; k < ITERATIONS; k++) {
        unsigned char hash1[SHA256_DIGEST_LENGTH];
        CMT(sha256Views.keys[k][0], sha256Views.localViews[k][0], sha256Views.rs[k][0], &hash1);
        memcpy(sha256Views.as[k].h[0], &hash1, 32);
        CMT(sha256Views.keys[k][1], sha256Views.localViews[k][1], sha256Views.rs[k][1], &hash1);
        memcpy(sha256Views.as[k].h[1], &hash1, 32);
        CMT(sha256Views.keys[k][2], sha256Views.localViews[k][2], sha256Views.rs[k][2], &hash1);
        memcpy(sha256Views.as[k].h[2], &hash1, 32);
    }

    return sha256Views;

}

int verifySHA256(A a, int e, Z zz, unsigned char inputs[2][SECRET_LENGTH / 8]) {
    z z = zFromZ(zz);
    unsigned char* hash = malloc(SHA256_DIGEST_LENGTH);
    H(z.ke, z.ve, z.re, hash);

    if (memcmp(a.h[e], hash, 32) != 0) {
        printf("Failing at %d", __LINE__);
        return 1;
    }
    H(z.ke1, z.ve1, z.re1, hash);
    if (memcmp(a.h[(e + 1) % 3], hash, 32) != 0) {
        printf("Failing at %d", __LINE__);
        return 1;
    }
    free(hash);

    uint32_t* result = malloc(32);
    output(z.ve, result);
    if (memcmp(a.yp[e], result, 32) != 0) {
        printf("Failing at %d", __LINE__);
        return 1;
    }

    output(z.ve1, result);
    if (memcmp(a.yp[(e + 1) % 3], result, 32) != 0) {
        printf("Failing at %d", __LINE__);
        return 1;
    }

    free(result);

    unsigned char randomness[2][2912];
    getAllRandomness(z.ke, randomness[0]);
    getAllRandomness(z.ke1, randomness[1]);

    int* randCount = calloc(1, sizeof(int));
    int* countY = calloc(1, sizeof(int));

    uint32_t w[64][2];
    for (int j = 0; j < 16; j++) {
        w[j][0] = (z.ve.x[j * 4] << 24) | (z.ve.x[j * 4 + 1] << 16)
                  | (z.ve.x[j * 4 + 2] << 8) | z.ve.x[j * 4 + 3];
        w[j][1] = (z.ve1.x[j * 4] << 24) | (z.ve1.x[j * 4 + 1] << 16)
                  | (z.ve1.x[j * 4 + 2] << 8) | z.ve1.x[j * 4 + 3];
    }

    uint32_t s0[2], s1[2];
    uint32_t t0[2], t1[2];
    for (int j = 16; j < 64; j++) {
        //s0[i] = RIGHTROTATE(w[i][j-15],7) ^ RIGHTROTATE(w[i][j-15],18) ^ (w[i][j-15] >> 3);
        mpc_RIGHTROTATE2(w[j-15], 7, t0);
        mpc_RIGHTROTATE2(w[j-15], 18, t1);
        mpc_XOR2(t0, t1, t0);
        mpc_RIGHTSHIFT2(w[j-15], 3, t1);
        mpc_XOR2(t0, t1, s0);

        //s1[i] = RIGHTROTATE(w[i][j-2],17) ^ RIGHTROTATE(w[i][j-2],19) ^ (w[i][j-2] >> 10);
        mpc_RIGHTROTATE2(w[j-2], 17, t0);
        mpc_RIGHTROTATE2(w[j-2], 19, t1);
        mpc_XOR2(t0, t1, t0);
        mpc_RIGHTSHIFT2(w[j-2],10,t1);
        mpc_XOR2(t0, t1, s1);

        //w[i][j] = w[i][j-16]+s0[i]+w[i][j-7]+s1[i];

        if(mpc_ADD_verify(w[j-16], s0, t1, z.ve, z.ve1, randomness, randCount, countY) == 1) {
            printf("Failing at %d, iteration %d", __LINE__, j);
            return 1;
        }


        if(mpc_ADD_verify(w[j-7], t1, t1, z.ve, z.ve1, randomness, randCount, countY) == 1) {
            printf("Failing at %d, iteration %d", __LINE__, j);
            return 1;
        }
        if(mpc_ADD_verify(t1, s1, w[j], z.ve, z.ve1, randomness, randCount, countY) == 1) {
            printf("Failing at %d, iteration %d", __LINE__, j);
            return 1;
        }

    }



    uint32_t va[2] = { hA[0],hA[0] };
    uint32_t vb[2] = { hA[1],hA[1] };
    uint32_t vc[2] = { hA[2],hA[2] };
    uint32_t vd[2] = { hA[3],hA[3] };
    uint32_t ve[2] = { hA[4],hA[4] };
    uint32_t vf[2] = { hA[5],hA[5] };
    uint32_t vg[2] = { hA[6],hA[6] };
    uint32_t vh[2] = { hA[7],hA[7] };
    uint32_t temp1[3], temp2[3], maj[3];
    for (int i = 0; i < 64; i++) {
        //s1 = RIGHTROTATE(e,6) ^ RIGHTROTATE(e,11) ^ RIGHTROTATE(e,25);
        mpc_RIGHTROTATE2(ve, 6, t0);
        mpc_RIGHTROTATE2(ve, 11, t1);
        mpc_XOR2(t0, t1, t0);
        mpc_RIGHTROTATE2(ve, 25, t1);
        mpc_XOR2(t0, t1, s1);




        //ch = (e & f) ^ ((~e) & g);
        //temp1 = h + s1 + CH(e,f,g) + k[i]+w[i];

        //t0 = h + s1

        if(mpc_ADD_verify(vh, s1, t0, z.ve, z.ve1, randomness, randCount, countY) == 1) {
            printf("Failing at %d, iteration %d", __LINE__, i);
            return 1;
        }



        if(mpc_CH_verify(ve, vf, vg, t1, z.ve, z.ve1, randomness, randCount, countY) == 1) {
            printf("Failing at %d, iteration %d", __LINE__, i);
            return 1;
        }

        //t1 = t0 + t1 (h+s1+ch)
        if(mpc_ADD_verify(t0, t1, t1, z.ve, z.ve1, randomness, randCount, countY) == 1) {
            printf("Failing at %d, iteration %d", __LINE__, i);
            return 1;
        }



        t0[0] = k[i];
        t0[1] = k[i];
        if(mpc_ADD_verify(t1, t0, t1, z.ve, z.ve1, randomness, randCount, countY) == 1) {
            printf("Failing at %d, iteration %d", __LINE__, i);
            return 1;
        }



        if(mpc_ADD_verify(t1, w[i], temp1, z.ve, z.ve1, randomness, randCount, countY) == 1) {
            printf("Failing at %d, iteration %d", __LINE__, i);
            return 1;
        }

        //s0 = RIGHTROTATE(a,2) ^ RIGHTROTATE(a,13) ^ RIGHTROTATE(a,22);
        mpc_RIGHTROTATE2(va, 2, t0);
        mpc_RIGHTROTATE2(va, 13, t1);
        mpc_XOR2(t0, t1, t0);
        mpc_RIGHTROTATE2(va, 22, t1);
        mpc_XOR2(t0, t1, s0);

        //maj = (a & (b ^ c)) ^ (b & c);
        //(a & b) ^ (a & c) ^ (b & c)

        if(mpc_MAJ_verify(va, vb, vc, maj, z.ve, z.ve1, randomness, randCount, countY) == 1) {
            printf("Failing at %d, iteration %d", __LINE__, i);
            return 1;
        }

        //temp2 = s0+maj;
        if(mpc_ADD_verify(s0, maj, temp2, z.ve, z.ve1, randomness, randCount, countY) == 1) {
            printf("Failing at %d, iteration %d", __LINE__, i);
            return 1;
        }



        memcpy(vh, vg, sizeof(uint32_t) * 2);
        memcpy(vg, vf, sizeof(uint32_t) * 2);
        memcpy(vf, ve, sizeof(uint32_t) * 2);
        //e = d+temp1;
        if(mpc_ADD_verify(vd, temp1, ve, z.ve, z.ve1, randomness, randCount, countY) == 1) {
            printf("Failing at %d, iteration %d", __LINE__, i);
            return 1;
        }

        memcpy(vd, vc, sizeof(uint32_t) * 2);
        memcpy(vc, vb, sizeof(uint32_t) * 2);
        memcpy(vb, va, sizeof(uint32_t) * 2);
        //a = temp1+temp2;

        if(mpc_ADD_verify(temp1, temp2, va, z.ve, z.ve1, randomness, randCount, countY) == 1) {
            printf("Failing at %d, iteration %d", __LINE__, i);
            return 1;
        }
    }

    uint32_t hHa[8][3] = { { hA[0],hA[0],hA[0]  }, { hA[1],hA[1],hA[1] }, { hA[2],hA[2],hA[2] }, { hA[3],hA[3],hA[3] },
                           { hA[4],hA[4],hA[4] }, { hA[5],hA[5],hA[5] }, { hA[6],hA[6],hA[6] }, { hA[7],hA[7],hA[7] } };
    if(mpc_ADD_verify(hHa[0], va, hHa[0], z.ve, z.ve1, randomness, randCount, countY) == 1) {
        printf("Failing at %d", __LINE__);
        return 1;
    }
    if(mpc_ADD_verify(hHa[1], vb, hHa[1], z.ve, z.ve1, randomness, randCount, countY) == 1) {
#if VERBOSE
        printf("Failing at %d", __LINE__);
#endif
        return 1;
    }
    if(mpc_ADD_verify(hHa[2], vc, hHa[2], z.ve, z.ve1, randomness, randCount, countY) == 1) {
#if VERBOSE
        printf("Failing at %d", __LINE__);
#endif
        return 1;
    }
    if(mpc_ADD_verify(hHa[3], vd, hHa[3], z.ve, z.ve1, randomness, randCount, countY) == 1) {
#if VERBOSE
        printf("Failing at %d", __LINE__);
#endif
        return 1;
    }
    if(mpc_ADD_verify(hHa[4], ve, hHa[4], z.ve, z.ve1, randomness, randCount, countY) == 1) {
#if VERBOSE
        printf("Failing at %d", __LINE__);
#endif
        return 1;
    }
    if(mpc_ADD_verify(hHa[5], vf, hHa[5], z.ve, z.ve1, randomness, randCount, countY) == 1) {
#if VERBOSE
        printf("Failing at %d", __LINE__);
#endif
        return 1;
    }
    if(mpc_ADD_verify(hHa[6], vg, hHa[6], z.ve, z.ve1, randomness, randCount, countY) == 1) {
#if VERBOSE
        printf("Failing at %d", __LINE__);
#endif
        return 1;
    }
    if(mpc_ADD_verify(hHa[7], vh, hHa[7], z.ve, z.ve1, randomness, randCount, countY) == 1) {
#if VERBOSE
        printf("Failing at %d", __LINE__);
#endif
        return 1;
    }

    free(randCount);
    free(countY);

    return 0;
}
