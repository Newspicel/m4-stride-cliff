// Reproducer for a performance cliff on Apple M4: the same loop runs at different speeds depending
// on which loop ran before it, and much faster when its 8 store streams are not exactly 64 KiB apart.
//
// It runs five passes of a radix-8 Stockham FFT over 32768 complex f64 values (512 KiB per buffer),
// like RustFFT does, and prints the median time of every pass.
//
// The first pass comes in two versions that compute the same thing:
//   A: twiddles stored as two vectors {re, re} and {-im, im}, 32 bytes each
//   B: twiddles stored as plain complex numbers (16 bytes) and broadcast in registers
// The later passes are always the same code.
//
// Build: clang -O2 -o cliff cliff.c
// Usage: ./cliff <A|B> [--skew BYTES] [--dit] [--reps N]
//   --skew BYTES  move each of the 8 output streams of the later passes BYTES further apart
//   --dit         set the DIT (data independent timing) flag before running
#include <arm_neon.h>
#include <mach/mach_time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N 32768        // complex numbers per buffer
#define R 8            // radix of every pass
#define ROWS (N / R)   // 4096

typedef float64x2_t v2;

static long g_skew = 0;  // in doubles, added per output stream index j

static double now_us(void) {
    static mach_timebase_info_data_t tb;
    if (!tb.denom) mach_timebase_info(&tb);
    return (double)mach_absolute_time() * tb.numer / tb.denom / 1000.0;
}

static int cmp_double(const void *x, const void *y) {
    double a = *(const double *)x, b = *(const double *)y;
    return (a > b) - (a < b);
}

// x * -i for a forward FFT: (re, im) -> (im, -re)
static inline v2 rot90(v2 x, v2 sign) { return vmulq_f64(vextq_f64(x, x, 1), sign); }

// radix-8 butterfly
static inline void bf8(v2 *a, v2 sign, v2 h) {
    v2 u0 = vaddq_f64(a[0], a[4]), u1 = vaddq_f64(a[1], a[5]), u2 = vaddq_f64(a[2], a[6]), u3 = vaddq_f64(a[3], a[7]);
    v2 w0 = vsubq_f64(a[0], a[4]), w1 = vsubq_f64(a[1], a[5]), w2 = vsubq_f64(a[2], a[6]), w3 = vsubq_f64(a[3], a[7]);
    v2 e0 = vaddq_f64(u0, u2), e1 = vsubq_f64(u0, u2), e2 = vaddq_f64(u1, u3), e3 = rot90(vsubq_f64(u1, u3), sign);
    v2 E0 = vaddq_f64(e0, e2), E2 = vsubq_f64(e0, e2), E1 = vaddq_f64(e1, e3), E3 = vsubq_f64(e1, e3);
    v2 v1 = vmulq_f64(vaddq_f64(w1, rot90(w1, sign)), h);
    v2 v3 = vmulq_f64(vsubq_f64(rot90(w3, sign), w3), h);
    v2 t0 = vaddq_f64(w0, rot90(w2, sign)), t1 = vsubq_f64(w0, rot90(w2, sign));
    v2 t2 = vaddq_f64(v1, v3), t3 = rot90(vsubq_f64(v1, v3), sign);
    a[0] = E0; a[1] = vaddq_f64(t0, t2); a[2] = E1; a[3] = vaddq_f64(t1, t3);
    a[4] = E2; a[5] = vsubq_f64(t0, t2); a[6] = E3; a[7] = vsubq_f64(t1, t3);
}

// x * (re + i im) with re = (r, r) and im = (-i, i)
static inline v2 cmul_split(v2 x, v2 re, v2 im) {
    return vfmaq_f64(vmulq_f64(x, re), vextq_f64(x, x, 1), im);
}

// First pass, version A: split twiddle vectors loaded from memory.
__attribute__((noinline)) static void first_a(const double *in, double *out, const double *tw_split, v2 sign, v2 h) {
    for (int row = 0; row < ROWS; row++) {
        v2 a[8];
        for (int k = 0; k < 8; k++) a[k] = vld1q_f64(in + 2 * (row + k * ROWS));
        bf8(a, sign, h);
        const double *t = tw_split + row * 7 * 4;
        for (int j = 1; j < 8; j++) a[j] = cmul_split(a[j], vld1q_f64(t + 4 * (j - 1)), vld1q_f64(t + 4 * (j - 1) + 2));
        for (int j = 0; j < 8; j++) vst1q_f64(out + 2 * (8 * row + j), a[j]);
    }
}

// First pass, version B: plain complex twiddles, broadcast in registers.
__attribute__((noinline)) static void first_b(const double *in, double *out, const double *tw, v2 sign, v2 h) {
    const v2 neg_re = {-1.0, 1.0};
    for (int row = 0; row < ROWS; row++) {
        v2 a[8];
        for (int k = 0; k < 8; k++) a[k] = vld1q_f64(in + 2 * (row + k * ROWS));
        bf8(a, sign, h);
        const double *t = tw + row * 7 * 2;
        for (int j = 1; j < 8; j++) {
            v2 re = vld1q_dup_f64(t + 2 * (j - 1));
            v2 im = vmulq_f64(vld1q_dup_f64(t + 2 * (j - 1) + 1), neg_re);
            a[j] = cmul_split(a[j], re, im);
        }
        for (int j = 0; j < 8; j++) vst1q_f64(out + 2 * (8 * row + j), a[j]);
    }
}

// Later pass: `rows` groups of `stride` butterflies. Reads 8 streams `stride * rows` apart, writes 8
// streams `stride` apart (plus the optional skew). With rows == 1 (the last pass) both are 64 KiB.
__attribute__((noinline)) static void later(const double *in, double *out, int rows, int stride, const double *tw, v2 sign, v2 h) {
    const v2 neg_re = {-1.0, 1.0};
    for (int row = 0; row < rows; row++) {
        v2 re[8], im[8];
        if (rows > 1)
            for (int j = 1; j < 8; j++) {
                re[j] = vld1q_dup_f64(tw + 2 * (row * 7 + j - 1));
                im[j] = vmulq_f64(vld1q_dup_f64(tw + 2 * (row * 7 + j - 1) + 1), neg_re);
            }
        const double *ip = in + 2 * (stride * row);
        double *op = out + 2 * (stride * R * row);
        for (int q = 0; q < stride; q++) {
            v2 a[8];
            for (int k = 0; k < 8; k++) a[k] = vld1q_f64(ip + 2 * (q + stride * rows * k));
            bf8(a, sign, h);
            vst1q_f64(op + 2 * q, a[0]);
            for (int j = 1; j < 8; j++) {
                v2 v = rows > 1 ? cmul_split(a[j], re[j], im[j]) : a[j];
                vst1q_f64(op + 2 * (q + stride * j) + j * g_skew, v);
            }
        }
    }
}

int main(int argc, char **argv) {
    if (argc < 2 || (argv[1][0] != 'A' && argv[1][0] != 'B')) {
        fprintf(stderr, "usage: %s <A|B> [--skew BYTES] [--dit] [--reps N]\n", argv[0]);
        return 1;
    }
    char variant = argv[1][0];
    int reps = 3000, dit = 0;
    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "--skew") && i + 1 < argc) g_skew = atol(argv[++i]) / 8;
        else if (!strcmp(argv[i], "--dit")) dit = 1;
        else if (!strcmp(argv[i], "--reps") && i + 1 < argc) reps = atoi(argv[++i]);
    }
    if (dit) __asm__ volatile("msr DIT, #1");

    size_t bytes = (size_t)N * 16, slack = 16384 + 65536;
    double *src = aligned_alloc(16384, bytes + slack);
    double *buf = aligned_alloc(16384, bytes + slack);  // input of pass 0, output of the last pass
    double *s0 = (double *)((char *)aligned_alloc(16384, bytes + slack) + 1024);
    double *s1 = (double *)((char *)aligned_alloc(16384, bytes + slack) + 2048);
    double *tw_split = aligned_alloc(64, (size_t)ROWS * 7 * 32);
    double *tw = aligned_alloc(64, (size_t)ROWS * 7 * 16);
    for (int i = 0; i < 2 * N; i++) src[i] = (i * 0.37) - (int)(i * 0.37);
    for (size_t i = 0; i < (size_t)ROWS * 7 * 4; i++) tw_split[i] = 0.5 + (i % 7) * 0.01;
    for (size_t i = 0; i < (size_t)ROWS * 7 * 2; i++) tw[i] = 0.5 + (i % 5) * 0.01;
    v2 sign = {1.0, -1.0}, h = vdupq_n_f64(0.70710678118654752);

    // passes 1 to 4: (rows, stride)
    int prow[4] = {512, 64, 8, 1}, pstr[4] = {8, 64, 512, 4096};
    double *t = malloc(sizeof(double) * 5 * reps);
    for (int r = 0; r < reps; r++) {
        memcpy(buf, src, bytes);
        double t0 = now_us();
        if (variant == 'A') first_a(buf, s0, tw_split, sign, h);
        else first_b(buf, s0, tw, sign, h);
        t[r] = now_us() - t0;
        double *in = s0, *out = s1;
        for (int p = 0; p < 4; p++) {
            double *dst = p == 3 ? buf : out;
            double u0 = now_us();
            later(in, dst, prow[p], pstr[p], tw, sign, h);
            t[(p + 1) * reps + r] = now_us() - u0;
            double *tmp = in; in = out; out = tmp;
        }
    }
    printf("first pass %c, skew %4ld B, DIT %d   median us per pass:", variant, g_skew * 8, dit);
    for (int i = 0; i < 5; i++) {
        qsort(t + i * reps, reps, sizeof(double), cmp_double);
        printf(" %6.2f", t[i * reps + reps / 2]);
    }
    printf("\n");
    return 0;
}
