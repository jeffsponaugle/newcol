/*
 * collatzdig -- enumerate n such that the digits of n (with multiplicity)
 * are NOT all encountered in the Collatz trajectory n -> ... -> 1
 * (excluding the initial n itself).
 *
 *   5  is a member: 16 8 4 2 1 contains no 5.
 *   29 is a member: trajectory contains a 2 but never a 9.
 *   77 is a member: trajectory contains only one 7, two are needed.
 *   30 is not:      15 46 23 70 covers both the 3 and the 0.
 *
 * Build:  cc -O3 -march=native -pthread -o collatzdig collatzdig.c
 * Usage:  collatzdig [-t threads] [-b blocksize] [-p secs] [-q] start end
 *         collatzdig -x n            explain a single n
 *
 * Output: one term per line on stdout, in increasing order.
 * Progress / summary on stderr.
 */
#define _GNU_SOURCE
#include <errno.h>
#include <inttypes.h>
#include <pthread.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

typedef uint64_t u64;
typedef unsigned __int128 u128;

/* Switch to 128-bit arithmetic when an odd x exceeds this (test builds lower it). */
#ifndef OVERFLOW_LIMIT
#define OVERFLOW_LIMIT ((UINT64_MAX - 1) / 3)
#endif

/* ------------------------------------------------------------------ */
/* Packed digit counters.                                               */
/* A u64 holds ten 6-bit lanes, lane d (bits 6d..6d+5) = count of digit  */
/* d.  cnt4[c] is the packed digit count of the 4-digit chunk c written  */
/* with leading zeros.  The "need" state of an n is such a word; each    */
/* trajectory value's packed count is subtracted lane-wise, saturating   */
/* at zero, in a handful of branch-free ops.  need == 0 means covered.   */
/* Lane values never exceed 31, so bit 5 of each lane can act as a       */
/* borrow guard.                                                          */
/* ------------------------------------------------------------------ */
static u64 cnt4[10000];

#define LANE(v) ( (u64)(v) * 0x0041041041041041ULL )   /* v (<64) in every lane */
#define LANE_H  LANE(0x20)                              /* guard bits */
#define LANE_1F LANE(0x1F)

static void init_tables(void)
{
    for (int v = 0; v < 10000; v++) {
        int a = v; u64 p = 0;
        for (int i = 0; i < 4; i++) { p += (u64)1 << (6 * (a % 10)); a /= 10; }
        cnt4[v] = p;
    }
}

/* need := max(need - have, 0) lane-wise.  Requires every lane of have <= 31. */
static inline u64 lanes_sub(u64 need, u64 have)
{
    u64 d = (need | LANE_H) - have;          /* guard bit survives iff need >= have */
    return d & (((d & LANE_H) >> 5) * 0x1F); /* keep need-have where it survived, else 0 */
}

/* Packed digit count of a 64-bit value (no leading zeros). */
static inline u64 count64(u64 y)
{
    u64 have = 0;
    while (y >= 10000) { have += cnt4[y % 10000]; y /= 10000; }
    /* top chunk: remove the leading zeros that cnt4 counted */
    return have + cnt4[y] - ((y < 10) + (y < 100) + (y < 1000));
}

/* Consume the digits of a 64-bit trajectory value. Returns 1 when covered. */
static inline int consume64(u64 *need, u64 y)
{
    *need = lanes_sub(*need, count64(y));
    return *need == 0;
}

/* Consume the digits of a 128-bit value (rare: only for huge peaks).
 * Processed in 19-digit pieces so no lane of "have" exceeds 31. */
static int consume128(u64 *need, u128 y)
{
    const u64 P19 = 10000000000000000000ULL;      /* 10^19 */
    if (y < P19) return consume64(need, (u64)y);
    u64 lo = (u64)(y % P19), have = 0;
    for (int i = 0; i < 19; i++) { have += (u64)1 << (6 * (lo % 10)); lo /= 10; }
    *need = lanes_sub(*need, have);
    if (*need == 0) return 1;
    return consume128(need, y / P19);
}

/* 128-bit trajectory section. Returns 1 if covered, 0 if reached 1,
 * and 2 (with *px updated) once x fits comfortably in 64 bits again. */
static int run128(u64 *need, u128 *px)
{
    u128 x = *px;
    while (x != 1) {
        if (x & 1) x = 3 * x + 1; else x >>= 1;
        if (consume128(need, x)) return 1;
        if (x < ((u128)1 << 62)) { *px = x; return 2; }
    }
    return 0;
}

/* Returns 1 if n is a member of the sequence, 0 otherwise.
 * need = packed digit count of n (see count64). */
static int is_member_need(u64 n, u64 need)
{
    u64 x = n;
    while (x != 1) {
        if (x > OVERFLOW_LIMIT) {                     /* 3x+1 might overflow */
            u128 w = x;
            int r = run128(&need, &w);
            if (r != 2) return r == 0;
            x = (u64)w;
            continue;
        }
        x = (x & 1) ? 3 * x + 1 : x >> 1;             /* branch-free select */
        if (consume64(&need, x)) return 0;
    }
    return 1;
}

/* ------------------------------------------------------------------ */
/* Explain a single n (debug / verification).                           */
/* ------------------------------------------------------------------ */
static void print_need(u64 need)
{
    int any = 0;
    for (int d = 0; d < 10; d++) {
        unsigned c = (need >> (6 * d)) & 0x3F;
        if (c) { printf(" %d x%u", d, c); any = 1; }
    }
    if (!any) printf(" (none)");
}

static void explain(u64 n)
{
    u64 need = count64(n);
    printf("%" PRIu64 ": need", n); print_need(need); printf("\n");
    u128 x = n;
    int covered = 0;
    while (x != 1 && !covered) {
        if (x & 1) x = 3 * x + 1; else x >>= 1;
        covered = consume128(&need, x);
        char buf[64]; int p = 63; buf[p] = 0;
        u128 t = x; do { buf[--p] = '0' + (int)(t % 10); t /= 10; } while (t);
        printf("  %s   still need:", buf + p); print_need(need); printf("\n");
    }
    printf("%" PRIu64 " is %sa member\n", n, covered ? "NOT " : "");
}

/* ------------------------------------------------------------------ */
/* Threading with in-order output.                                      */
/* ------------------------------------------------------------------ */
struct slot {
    u64  *terms;
    size_t count;
    int   ready;
};

static u64 g_start, g_end, g_block, g_nblocks;
static _Atomic u64 g_next_block;
static _Atomic u64 g_found;
static _Atomic u64 g_printed_blocks;
static pthread_mutex_t g_mu = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t  g_cv = PTHREAD_COND_INITIALIZER;
static struct slot *g_ring;
static size_t g_ring_size;
static u64 g_printed;          /* next block index to print (under g_mu) */
static int g_quiet;
static int g_status_tty;       /* status line uses \r on a terminal */
static int g_clear_for_terms;  /* stdout shares the terminal: wipe status before printing terms */

static void *worker(void *arg)
{
    (void)arg;
    size_t cap = 64;
    u64 *buf = malloc(cap * sizeof *buf);
    if (!buf) { perror("malloc"); exit(1); }

    for (;;) {
        u64 b = atomic_fetch_add(&g_next_block, 1);
        if (b >= g_nblocks) break;

        u64 lo = g_start + b * g_block;
        u64 hi = lo + g_block; if (hi > g_end) hi = g_end;   /* [lo,hi) */
        size_t cnt = 0;
        /* digits of n = digits of q=n/10000 (constant for 10000 consecutive n)
         * plus the 4-digit chunk r=n%10000, so "need" costs one table lookup */
        u64 q = lo / 10000, r = lo % 10000;
        u64 hi_need = q ? count64(q) : 0;
        for (u64 n = lo; n < hi; n++) {
            u64 need = q ? hi_need + cnt4[r] : count64(n);   /* n<10000: no leading zeros */
            if (++r == 10000) { r = 0; q++; hi_need = count64(q); }
            if (is_member_need(n, need)) {
                if (cnt == cap) {
                    cap *= 2;
                    buf = realloc(buf, cap * sizeof *buf);
                    if (!buf) { perror("realloc"); exit(1); }
                }
                buf[cnt++] = n;
            }
        }

        /* hand the block's results to the ordered printer */
        u64 *out = NULL;
        if (cnt) {
            out = malloc(cnt * sizeof *out);
            if (!out) { perror("malloc"); exit(1); }
            memcpy(out, buf, cnt * sizeof *out);
        }
        pthread_mutex_lock(&g_mu);
        while (b - g_printed >= g_ring_size)
            pthread_cond_wait(&g_cv, &g_mu);
        struct slot *sl = &g_ring[b % g_ring_size];
        sl->terms = out; sl->count = cnt; sl->ready = 1;
        /* flush every consecutive finished block */
        for (;;) {
            struct slot *f = &g_ring[g_printed % g_ring_size];
            if (!f->ready) break;
            if (f->count && g_clear_for_terms) fputs("\r\033[K", stderr);
            for (size_t i = 0; i < f->count; i++)
                printf("%" PRIu64 "\n", f->terms[i]);
            if (f->count) { fflush(stdout); atomic_fetch_add(&g_found, f->count); }
            free(f->terms); f->terms = NULL; f->count = 0; f->ready = 0;
            g_printed++;
            atomic_store(&g_printed_blocks, g_printed);
        }
        pthread_cond_broadcast(&g_cv);
        pthread_mutex_unlock(&g_mu);
    }
    free(buf);
    return NULL;
}

static double now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* Parse an unsigned integer with optional decimal point, exponent and
 * scale suffix, exactly:  123  1e14  2.5e9  10M  1.5B  3T  0x1000
 * Suffixes (case-insensitive): K=1e3 M=1e6 B=G=1e9 T=1e12 P=1e15 */
/* 12345678901 -> "12.35B" */
static const char *human(u64 v, char *buf, size_t n)
{
    if      (v >= 1000000000000000ULL) snprintf(buf, n, "%.3fP", v / 1e15);
    else if (v >= 1000000000000ULL)    snprintf(buf, n, "%.3fT", v / 1e12);
    else if (v >= 1000000000ULL)       snprintf(buf, n, "%.3fB", v / 1e9);
    else if (v >= 1000000ULL)          snprintf(buf, n, "%.3fM", v / 1e6);
    else if (v >= 1000ULL)             snprintf(buf, n, "%.3fK", v / 1e3);
    else                               snprintf(buf, n, "%" PRIu64, v);
    return buf;
}

/* seconds -> "1:02:03" or "3d 04:05:06" */
static const char *duration(double secs, char *buf, size_t n)
{
    if (secs < 0 || secs != secs || secs > 1e9) { snprintf(buf, n, "--:--:--"); return buf; }
    u64 t = (u64)(secs + 0.5), d = t / 86400; t %= 86400;
    if (d) snprintf(buf, n, "%" PRIu64 "d %02" PRIu64 ":%02" PRIu64 ":%02" PRIu64, d, t / 3600, t % 3600 / 60, t % 60);
    else   snprintf(buf, n, "%" PRIu64 ":%02" PRIu64 ":%02" PRIu64, t / 3600, t % 3600 / 60, t % 60);
    return buf;
}

static u64 parse_u64(const char *s)
{
    const char *p = s;
    u128 mant = 0;
    int frac = 0, ndig = 0, exp10 = 0;

    if (s[0] == '0' && (s[1] == 'x' || s[1] == 'X')) {       /* hex */
        char *end; errno = 0;
        unsigned long long v = strtoull(s, &end, 16);
        if (errno || *end || end == s + 2) goto bad;
        return v;
    }
    for (; *p >= '0' && *p <= '9'; p++, ndig++) {
        mant = mant * 10 + (*p - '0');
        if (mant > ((u128)1 << 100)) goto bad;
    }
    if (*p == '.') {
        for (p++; *p >= '0' && *p <= '9'; p++, ndig++, frac++) {
            mant = mant * 10 + (*p - '0');
            if (mant > ((u128)1 << 100)) goto bad;
        }
    }
    if (ndig == 0) goto bad;
    if (*p == 'e' || *p == 'E') {
        char *end;
        long e = strtol(p + 1, &end, 10);
        if (end == p + 1 || e > 40 || e < -40) goto bad;
        exp10 += (int)e; p = end;
    }
    switch (*p) {
    case 'k': case 'K': exp10 += 3;  p++; break;
    case 'm': case 'M': exp10 += 6;  p++; break;
    case 'b': case 'B':
    case 'g': case 'G': exp10 += 9;  p++; break;
    case 't': case 'T': exp10 += 12; p++; break;
    case 'p': case 'P': exp10 += 15; p++; break;
    default: break;
    }
    if (*p) goto bad;
    exp10 -= frac;
    for (; exp10 < 0; exp10++) {                 /* e.g. 1.25K -> must stay integral */
        if (mant % 10) goto bad;
        mant /= 10;
    }
    for (; exp10 > 0; exp10--) {
        mant *= 10;
        if (mant > UINT64_MAX) goto bad;
    }
    return (u64)mant;
bad:
    fprintf(stderr, "bad number: %s (use e.g. 123, 1e14, 2.5e9, 10M, 1.5B, 3T)\n", s);
    exit(2);
}

static void usage(void)
{
    fprintf(stderr,
        "usage: collatzdig [-t threads] [-b blocksize] [-p progress_secs] [-q] start end\n"
        "       collatzdig -x n\n"
        "  Enumerates n in [start,end] whose digits (with multiplicity) are not all\n"
        "  encountered in the Collatz trajectory from n down to 1.\n"
        "  -t N   worker threads (default: number of CPUs)\n"
        "  -b N   numbers per work block (default 65536)\n"
        "  -p S   status line refresh interval in seconds on stderr (default 1, 0=off)\n"
        "  -q     no summary/progress on stderr\n"
        "  -x n   print the trajectory of n with the digits still needed\n"
        "  Numbers accept scientific notation and suffixes: 1e14, 2.5e9, 10M, 1.5B, 3T, 1P\n");
    exit(2);
}

int main(int argc, char **argv)
{
    long threads = sysconf(_SC_NPROCESSORS_ONLN);
    if (threads < 1) threads = 1;
    u64 block = 65536;
    double prog = 1.0;
    int opt;

    init_tables();

    while ((opt = getopt(argc, argv, "t:b:p:qx:h")) != -1) {
        switch (opt) {
        case 't': threads = (long)parse_u64(optarg); if (threads < 1) threads = 1; break;
        case 'b': block = parse_u64(optarg); if (block < 1) block = 1; break;
        case 'p': prog = atof(optarg); break;
        case 'q': g_quiet = 1; break;
        case 'x': explain(parse_u64(optarg)); return 0;
        default: usage();
        }
    }
    if (argc - optind != 2) usage();
    g_start = parse_u64(argv[optind]);
    g_end   = parse_u64(argv[optind + 1]);
    if (g_start < 1) g_start = 1;           /* 0 never reaches 1 */
    if (g_end < g_start) return 0;
    g_end += 1;                              /* make [start,end] -> [start,end) */
    g_block   = block;
    g_nblocks = (g_end - g_start + g_block - 1) / g_block;
    g_ring_size = (size_t)threads * 8;
    g_ring = calloc(g_ring_size, sizeof *g_ring);
    if (!g_ring) { perror("calloc"); return 1; }

    g_status_tty = !g_quiet && prog > 0 && isatty(2);
    g_clear_for_terms = g_status_tty && isatty(1);

    pthread_t *th = malloc(threads * sizeof *th);
    double t0 = now();
    for (long i = 0; i < threads; i++)
        if (pthread_create(&th[i], NULL, worker, NULL)) { perror("pthread_create"); return 1; }

    if (!g_quiet && prog > 0) {
        double last = t0; u64 last_n = g_start;
        const u64 total = g_end - g_start;
        for (;;) {
            usleep(100000);
            u64 done_blocks = atomic_load(&g_printed_blocks);
            int finished = done_blocks >= g_nblocks;
            double t = now();
            if (finished || t - last >= prog) {
                u64 upto = g_start + done_blocks * g_block;
                if (upto > g_end) upto = g_end;
                u64 done = upto - g_start;
                double elapsed = t - t0;
                double inst = (t > last) ? (double)(upto - last_n) / (t - last) : 0;
                double avg  = elapsed > 0 ? (double)done / elapsed : 0;
                double eta  = avg > 0 ? (double)(total - done) / avg : -1;
                char hb[32], eb[32], tb[32];
                fprintf(stderr, "%s%6.2f%%  at %-9s found %-6" PRIu64 " %7.1f M/s  elapsed %s  ETA %s%s",
                        g_status_tty ? "\r\033[K" : "",
                        100.0 * (double)done / (double)total,
                        human(upto, hb, sizeof hb), atomic_load(&g_found), inst / 1e6,
                        duration(elapsed, eb, sizeof eb), duration(eta, tb, sizeof tb),
                        g_status_tty ? "" : "\n");
                fflush(stderr);
                last = t; last_n = upto;
            }
            if (finished) break;
        }
        if (g_status_tty) fputs("\r\033[K", stderr);
    }
    for (long i = 0; i < threads; i++) pthread_join(th[i], NULL);
    double t1 = now();
    if (!g_quiet)
        fprintf(stderr, "done: [%" PRIu64 ", %" PRIu64 "]  %" PRIu64 " terms  %.2fs  (%.1f M/s, %ld threads)\n",
                g_start, g_end - 1, atomic_load(&g_found), t1 - t0,
                (double)(g_end - g_start) / (t1 - t0) / 1e6, threads);
    free(th); free(g_ring);
    return 0;
}
