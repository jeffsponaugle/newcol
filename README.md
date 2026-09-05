# collatzdig

Enumerates n such that the digits of n (counted with multiplicity) are **not** all
encountered in the Collatz trajectory n -> ... -> 1, excluding the initial n.

    5  is a member: 16 8 4 2 1 contains no 5
    29 is a member: the 2 is covered but a 9 never appears
    77 is a member: only one 7 appears (in 17); two are needed
    30 is not:      15 46 23 70 covers both the 3 and the 0

The first terms are 1, 2, 3, 4, 5, 8, 9, 10, 13, 16, 17, 29, 32, 49, 53, 69, 70, 77, ...
There are 216 terms below 2e7, 229 below 1e10, and the density keeps falling.
There are 231 terms below 1e12 and the last term I have found to date is 15,078,777,741

## Build and test

    make          # builds ./collatzdig
    make test     # checks output against the 216 posted terms through 20,000,000

Requires a C11 compiler with `unsigned __int128` (clang or gcc) and pthreads.

## Run

    ./collatzdig 1 1000000000            # terms in [1, 10^9], one per line on stdout
    ./collatzdig -t 8 -p 10 1 1e9        # 8 threads, status every 10 s on stderr
    ./collatzdig 10B 1T                  # suffixes: K M B/G T P
    ./collatzdig -x 77                   # show trajectory with digits still needed

Options: `-t` threads (default: all CPUs), `-b` block size (default 65536),
`-p` status refresh interval in seconds (default 1, 0 = off), `-q` quiet.
Output is always in increasing order, so a run can be appended to an earlier one.

While running, a status line on stderr shows percent done, position, terms found,
current rate, elapsed time and ETA. On a terminal it refreshes in place; when
stderr is redirected to a file it writes one line per interval instead.

Numbers (start, end, block size, `-x`) accept plain integers, hex (`0x1000`),
scientific notation (`1e14`, `2.5e9`) and scale suffixes, case-insensitive:
`K`=1e3, `M`=1e6, `B`/`G`=1e9, `T`=1e12, `P`=1e15, e.g. `10M`, `1.5B`, `3T`.
Values are parsed exactly in integer arithmetic; anything non-integral is rejected.

## How it works

### The test

For each n, keep a count of how many of each digit still need to be seen, then walk
the Collatz trajectory and subtract the digits of every value encountered. If every
count reaches zero the trajectory covers n and n is not a member. If the walk reaches
1 with something still outstanding, n is a member.

The walk stops the moment coverage is complete. Measured near n = 4e12, the average n
is decided after 4.4 steps, 99% within 10 steps, and only members (very rare) need the
full trajectory of a few hundred steps. So the program's cost is almost entirely the
first handful of steps, and everything below is about making those steps cheap.

### Packed digit counters (SWAR)

The ten "still needed" counts live in one 64-bit word, `need`, as ten 6-bit lanes:
lane d occupies bits 6d..6d+5 and holds the count for digit d. A trajectory value is
converted to the same packed form, `have`, and the update is a lane-wise saturating
subtraction, `need = max(need - have, 0)` per lane, performed on all ten lanes at once
with plain integer arithmetic (SIMD within a register):

    d    = (need | H) - have          H = 0x20 in every lane (bit 5, the guard bit)
    need = d & (((d & H) >> 5) * 0x1F)

Setting the guard bit before subtracting guarantees no lane borrows from its neighbour
(every lane of `have` is at most 31). After the subtraction the guard bit is still set
exactly in the lanes where `need >= have`. Shifting those surviving bits down to lane
bit 0 and multiplying by 0x1F spreads each into a 5-bit mask that keeps `need - have`
in those lanes and zeroes the lanes where `have` exceeded `need`. The whole update is
five ALU operations with no branches and no loops, and `need == 0` is the coverage test.

Six-bit lanes leave headroom: `need` lanes are at most the digit count of n (up to 20
for a 64-bit n) and `have` lanes are at most 19 or 20 digits, all below the 32 that
the guard bit requires.

### Digits four at a time

Converting a value into packed form is done in chunks of four decimal digits. A
10,000-entry table, `cnt4`, maps a chunk (with leading zeros) to its packed digit
count, so a 13-digit value costs three `% 10000` / `/ 10000` steps and four table
loads that are simply added together. The table is 80 KB and stays in L1/L2 cache.
The most significant chunk is written without leading zeros, so the zeros that `cnt4`
counted for a short top chunk are subtracted from lane 0 with a comparison count,
`(y < 10) + (y < 100) + (y < 1000)`, rather than a second table.

### Branch-free Collatz step

The odd/even decision is essentially a coin flip, so a conditional branch would be
mispredicted half the time at roughly 15 cycles each. The step is written as a select,
`x = (x & 1) ? 3*x + 1 : x >> 1`, which the compiler emits as a conditional move.
The only branches left per step are the coverage test and a rarely taken overflow check.

### Incremental digit counts of n

The digit counts of n itself would otherwise cost a 13-step division chain per n.
Within a run of 10,000 consecutive numbers the upper digits n / 10000 are constant, so
each worker computes their packed count once per 10,000 and adds `cnt4[n % 10000]`
for each n: one table lookup. Numbers below 10,000 are handled directly so leading
zeros are not counted.

### Overflow handling

Trajectory arithmetic is 64-bit. Before a step, if x is large enough that 3x + 1
could overflow (`OVERFLOW_LIMIT`), the trajectory continues in 128-bit arithmetic
until it drops back below 2^62, then returns to the fast path. Digits of a 128-bit
value are consumed in 19-digit pieces so no lane of `have` exceeds 31. This path is
extremely rare: no start value below about 8.5e9 even approaches 2^64. Because that
makes it hard to exercise naturally, `make collatzdig_wide` builds a test binary with
the threshold lowered to 2^40, which reproduces the reference list through the same
code (`make test` runs it).

### Threads and ordered output

Workers take fixed-size blocks of n from an atomic counter, so load balances itself
regardless of how uneven the blocks are. Each finished block's terms are handed to a
small ring buffer keyed by block index, and whichever worker completes the next
block in sequence flushes all consecutive ready blocks to stdout. Output is therefore
in increasing order even though blocks complete out of order. A worker whose block is
too far ahead of the printed frontier waits on a condition variable; the worker holding
the frontier block never waits, so the scheme cannot deadlock.

### What did not help

- **Caching trajectories of smaller numbers.** Since n is nearly always decided within
  a few steps while the trajectory values are still the size of n, a cache indexed by
  small values would almost never be consulted before the decision. It would cost memory
  bandwidth and return nothing.
- **Splitting values at 10^8** to reduce the low half with independent 32-bit
  divisions instead of one 64-bit chain. It measured within noise of the simple loop,
  so the simpler code stayed.
- **A presence-bitmask precheck** (the original design) that skipped chunks sharing no
  digit with the needed set. Early in a trajectory most chunks do intersect, so the
  branchy per-digit fallback ran most of the time; the SWAR update replaced both.

## Performance

Measured on an Apple M-series 14-core machine, n near 4e12, identical output:

| version                                  | 1 thread  | 14 threads |
|------------------------------------------|-----------|------------|
| original (bitmask + per-digit loop)      | 8.7 M/s   | 86 M/s     |
| packed counters, branch-free step        | 36 M/s    | 385 M/s    |
| + incremental digit counts of n          | 50 M/s    | 580 M/s    |

On a Mac Studio Ultra with 24 cores it gets about 900 M/s.

Rates fall slowly as n grows because values have more digits and take slightly longer
to cover. At 1e13 to 1e14 expect a few hundred million n per second on such a machine.

## Verification

- Output for 1..2e7 matches the 216 terms posted to seqfan exactly.
- The forced-128-bit build reproduces the same list.
- Odd block sizes and unaligned start values give identical output to the original
  straightforward implementation.
