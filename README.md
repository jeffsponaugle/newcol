# collatzdig

Enumerates n such that the digits of n (counted with multiplicity) are **not** all
encountered in the Collatz trajectory n -> ... -> 1, excluding the initial n.

    5  is a member: 16 8 4 2 1 contains no 5
    29 is a member: the 2 is covered but a 9 never appears
    77 is a member: only one 7 appears (in 17); two are needed
    30 is not:      15 46 23 70 covers both the 3 and the 0

## Build and test

    make          # builds ./collatzdig
    make test     # checks output against the 216 posted terms through 20,000,000

## Run

    ./collatzdig 1 1000000000            # terms in [1, 10^9], one per line on stdout
    ./collatzdig -t 8 -p 10 1 1e9        # 8 threads, progress every 10 s on stderr
    ./collatzdig 10B 1T                  # suffixes: K M B/G T P
    ./collatzdig -x 77                   # show trajectory with digits still needed

Options: `-t` threads (default: all CPUs), `-b` block size, `-p` status refresh
interval in seconds (default 1, 0 = off), `-q` quiet. Output is always in increasing order.

While running, a status line on stderr shows percent done, position, terms found,
current rate, elapsed time and ETA. On a terminal it refreshes in place; when
stderr is redirected to a file it writes one line per interval instead.

Numbers (start, end, block size, `-x`) accept plain integers, hex (`0x1000`),
scientific notation (`1e14`, `2.5e9`) and scale suffixes, case-insensitive:
`K`=1e3, `M`=1e6, `B`/`G`=1e9, `T`=1e12, `P`=1e15, e.g. `10M`, `1.5B`, `3T`.
Values are parsed exactly in integer arithmetic; anything non-integral is rejected.

## Implementation notes

- The digit counts of n live in one 64-bit word as ten 6-bit lanes ("need").
  Each trajectory value is reduced 4 decimal digits at a time through a 10,000-entry
  table of packed digit counts, then subtracted from "need" lane-wise with a
  saturating branch-free SWAR step. need == 0 means every digit was covered, and the
  test stops there. The odd/even Collatz step is a branch-free select.
- Since the upper digits of n are constant across 10,000 consecutive n, a worker
  maintains that part of "need" incrementally and adds one table lookup per n.
- Typical n is decided within about 4 steps while its trajectory values are still the
  size of n, so caching trajectories of smaller numbers would almost never be hit;
  that is why the work went into per-step cost instead.
- Arithmetic is 64-bit, falling back to 128-bit for the rare trajectory that would
  overflow (`OVERFLOW_LIMIT`). `make collatzdig_wide` builds a test binary that
  forces the 128-bit path above 2^40.
- Work is handed out in blocks via an atomic counter; a small ring buffer re-orders
  finished blocks so stdout is printed in sequence order.
- Throughput on an M-series 14-core near n = 4e12: roughly 550-600 million n per second
  (about 50 million per core).
