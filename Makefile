CC      ?= cc
CFLAGS  ?= -O3 -march=native -Wall -Wextra -pthread

collatzdig: collatzdig.c
	$(CC) $(CFLAGS) -o $@ $<

# test build: uses 128-bit arithmetic for anything above 2^40 (exercises the wide path)
collatzdig_wide: collatzdig.c
	$(CC) $(CFLAGS) -DOVERFLOW_LIMIT='(1ULL<<40)' -o $@ $<

test: collatzdig collatzdig_wide
	./collatzdig -q 1 20000000 | diff expected_lines.txt - && echo "ok: matches posted list"
	./collatzdig_wide -q 1 20000000 | diff expected_lines.txt - && echo "ok: wide path matches too"

clean:
	rm -f collatzdig collatzdig_wide

.PHONY: test clean
