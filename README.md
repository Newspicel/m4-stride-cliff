# m4-stride-cliff

A small C program that shows a performance cliff I found on Apple M4 while making
[RustFFT](https://github.com/ejmahler/RustFFT) faster ([PR #207](https://github.com/ejmahler/RustFFT/pull/207)).

It runs the five passes of a radix-8 FFT over 32768 complex numbers (512 KiB per buffer), the
same way RustFFT does, and prints the median time of every pass. The last pass reads 8 places in
memory and writes 8 places in memory, each exactly 64 KiB apart.

## Run it

You need a Mac with Apple Silicon and the Xcode command line tools.

```
git clone https://github.com/Newspicel/m4-stride-cliff
cd m4-stride-cliff
./run.sh
```

It takes about a minute. `./cliff A`, `./cliff B --skew 128` and `./cliff A --dit` run a single case.

## What it shows

On my M4 Max the last pass gets about 30% faster if the 8 places it writes to are moved 128 bytes
further apart from each other (`--skew 128`). On my M1 the same change does almost nothing.

| last pass, median | M4 Max | M1 |
|---|---|---|
| first pass A | 19.0 µs | 17.0 µs |
| first pass B | 18.0 µs | 17.0 µs |
| first pass A or B, `--skew 128` | 13.3 µs | 16.6 µs |

So without the skew the M4 Max is slower than the M1 in this pass, and with it, it's a lot faster.

There are two versions of the first pass. They compute the same thing, A just stores its twiddle
factors in a different layout. The later passes are exactly the same code in both cases, but on the
M4 they're about 5% slower after version A. On the M1 there's no difference.

## What I think is going on

With 64 KiB between them, all 8 writes land in the same set of the L1 cache, and that set only has 8
places. The CPU's prefetcher (the part that loads memory early because it expects you to need it)
has a hard time there. In the CPU performance counters the slow cases have a lot more L1 load misses
for the same instructions.

## The bigger version in RustFFT

In RustFFT itself the effect is bigger and stranger. On my M4 Max the last pass of a 32K-point f64
FFT takes either 15.5 µs or 18.8 µs, with the same code, the same data and bit-identical results.
Once it's slow it stays slow for a few milliseconds. My M1 doesn't do this (under 2%).

What I found out about it:

- The fast case only happens when the first pass and the last pass run through the very same
  function. The first pass reads the user's buffer and the last pass writes that same buffer, in
  the same pattern of 8 streams 64 KiB apart.
- Pretty much any change puts it into the slow case: a separate function for the first pass, a
  copy of the same code at another address, writing the last pass into a different buffer, a
  memory barrier (`dsb sy`) between the passes, or setting the DIT flag. DIT is known to turn off
  the data memory-dependent prefetcher and the load value predictor on the M3, nobody seems to have
  published what it does on the M4.
- Software prefetch (`prfm`) in the last pass doesn't help at all, at any distance. Neither does
  issuing the next loads before the current stores. So it's not just the prefetcher missing data
  or loads waiting behind stores.
- In the slow case the counters show about 25% more L1 load misses and about 50% more cycles where
  the oldest load waits for one. The clock and the cache latencies are the same in both cases.
- The kernel doesn't seem to be involved: macOS doesn't switch prefetcher or predictor settings
  per thread, so this has to come from the firmware setup or the hardware itself.

So the fast case looks like some state that the CPU learns in the first pass, tied to the
instruction addresses, which then makes the last pass faster. I don't know which part of the CPU
that is. This program doesn't show that part of the effect, only the smaller difference between
the two first passes and the big effect of the skew.

## Help

If you have an M2, M3, M4 (non-Max), M5 or an A-series device, I'd like to see your results. Please
open an issue with the output of `./run.sh`. Same if you know more about how the M4 prefetchers work.
