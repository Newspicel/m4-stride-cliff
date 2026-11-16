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
has a hard time there, and on the M4 it seems to depend on what ran before how well it copes. In the
CPU performance counters the slow cases have a lot more L1 load misses for the same instructions.
Apple doesn't document how its prefetchers work, so this is my best guess and not proven.

In RustFFT itself the effect is bigger: the same last pass took 15.5 µs or 18.8 µs depending on the
code of the first pass, and stayed slow for a few milliseconds afterwards. There, setting the DIT
flag (which on M3 is known to turn off the data memory-dependent prefetcher and the load value
predictor) also made the fast case slow in most of my runs. This program doesn't show the DIT part,
I haven't found out yet which detail of the real code is needed for it.

## Help

If you have an M2, M3, M4 (non-Max), M5 or an A-series device, I'd like to see your results. Please
open an issue with the output of `./run.sh`. Same if you know more about how the M4 prefetchers work.
