# Blinded Sound /i/ Cross-Vowel Primitive-Closure Replication — Frozen Protocol

Frozen before decoding or analyzing any healthy-control /i/ waveform.

## Purpose
Cross-vowel replication of the primitive harmonic phase-closure signal found in the frozen healthy-control /a/ test, using the corresponding healthy-control /i/ recordings from the same frozen external database.

This is a within-speaker/cross-vowel replication, not an independent-speaker cohort.

## External dataset
Source repository: Mak-Sim/Minsk2020_ALS_database
Frozen source commit: 4fa0cfcde7c919e1b68c99fdecd83ff9ac7334de

## Frozen cohort
All 33 healthy-control files matching HC/*_i.wav at the frozen source commit:

HC/002_i.wav
HC/004_i.wav
HC/006_i.wav
HC/016_i.wav
HC/024_i.wav
HC/026_i.wav
HC/028_i.wav
HC/042_i.wav
HC/049_i.wav
HC/053_i.wav
HC/061_i.wav
HC/063_i.wav
HC/065_i.wav
HC/077_i.wav
HC/081_i.wav
HC/085_i.wav
HC/086_i.wav
HC/089_i.wav
HC/097_i.wav
HC/099_i.wav
HC/107_i.wav
HC/109_i.wav
HC/111_i.wav
HC/113_i.wav
HC/115_i.wav
HC/117_i.wav
HC/119_i.wav
HC/121_i.wav
HC/123_i.wav
HC/125_i.wav
HC/127_i.wav
HC/129_i.wav
HC/131_i.wav

No file may be removed or replaced for statistical reasons.

## Frozen analysis
Use exactly the same signal-processing and statistical pipeline as the blinded /a/ validation:

1. PCM WAV, mono, 44,100 Hz.
2. Central 1.0-s segment.
3. 120-ms frames, 60-ms hop.
4. F0 by FFT autocorrelation in 70-350 Hz.
5. Frame accepted only if normalized autocorrelation peak >= 0.60.
6. Complex harmonic coefficients c_n for n=1,...,20 by Hann-windowed complex demodulation at n*F0.
7. At least 8 accepted frames required.

For each accepted frame f and harmonic center n=2,...,19:

  kappa_(f,n) = arg(c_(f,n+1)) + arg(c_(f,n-1)) - 2 arg(c_(f,n))

  w_(f,n) = |c_(f,n-1)| |c_(f,n)|^2 |c_(f,n+1)|

Pool over all accepted frame/harmonic triplets:

  M1 = [sum w exp(i*kappa)] / [sum w]

Primary recording statistic:

  S1 = Re(M1)

Directional prediction:

  S1 > S1_null.

No tuning, alternate phase convention, harmonic-range change, adaptive smoothing, or threshold change is permitted after waveform opening.

## Frozen null
B=4,999 surrogates per valid recording.

For every accepted frame, each surrogate preserves:
- all 20 harmonic magnitudes exactly;
- the exact complex fundamental c_1 exactly;

and independently randomizes the phases of c_2,...,c_20 uniformly on [0,2*pi).

Per-recording one-sided Monte Carlo p-value:

  p_i = [1 + #{S1_null >= S1_obs}] / 5000

## Seeds
Lexicographically sort the 33 frozen filenames.
For 1-based recording index r:

  seed_r = 20261012 + 100003*r

## Primary endpoint
For N valid recordings, Fisher combine the per-recording p-values:

  X = -2 sum log(p_i)

with 2N chi-square degrees of freedom.

Replication is positive iff:
1. N >= 25, and
2. p_Fisher < 0.05.

This is the sole primary success criterion.

## Secondary/descriptive
- raw p<0.05 count;
- BH q<0.05 count;
- mean and median S1;
- |M1| and arg(M1);
- median F0;
- accepted-frame count;
- effective closure-weight count.

No secondary statistic may replace S1.
