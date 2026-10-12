# Blinded Sound Primitive-Closure Validation — Frozen Protocol

Frozen before reading or decoding any selected audio waveform.

## Purpose
Test whether natural sustained human voice shows the same primitive harmonic phase-closure alignment previously derived and independently validated on bacterial GC-skew signals.

This sound test evaluates the generic spectral closure object only. A positive result does not establish identical microscopic physics between sound, light, and replication.

## External dataset
Source repository: Mak-Sim/Minsk2020_ALS_database
Frozen source commit: 4fa0cfcde7c919e1b68c99fdecd83ff9ac7334de

The source publication/database describes 128 sustained vowel phonations from 64 speakers, recorded as 44.1 kHz, 16-bit PCM. The present primary cohort uses only healthy-control vowel /a/ recordings.

## Primary cohort
All 33 files under the frozen source commit matching HC/*_a.wav, fixed before waveform inspection:

HC/002_a.wav
HC/004_a.wav
HC/006_a.wav
HC/016_a.wav
HC/024_a.wav
HC/026_a.wav
HC/028_a.wav
HC/042_a.wav
HC/049_a.wav
HC/053_a.wav
HC/061_a.wav
HC/063_a.wav
HC/065_a.wav
HC/077_a.wav
HC/081_a.wav
HC/085_a.wav
HC/086_a.wav
HC/089_a.wav
HC/097_a.wav
HC/099_a.wav
HC/107_a.wav
HC/109_a.wav
HC/111_a.wav
HC/113_a.wav
HC/115_a.wav
HC/117_a.wav
HC/119_a.wav
HC/121_a.wav
HC/123_a.wav
HC/125_a.wav
HC/127_a.wav
HC/129_a.wav
HC/131_a.wav

The HC/*_i.wav recordings are untouched by this primary test and are reserved for a later cross-vowel replication.

## Audio/provenance gate
A primary file is valid only if:
1. It decodes as PCM WAV.
2. Sampling rate is exactly 44,100 Hz.
3. It is mono.
4. Duration is at least 1.0 s.
5. A central 1.0-s segment can be extracted.
6. At least 8 analysis frames pass the frozen pitch-quality gate below.

No invalid file is replaced. The cohort-level test is considered interpretable only if at least 25 of 33 files are valid.

## Frozen signal analysis
For each WAV:
1. Convert PCM samples to floating point.
2. Extract exactly the central 1.0 s.
3. Subtract the segment mean.
4. Divide into 120-ms frames with 60-ms hop.
5. For each frame, subtract its mean and estimate F0 by FFT autocorrelation.
6. Search lags corresponding to 70-350 Hz.
7. Choose the maximum autocorrelation lag in that range and refine it by parabolic interpolation.
8. Accept the frame only if the autocorrelation peak divided by zero-lag autocorrelation is at least 0.60.
9. For each accepted frame estimate complex harmonic coefficients c_n, n=1,...,20, at frequencies n*F0 using Hann-windowed complex demodulation:
   c_n = sum_t h_t x_t exp(-i 2*pi*n*F0*t) / sum_t h_t.

No harmonic-phase smoothing, phase unwrapping, tuning, or post-opening threshold change is allowed.

## Primitive closure statistic
For every accepted frame f and harmonic center n=2,...,19:

  kappa_(f,n) = arg(c_(f,n+1)) + arg(c_(f,n-1)) - 2 arg(c_(f,n))

and

  w_(f,n) = |c_(f,n-1)| |c_(f,n)|^2 |c_(f,n+1)|.

Pool all accepted frame/harmonic triplets in one recording:

  M1 = [sum_(f,n) w_(f,n) exp(i*kappa_(f,n))] / [sum_(f,n) w_(f,n)]

Primary per-recording statistic:

  S1 = Re(M1).

The prespecified directional prediction is positive closure alignment:

  S1 > S1_null,

corresponding to primitive closure phases preferentially near kappa=0.

## Frozen phase-randomized null
B = 4,999 surrogates per valid recording.

For every accepted frame separately, every surrogate preserves:
- all 20 observed harmonic magnitudes exactly;
- the exact complex fundamental c_1 exactly;

and randomizes phases of c_2,...,c_20 independently and uniformly on [0,2*pi).

The same observed weights are used for observed and null statistics.

Per-recording one-sided Monte Carlo p-value:

  p_i = [1 + #{S1_null >= S1_obs}] / 5000.

Minimum p is 0.0002.

## Deterministic random seeds
Base seed = 20261012.
After sorting the 33 frozen paths lexicographically, recording index r is 1,...,33 and uses:

  seed_r = 20261012 + 100003*r.

## Primary cohort-level endpoint
For N valid recordings, combine per-recording p-values by Fisher's method:

  X = -2 sum_i log(p_i)

with 2N chi-square degrees of freedom.

The blinded sound validation is positive iff:
1. N >= 25 valid recordings, and
2. p_Fisher < 0.05.

This is the sole primary success criterion.

## Secondary/descriptive outputs
- number of raw p<0.05 and BH q<0.05 recordings;
- mean and median S1;
- |M1| and arg(M1);
- estimated median F0 per recording;
- accepted-frame count;
- effective closure-weight count;
- per-recording and cohort distributions.

No secondary statistic may replace S1 after waveform opening.

## Reserved test
All healthy-control /i/ recordings remain reserved and are not analyzed in this run.
