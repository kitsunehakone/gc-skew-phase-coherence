# Blinded Sound Primitive-Closure Validation — Zenodo Vowel Cohort

Frozen before downloading or decoding any selected waveform.

## Purpose
Test whether natural sustained human voice shows positive primitive harmonic phase-closure alignment using the same gauge-invariant closure statistic previously derived for GC-skew and coherent harmonic light.

A positive result supports a shared closure architecture. It does not by itself establish identical microscopic physics.

## External dataset
Zenodo record: 7098345
DOI: 10.5281/zenodo.7098345
Dataset description: sustained vowel /a/ phonation recordings, one recording per participant.

## Primary cohort
Use every vowel-phonation WAV in the Zenodo record whose filename does not contain "_phone".
The frozen primary filenames are:

P01_VWL.wav
P03_VWL.wav
P04_VWL.wav
P06_VWL.wav
P09_VWL.wav
P10_VWL.wav
P13_VWL.wav
P15_VWL.wav
P17_VWL.wav
P20_VWL.wav
P21_VWL.wav
P22_VWL.wav
P24_VWL.wav
P26_VWL.wav
P27_VWL.wav
P29_VWL.wav
P32_VWL.wav
P33_VWL.wav
P36_VWL.wav
P37_VWL.wav
P38_VWL.wav
P42_VWL.wav
P49_VWL.wav

No waveform is removed for statistical reasons.

## Audio validity gate
A file is valid only if:
1. It decodes as WAV.
2. Sampling rate is at least 8,000 Hz.
3. After channel averaging, duration is at least 1.0 s.
4. A central 1.0-s segment can be extracted.
5. At least 8 analysis frames pass the pitch-quality gate.

No invalid file is replaced.
The cohort test is considered interpretable only if at least 18 of the 23 frozen files are valid.

## Frozen signal analysis
For each valid WAV:
1. Convert samples to floating point.
2. Average channels if needed.
3. Extract exactly the central 1.0 s.
4. Subtract the segment mean.
5. Divide into 120-ms frames with 60-ms hop.
6. Estimate F0 in each frame by FFT autocorrelation.
7. Search F0 only in 70-350 Hz.
8. Choose the maximum autocorrelation peak in that range and refine its lag by parabolic interpolation.
9. Accept a frame only if peak autocorrelation divided by zero-lag autocorrelation is at least 0.60.
10. Estimate complex harmonic coefficients c_n for n=1,...,20 using Hann-windowed complex demodulation at n*F0:
    c_n = sum_t h_t x_t exp(-i 2*pi*n*F0*t) / sum_t h_t.

No phase smoothing, phase unwrapping, adaptive harmonic count, post-opening frequency-range change, or threshold tuning is allowed.

## Primitive closure
For every accepted frame f and harmonic center n=2,...,19:

  kappa_(f,n) = arg(c_(f,n+1)) + arg(c_(f,n-1)) - 2 arg(c_(f,n))

  w_(f,n) = |c_(f,n-1)| |c_(f,n)|^2 |c_(f,n+1)|

Pool all frame/harmonic triplets in a recording:

  M1 = [sum_(f,n) w_(f,n) exp(i*kappa_(f,n))] / [sum_(f,n) w_(f,n)]

Primary per-recording statistic:

  S1 = Re(M1)

The frozen directional prediction is positive closure alignment:
  S1 > S1_null
corresponding to primitive closure phases preferentially near kappa=0.

## Frozen null
B=4,999 phase-randomized surrogates per valid recording.

For every accepted frame, each surrogate preserves:
- the exact 20 harmonic magnitudes;
- the exact complex fundamental c_1;

and independently randomizes phases of c_2,...,c_20 uniformly on [0,2*pi).

Observed weights are retained exactly in every surrogate.

Per-recording one-sided Monte Carlo p-value:

  p_i = [1 + #{S1_null >= S1_obs}] / 5000

Minimum p=0.0002.

## Deterministic seeds
Sort the frozen filenames lexicographically.
For 1-based index r:

  seed_r = 20261012 + 100003*r

## Primary cohort endpoint
For N valid recordings:

  X = -2 sum_i log(p_i)

with 2N chi-square degrees of freedom.

The blinded sound validation is positive iff:
1. N >= 18, and
2. p_Fisher < 0.05.

This is the sole primary success criterion.

## Secondary outputs
- per-recording S1 and p_i;
- BH q-values;
- mean/median S1;
- |M1| and arg(M1);
- median accepted-frame F0;
- accepted-frame count;
- effective closure-weight count.

No secondary statistic may replace S1 after unblinding.
