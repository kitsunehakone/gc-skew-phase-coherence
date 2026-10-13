# Frozen protocol: prospective K4 harmonic-coherence test on untouched ALS vowel data

Frozen before waveform download/decoding or closure-outcome inspection.

## Source
Public dataset: Mak-Sim/Minsk2020_ALS_database
Frozen source commit: 4fa0cfcde7c919e1b68c99fdecd83ff9ac7334de

The earlier sound analyses used the 33 healthy-control (HC) recordings. This test uses the previously untouched ALS patient recordings.

## Primary dataset
Use ALS vowel /a/ recordings only.

Sorted ALS IDs:
008, 020, 021, 022, 024, 025, 027, 028, 031, 032, 039, 042, 046, 048, 052, 055, 058, 062, 064, 068, 072, 076, 078, 080, 084, 092, 094, 096, 098, 100, 102.

Deterministic alternating split before decoding:

Calibration IDs (positions 1,3,5,...):
008, 021, 024, 027, 031, 039, 046, 052, 058, 064, 072, 078, 084, 094, 098, 102.

Validation IDs (positions 2,4,6,...):
020, 022, 025, 028, 032, 042, 048, 055, 062, 068, 076, 080, 092, 096, 100.

No reassignment after outcomes are opened.

## Signal extraction
For each recording:
- use the central 1.0 s;
- 120 ms Hann frames with 60 ms hop;
- estimate F0 by normalized autocorrelation over 70–350 Hz;
- accept frame only if normalized autocorrelation peak >= 0.60;
- compute complex harmonics c_n, n=1..4, by Hann-weighted complex demodulation at n*F0;
- require 4*F0 < 0.9*Nyquist;
- no phase unwrapping across frames/recordings.

A recording is valid with at least 8 accepted frames.

## K4 closures
phi_n = arg(c_n)

kappa1 = wrap(phi_1 - 2 phi_2 + phi_3)
kappa2 = wrap(phi_2 - 2 phi_3 + phi_4)
chi    = wrap(kappa1 + kappa2)
       = wrap(phi_1 - phi_2 - phi_3 + phi_4)

## Frozen minimal model
The leading equal-edge nonchiral K4 model is

p_j(kappa1,kappa2) ∝ exp{ j [cos(kappa1) + cos(kappa2) + 2 cos(kappa1+kappa2)] }, j >= 0.

This encodes the frozen 1:1:2 closure-coupling prediction.
No phase offsets are fitted.
One global j is fitted on calibration data only.
No per-recording or per-closure coupling is allowed.

## Calibration
Fit j by maximum likelihood using accepted frames from calibration recordings only.
Use deterministic torus quadrature for the normalizer and predicted moments.
Restrict j to [0,20].

## Primary held-out quantities
On validation recordings compute:
M1   = mean cos(kappa1)
M2   = mean cos(kappa2)
Mchi = mean cos(chi)
Lval = mean held-out log density under calibration-fitted j.

Also report recording-clustered means as robustness diagnostics.

## Frozen primary success rule
The minimal one-j K4 model is provisionally supported only if all are true:
1. fitted j > 0;
2. held-out mean log density > uniform-torus log density, -2 log(2*pi);
3. Mchi > M1 and Mchi > M2 on held-out validation data;
4. Euclidean discrepancy between held-out empirical moment vector (M1,M2,Mchi) and model-predicted vector is <= 0.15.

If any condition fails, record the primary test as a failure of the minimal one-j K4 closure model for this dataset.

## Secondary diagnostics
These do not alter the primary verdict:
- fit unconstrained calibration model exp[j1 cos(k1)+j2 cos(k2)+j3 cos(k1+k2)];
- report fitted j1:j2:j3 and held-out likelihood difference versus frozen one-j model;
- report circular means/resultant lengths and valid-recording/frame counts;
- report small-fluctuation symmetric/antisymmetric variance ratio if concentration is sufficient;
- reserve ALS /i/ recordings for a later cross-vowel replication; do not inspect them during the primary /a/ test.
