# A 12b 4GS/s Timing-Skew-Free Pipelined ADC with Branched-Interleaving Fully Charge-Domain Residue Transfer

This repository provides the measured raw data and MATLAB code used to reproduce the low-frequency (LF) and high-frequency (HF) FFT spectra reported in Figure 4 of the paper. Only these two dynamic-performance records are included.

The ADC uses one common full-rate 4 GS/s input sampler followed by a 1-2-4 branched-interleaved backend. The four paths processed by this code are the four 1 GS/s STG2/STG3 backend channels, not four independent front-end samplers. The architecture is therefore timing-skew-free at the ADC input, and the released processing does not perform input-channel timing-skew calibration.

The processing applies only linear bit-weight, interstage-gain, and channel-offset calibration. No nonlinear calibration, lookup table, or polynomial correction is used. The implementation references corresponding MATLAB ADCToolbox algorithms, and all files required at runtime are included locally.

## Running the Project

Open MATLAB, change to the project directory, and run:

```matlab
run_data
```

The script fits the foreground-calibration coefficients from the LF record, freezes those coefficients, applies them to both records, calculates SNDR and SFDR, and displays the two FFT spectra. The measured data, fitted weights, channel-calibration terms, reconstructed streams, FFT arrays, and metrics remain available in the MATLAB workspace.

## Project Structure

```text
adc_calibration/
|- run_data.m
|- LICENSE
|- THIRD_PARTY_NOTICES.md
|- data/
|  |- lf_measured_raw.mat
|  `- hf_measured_raw.mat
`- src/
   |- wcalsin.m
   |- sinfit.m
   |- findfreq.m
   |- findbin.m
   |- fit_sine.m
   |- fft_metrics.m
   |- alias_frequency.m
   |- local_rms.m
   |- rms.m
   `- plot_fft.m
```

## Measured Raw Data

Both datasets were measured from the same chip at room temperature. Each file contains an 8192-sample coherent record.

### `data/lf_measured_raw.mat`

Measured low-input-frequency data used for the one-time foreground calibration and the LF spectrum.

```text
DATA.bits                         8192 x 20 uint8 binary decision stream
DATA.config.Fs_Hz                 4 GHz sampling rate
DATA.config.Fin_Hz                170.41015625 MHz input frequency
DATA.config.normalized_frequency  Fin/Fs
DATA.config.tone_bin              349
```

### `data/hf_measured_raw.mat`

Measured high-input-frequency validation data used for the HF spectrum. No weights are refitted from this dataset; the coefficients obtained from the LF record are applied directly.

```text
DATA.bits                         8192 x 20 uint8 binary decision stream
DATA.config.Fs_Hz                 4 GHz sampling rate
DATA.config.Fin_Hz                1957.03125 MHz input frequency
DATA.config.normalized_frequency  Fin/Fs
DATA.config.tone_bin              4008
```

Both records satisfy:

```matlab
normalized_frequency * 8192 == tone_bin
```

## Raw Decision-Column Mapping

The 20 columns are internal stage decisions, not a 20-bit ADC output word:

```text
bits(:,1:8)    STG1 3-bit flash complete 8-level thermometer-code representation
bits(:,9:12)   STG2 4-bit SAR outputs
bits(:,13:20)  STG3 7-bit SAR outputs plus 1 in-stage redundant decision bit
```

## Calibration Workflow

### 1. Split the Four Branched Backend Paths

The full-rate record is separated by its four-sample backend cycle:

```matlab
channel_1 = bits(1:4:end,:);
channel_2 = bits(2:4:end,:);
channel_3 = bits(3:4:end,:);
channel_4 = bits(4:4:end,:);
```

Each path contains 2048 samples. Its calibration frequency is obtained by folding `4*Fin/Fs` into the first Nyquist zone. This split identifies the four STG2/STG3 backend paths and does not imply four independent input samplers.

### 2. Fit Four Independent 20-Column Weight Sets

`wcalsin.m` performs a sine fit independently for each backend path. Its linear least-squares model uses the 20 decision columns and sine/cosine basis functions, producing one 20-element weight vector per path:

```text
weight(1,:)  backend path 1 weights
weight(2,:)  backend path 2 weights
weight(3,:)  backend path 3 weights
weight(4,:)  backend path 4 weights
```

The fitted weights correct linear CDAC weight mismatch and absorb the linear interstage-gain terms visible in the end-to-end decisions. Although STG1 is a common full-rate stage, separate end-to-end fits are retained so that each backend path has its own linear reconstruction coefficients. Only the LF dataset is used to fit these coefficients.

### 3. Correct Backend-Path Offset and Linear Gain

Each path is reconstructed using its fitted weights:

```matlab
weighted = bits * weight.';
```

`fit_sine.m` obtains the residual DC offset and sine amplitude of each LF path:

```matlab
gain = 1 / fitted_amplitude;
postcal = (weighted - fitted_offset) * gain;
```

The LF-derived offset and gain terms are frozen and applied to the HF record together with the LF-derived bit weights.

### 4. Merge the Four Backend Paths

The calibrated paths are placed back in the original full-rate sample order:

```matlab
postcal(channel:4:end) = channel_postcal{channel};
```

## FFT, SNDR, and SFDR

`fft_metrics.m` performs the following operations on each merged stream:

1. Removes the DC mean.
2. Computes an 8192-point FFT.
3. Generates a one-sided spectrum from DC to Nyquist.
4. Locates the coherent input tone using `tone_bin`.
5. Excludes DC and the fundamental and treats the remaining spectral power as noise and distortion.
6. Uses the largest non-fundamental spectral line to calculate SFDR.

```text
SNDR = fundamental power / all spectral power excluding DC and the fundamental
SFDR = fundamental power / largest spur power
```

The plots use `Fs = 4 GS/s`, 8192 points, a frequency axis from `0` to `2000 MHz`, and an amplitude axis from `-120` to `0 dBFS`.

## Main Workspace Variables

After `run_data.m` finishes, the MATLAB workspace contains:

```text
LF / HF                         measured input structures
lf_bits / hf_bits               measured 8192 x 20 decision streams
weight                          four independent 20-element backend-path weight vectors
channel_fixed_offset            LF-fitted residual offset for each backend path
channel_fixed_amplitude         LF-fitted amplitude for each backend path
channel_gain                    linear amplitude-normalization gain for each backend path
lf_postcal / hf_postcal         calibrated and merged 8192-sample streams
lf_fft_db / hf_fft_db           one-sided FFT spectra
lf_sndr_db / lf_sfdr_db         LF dynamic-performance metrics
hf_sndr_db / hf_sfdr_db         HF dynamic-performance metrics
result_figure                   handle of the two-panel FFT figure
```

## License and Third-Party Attribution

This project is released under the MIT License. See `LICENSE` for the full license text.

Several calibration, sine-fitting, frequency-search, and bin-search algorithms reference or adapt ADCToolbox, which is also distributed under the MIT License. See `THIRD_PARTY_NOTICES.md` for its copyright, source, license, and citation information.
