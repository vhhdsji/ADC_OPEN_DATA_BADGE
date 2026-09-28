# A 12b 4GS/s Timing-Skew-Free Pipelined ADC with Branched-Interleaving Fully Charge-Domain Residue Transfer

This repository provides the measured raw data and MATLAB code used to calibrate the raw decisions and display the low-frequency (LF) and high-frequency (HF) FFT spectra reported in Figure 4 of the paper. It also includes a plotting report for Figure 5: the four measurement plots (input-frequency response, supply-voltage sweep, temperature sweep, and dynamic range) together with the power breakdown. These Figure 5 displays use five supplied MAT datasets.

The ADC uses one common full-rate 4 GS/s input sampler followed by a 1-2-4 branched-interleaved backend. The four paths processed by this code are the four 1 GS/s STG2/STG3 backend channels, not four independent front-end samplers. The architecture is therefore timing-skew-free at the ADC input, and the released processing does not perform input-channel timing-skew calibration.

The raw-data processing in `calibrate_figure4_fft.m` applies only linear bit-weight, interstage-gain, and channel-offset calibration. No nonlinear calibration, lookup table, or polynomial correction is used. The implementation references corresponding MATLAB ADCToolbox algorithms, and all files required at runtime are included locally. The additional `plot_figure5_report.m` function only displays values already stored in the five plotting datasets; it does not implement calibration, fitting, FFT analysis, or metric extraction.

## Running the Figure 4 Calibration and Figure 5 Report

Open MATLAB and change to the repository root. Each entry point runs independently:

| Paper figure | MATLAB entry point | Purpose |
| --- | --- | --- |
| Figure 4 | `calibrate_figure4_fft` | Process the LF/HF raw decision streams, apply linear calibration, calculate FFT metrics, and display the two spectra |
| Figure 5 | `plot_figure5_report` | Display the four measurement plots and power breakdown from supplied plotting datasets; no calibration or signal-processing algorithm is executed |

### Figure 4: LF/HF Calibration and FFT

Run:

```matlab
calibrate_figure4_fft
```

The script fits the foreground-calibration coefficients from the LF record, freezes those coefficients, applies them to both records, calculates SNDR and SFDR, and displays the two FFT spectra. The measured data, fitted weights, channel-calibration terms, reconstructed streams, FFT arrays, and metrics remain available in the MATLAB workspace.

### Figure 5: Measurement Plot Report

To display the four Figure 5 measurement plots and power breakdown, run:

```matlab
plot_figure5_report
```

This function loads the five plotting datasets from `data/`, draws the supplied values and annotations, and exports `figure5_report.fig`, `.png`, `.svg`, and `.pdf` into `plots/` (created automatically). It does not require running `calibrate_figure4_fft` first. The MATLAB installation must support `xline`, `savefig`, and the requested export formats.

## Project Structure

```text
adc_calibration/
|- calibrate_figure4_fft.m
|- plot_figure5_report.m
|- LICENSE
|- THIRD_PARTY_NOTICES.md
|- data/
|  |- lf_measured_raw.mat
|  |- hf_measured_raw.mat
|  |- frequency.mat
|  |- supply_voltage.mat
|  |- temperature.mat
|  |- dynamic_range_170MHz.mat
|  `- power_breakdown.mat
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

## Figure 5 Plotting Datasets

These five MAT files contain plotting values and presentation metadata, such as axis limits, ticks, labels, legend geometry, and annotation positions. They are used only by `plot_figure5_report.m` for visualization; they do not provide or execute a signal-processing or calibration algorithm.

| File | Displayed content |
| --- | --- |
| `data/frequency.mat` | SNDR and SFDR versus input frequency |
| `data/supply_voltage.mat` | Multi-chip SNDR versus supply voltage, showing stored one-time-calibration and recalibration results |
| `data/temperature.mat` | Multi-chip SNDR versus temperature, showing stored one-time-calibration and recalibration results |
| `data/dynamic_range_170MHz.mat` | SNDR/SFDR versus signal power at a 170 MHz input, dynamic-range annotation, and a full-scale inset |
| `data/power_breakdown.mat` | Power breakdown by circuit block and displayed total power |

Curve values and calibration-result labels are read from the MAT files. The plotting code performs only drawing-related operations, including unit conversion and power-sector proportions; it does not recompute the underlying measurement results.

## Figure 4 Measured Raw Data

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

After `calibrate_figure4_fft.m` finishes, the MATLAB workspace contains:

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
