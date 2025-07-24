use core::f64;

use ndarray::{s, Array1};

use crate::{error::Error, BIN_SHIFT};

/// According to the original authors, after binning the experimental spectrum, the maximum intensity is normalized
/// over a number of fixed windows.
const NUM_WINDOWS_FOR_NORMALIZATION: u8 = 10;

/// Calculates the binned position for a given m/z value.
/// This is used to determine which bin a given m/z value falls into.
///
/// Arguments:
/// * `mz` - The m/z value to be binned.
/// * `bin_offset_mz` - The offset in m/z to be applied before bin
/// * `bin_size` - The size of the bins to be used for binning the spectrum.
/// * `shift` - The number of bins to add to the start and end of the
///
pub fn calc_binned_position(mz: f64, bin_offset_mz: f64, bin_size: f64, shift: usize) -> usize {
    ((mz + bin_offset_mz) / bin_size).floor() as usize + shift - 1
}

/// Bins the theoretical spectrum based on the m/z values. Add bins for the shifts which avoids resizing the array later.
///
/// Arguments:
/// * `mz` - The m/z values of the theoretical spectrum.
/// * `bin_size` - The size of the bins to be used for binning the spectrum.
/// * `bin_offset` - The offset in monoisotopic mass. This is not used in the binning, but it is included for consistency with the experimental spectrum binning.
/// * `mz_max` - The maximum m/z value to consider for the spectrum. If `None`, the maximum m/z value from the input will be used.
/// * `use_flanking_peaks` - Whether to use flanking peaks in the binning.
/// * `apply_shift` - Whether to apply the BIN_SHIFT to the binned spectrum
///
pub fn theoretical_spectrum_binning(
    mz: &Array1<f64>,
    bin_size: f64,
    bin_offset: f64,
    mz_max: Option<f64>,
    use_flanking_peaks: bool,
    apply_shift: bool,
) -> Result<Array1<f64>, Error> {
    if mz.is_empty() {
        return Err(Error::EmptyTheoreticalSpectrum);
    }

    // bins_filled = np.zeros(math.ceil(np.max(mz_array) / bin_width) + 1)
    let mz_max = match mz_max {
        Some(max) => max,
        None => *mz.last().unwrap(),
    };

    let number_of_bins = (mz_max / bin_size).ceil() as usize;
    let shift = if apply_shift { BIN_SHIFT } else { 0 };
    let mut bins: Array1<f64> = Array1::zeros(number_of_bins + 2 * shift + 1);

    //  for mass in mz_array:
    for mz in mz.iter() {
        // index = int(mass // bin_width)
        let index = (*mz / bin_size + bin_offset).floor() as usize + shift - 1;
        // bins_filled[index] = 50.0
        bins[index] = 50.0;

        if use_flanking_peaks {
            // if index - 1 != -1:
            //     bins_filled[index - 1] = max(bins_filled[index - 1], 25.0)
            if index > 0 {
                bins[index - 1] = bins[index - 1].max(25.0);
            }
            if index + 1 < bins.len() {
                bins[index + 1] = bins[index + 1].max(25.0);
            }
        }
    }

    Ok(bins)
}

/// Bins and normalizes the experimental spectrum.
///
/// Arguments:
/// * `mz` - The m/z values of the experimental spectrum.
/// * `intensities` - The intensity values of the experimental spectrum.
/// * `bin_size` - The size of the bins to be used for binning the spectrum.
/// * `bin_offset` - The offset in monoisotopic mass
/// * `use_flanking_peaks` - Whether to use flanking peaks in the binning.
///
pub fn experimental_spectrum_binning(
    mz: &Array1<f64>,
    intensities: &Array1<f64>,
    bin_size: f64,
    bin_offset: f64,
    use_flanking_peaks: bool,
) -> Result<Array1<f64>, Error> {
    if mz.len() != intensities.len() {
        return Err(Error::ExperimentalSpectrumShape(
            mz.len(),
            intensities.len(),
        ));
    }

    // bins_filled = np.zeros(math.ceil(np.max(mz_array) / bin_width) + 1)
    let mz_max = mz.last().unwrap();
    let number_of_bins = (mz_max / bin_size).ceil() as usize;
    let mut binned_spectrum: Array1<f64> = Array1::zeros(number_of_bins + 2 * BIN_SHIFT + 1);

    // intensity_array = np.sqrt(intensity_array)
    // Doing this when adding to the bins, so we don't need to create a new array.

    // for mass, intensity in zip(mz_array, intensity_array):
    //     index = int(mass // bin_width)
    //     bins_filled[index] = max(bins_filled[index], intensity)
    for (mz, intensity) in mz.iter().zip(intensities.iter()) {
        let index = (*mz / bin_size + bin_offset).floor() as usize + BIN_SHIFT - 1;
        binned_spectrum[index] = binned_spectrum[index].max(intensity.sqrt());
    }

    // Normalization

    // highest_ion = bins_filled.size
    // num_wins = 10
    // win_size = int(highest_ion/num_wins) + 1
    let windows_size =
        (number_of_bins as f64 / NUM_WINDOWS_FOR_NORMALIZATION as f64).floor() as usize + 1;

    // norm_bins = np.array([])
    // let mut normalized_bins: Array1<f64> = Array1::zeros(0);

    // for i in range(0, len(bins_filled), win_size):
    //     win = bins_filled[i:i + win_size]
    //     if np.max(win) != 0:
    //         win = 50 * (win  / np.max(win))
    //     norm_bins = np.append(norm_bins, win)
    let binned_spec_start = BIN_SHIFT;
    let binned_spec_end = binned_spec_start + number_of_bins;
    for window_start in (binned_spec_start..binned_spec_end).step_by(windows_size) {
        let window_end = (window_start + windows_size).min(binned_spec_end);
        let mut window = binned_spectrum.slice_mut(s![window_start..window_end]);
        let window_max: f64 = window.iter().fold(f64::NEG_INFINITY, |acc, &x| acc.max(x));

        if window_max > 0.0 {
            window /= window_max;
            window *= 50.0;
        }
    }

    if use_flanking_peaks {
        for mz in mz.iter() {
            let index = (*mz / bin_size + bin_offset).floor() as usize + BIN_SHIFT - 1;
            let flanking_intensity = binned_spectrum[index] / 2.0;
            if index > 0 {
                binned_spectrum[index - 1] = flanking_intensity;
            }

            if index < binned_spectrum.len() - 1 {
                binned_spectrum[index + 1] = flanking_intensity;
            }
        }
    }

    // del bins_filled
    // drop(bins); // we're doing it in place

    // return norm_bins
    Ok(binned_spectrum)
}
