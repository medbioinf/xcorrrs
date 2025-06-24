use core::f64;

use ndarray::{s, Array1};

use crate::error::Error;

/// According to the original authors, after binning the experimental spectrum, the maximum intensity is normalized
/// over a number of fixed windows.
const NUM_WINDOWS_FOR_NORMALIZATION: u8 = 10;

/// Bins the theoretical spectrum based on the m/z values. Add bins for the shifts which avoids resizing the array later.
///
/// Arguments:
/// * `mz` - The m/z values of the theoretical spectrum.
/// * `bin_size` - The size of the bins to be used for binning the spectrum.
/// * `shift` - The number of bins to add to the start and end of the binned spectrum. (Adding the shift avoids resizing the array later.)
/// * `mz_max` - The maximum m/z value to consider for the spectrum. If `None`, the maximum m/z value from the input will be used.
///
pub fn theoretical_spectrum_binning(
    mz: &Array1<f64>,
    bin_size: f64,
    shift: usize,
    mz_max: Option<f64>,
) -> Result<Array1<f64>, Error> {
    if mz.is_empty() {
        return Err(Error::EmptyTheoreticalSpectrum);
    }

    // bins_filled = np.zeros(math.ceil(np.max(mz_array) / bin_width) + 1)
    let mz_max = match mz_max {
        Some(max) => max,
        None => *mz.last().unwrap(),
    };
    let number_of_bins = (mz_max / bin_size + 1.0).ceil() as usize;
    let mut bins: Array1<f64> = Array1::zeros(number_of_bins + 2 * shift);

    //  for mass in mz_array:
    for mz in mz.iter() {
        // index = int(mass // bin_width)
        let index = (*mz / bin_size).floor() as usize + shift;
        // bins_filled[index] = 50.0
        bins[index] = 50.0;
        // if index - 1 != -1:
        //     bins_filled[index - 1] = max(bins_filled[index - 1], 25.0)
        if index > 0 {
            bins[index - 1] = bins[index - 1].max(25.0);
        }
        bins[index + 1] = bins[index + 1].max(25.0);
    }

    Ok(bins)
}

/// Bins and normalizes the experimental spectrum.
///
/// Arguments:
/// * `mz` - The m/z values of the experimental spectrum.
/// * `intensities` - The intensity values of the experimental spectrum.
/// * `bin_size` - The size of the bins to be used for binning the spectrum.
/// * `shift` - The number of bins to add to the start and end of the binned spectrum. (Adding the shift avoids resizing the array later.)
///
pub fn experimental_spectrum_binning(
    mz: &Array1<f64>,
    intensities: &Array1<f64>,
    bin_size: f64,
    shift: usize,
) -> Result<Array1<f64>, Error> {
    if mz.len() != intensities.len() {
        return Err(Error::ExperimentalSpectrumShape(
            mz.len(),
            intensities.len(),
        ));
    }

    // bins_filled = np.zeros(math.ceil(np.max(mz_array) / bin_width) + 1)
    let mz_max = mz.last().unwrap();
    let number_of_bins = (mz_max / bin_size + 1.0).ceil() as usize;
    let mut binned_spectrum: Array1<f64> = Array1::zeros(number_of_bins + 2 * shift);

    // intensity_array = np.sqrt(intensity_array)
    // Doing this when adding to the bins, so we don't need to create a new array.

    // for mass, intensity in zip(mz_array, intensity_array):
    //     index = int(mass // bin_width)
    //     bins_filled[index] = max(bins_filled[index], intensity)
    for (mz, intensity) in mz.iter().zip(intensities.iter()) {
        let index = (*mz / bin_size).floor() as usize + shift;
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
    let binned_spec_start = shift;
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

    // del bins_filled
    // drop(bins); // we're doing it in place

    // return norm_bins
    Ok(binned_spectrum)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_theoretical_spectrum_binning() {
        let mz = Array1::from(vec![10.0, 20.0, 30.0, 40.0, 50.0, 60.0]);
        let bin_size = 5.0;
        let result = theoretical_spectrum_binning(&mz, bin_size, 2, None).unwrap();
        assert_eq!(result.len(), 13);
        assert_eq!(
            result,
            Array1::from(vec![
                0.0, 25.0, 50.0, 25.0, 50.0, 25.0, 50.0, 25.0, 50.0, 25.0, 50.0, 25.0, 50.0,
            ])
        );
    }

    #[test]
    fn test_theoretical_spectrum_binning_with_shift() {
        let shift: usize = 2;
        let mz = Array1::from(vec![10.0, 20.0, 30.0, 40.0, 50.0, 60.0]);
        let bin_size = 5.0;
        let result = theoretical_spectrum_binning(&mz, bin_size, shift, None).unwrap();
        assert_eq!(result.len(), 17);
        assert_eq!(
            result,
            Array1::from(vec![
                0.0, 0.0, 0.0, 25.0, 50.0, 25.0, 50.0, 25.0, 50.0, 25.0, 50.0, 25.0, 50.0, 25.0,
                50.0, 25.0, 0.0
            ])
        );
    }

    #[test]
    fn test_experimental_spectrum_binning() {
        let mz = Array1::from(vec![
            10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0,
            140.0, 150.0, 160.0, 170.0, 180.0, 190.0, 200.0,
        ]);
        let inten = Array1::from(vec![
            21.0, 20.0, 99.0, 50.0, 17.0, 79.0, 89.0, 35.0, 38.0, 85.0, 14.0, 19.0, 56.0, 39.0,
            39.0, 17.0, 98.0, 89.0, 73.0, 74.0,
        ]);
        let bin_size = 5.0;
        let result = experimental_spectrum_binning(&mz, &inten, bin_size, 0).unwrap();

        assert_eq!(result.len(), 41);
        assert_eq!(
            result,
            Array1::from(vec![
                0.0,
                0.0,
                50.0,
                0.0,
                48.79500364742667,
                0.0,
                50.0,
                0.0,
                35.53345272593507,
                0.0,
                21.852416110985086,
                0.0,
                47.10733619719444,
                0.0,
                50.0,
                0.0,
                47.98574349686966,
                0.0,
                50.0,
                0.0,
                50.0,
                0.0,
                20.291986247835695,
                0.0,
                23.63944858518838,
                0.0,
                50.0,
                0.0,
                41.72614801981401,
                0.0,
                31.542003094028026,
                0.0,
                20.824828195876073,
                0.0,
                50.0,
                0.0,
                50.0,
                0.0,
                45.28312928401491,
                0.0,
                50.0
            ])
        );
    }

    #[test]
    fn test_experimental_spectrum_binning_with_shift() {
        let mz = Array1::from(vec![
            10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0,
            140.0, 150.0, 160.0, 170.0, 180.0, 190.0, 200.0,
        ]);
        let inten = Array1::from(vec![
            21.0, 20.0, 99.0, 50.0, 17.0, 79.0, 89.0, 35.0, 38.0, 85.0, 14.0, 19.0, 56.0, 39.0,
            39.0, 17.0, 98.0, 89.0, 73.0, 74.0,
        ]);
        let bin_size = 5.0;
        let result = experimental_spectrum_binning(&mz, &inten, bin_size, 2).unwrap();

        assert_eq!(result.len(), 45);
        assert_eq!(
            result,
            Array1::from(vec![
                0.0,
                0.0,
                0.0,
                0.0,
                50.0,
                0.0,
                48.79500364742667,
                0.0,
                50.0,
                0.0,
                35.53345272593507,
                0.0,
                21.852416110985086,
                0.0,
                47.10733619719444,
                0.0,
                50.0,
                0.0,
                47.98574349686966,
                0.0,
                50.0,
                0.0,
                50.0,
                0.0,
                20.291986247835695,
                0.0,
                23.63944858518838,
                0.0,
                50.0,
                0.0,
                41.72614801981401,
                0.0,
                31.542003094028026,
                0.0,
                20.824828195876073,
                0.0,
                50.0,
                0.0,
                50.0,
                0.0,
                45.28312928401491,
                0.0,
                50.0,
                0.0,
                0.0,
            ])
        );
    }
}
