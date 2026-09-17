use std::ops::Deref;

use serde::{Deserialize, Serialize};

/// Ion types used in configuration to select which ion series to use
/// for the fragmentation model.
///
pub use crate::utils::Ion;
use crate::utils::{FragmentIonConfiguration, HasFragmentIonConfiguration};

/// Configuration for the xcorr algorithm.
/// The configuration is similar to Comet's configuration parameters.
/// Although it is much simpler as it only contains parameters related to the
/// fragmentation model and the experimental spectrum processing.
///
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Configuration {
    /// Same as: https://uwpr.github.io/Comet/parameters/parameters_202501/fragment_bin_tol.html
    pub bin_size: f64,
    /// Same as: https://uwpr.github.io/Comet/parameters/parameters_202501/fragment_bin_offset.html
    pub bin_offset: f64,
    /// Same as: https://uwpr.github.io/Comet/parameters/parameters_202502/theoretical_fragment_ions.html
    pub use_flanking_peaks: bool,
    /// Same as: https://uwpr.github.io/Comet/parameters/parameters_202502/minimum_fragment_intensity.html
    pub minimum_intensity: f64,
    /// Fragment-ion selection shared across scoring algorithms: max fragment charge, ion
    /// series (`max_fragment_charge`/`ions`/`use_neutral_loss_ions`, same as Comet's
    /// https://uwpr.github.io/Comet/parameters/parameters_202502/max_fragment_charge.html and
    /// https://uwpr.github.io/Comet/parameters/parameters_202502/use_NL_ions.html).
    #[serde(flatten)]
    pub common: FragmentIonConfiguration,
    /// Same as: https://uwpr.github.io/Comet/parameters/parameters_202502/clear_mz_range.html
    pub clear_mz_range: Option<(f64, f64)>,
}

impl Deref for Configuration {
    type Target = FragmentIonConfiguration;

    fn deref(&self) -> &Self::Target {
        &self.common
    }
}

impl HasFragmentIonConfiguration for Configuration {
    fn fragment_ion_configuration(&self) -> &FragmentIonConfiguration {
        &self.common
    }
}

impl Default for Configuration {
    /// Default is for high-resolution MS2 data (e.g. Orbitrap) and CID/HCD
    fn default() -> Self {
        Self {
            bin_size: 0.02,
            bin_offset: 0.0,
            use_flanking_peaks: true,
            minimum_intensity: 0.0,
            common: FragmentIonConfiguration::default(),
            clear_mz_range: None,
        }
    }
}

/// Finalized configuration for the xcorr algorithm.
/// This configuration is created from the `Configuration` struct and contains
/// a `FragmentationModel` that is built based on the selected ions.
///
pub type FinalizedConfiguration = crate::utils::Finalized<Configuration>;
