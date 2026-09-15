use std::ops::Deref;

use mzannotate::annotation::model::FragmentationModel;
use serde::{Deserialize, Serialize};

/// Ion types used in configuration to select which ion series to use
/// for the fragmentation model.
///
pub use crate::utils::Ion;

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
    /// Same as: https://uwpr.github.io/Comet/parameters/parameters_202502/max_fragment_charge.html
    pub max_fragment_charge: usize,
    /// Ions to use in the fragmentation model.
    pub ions: Vec<Ion>,
    /// Same as: https://uwpr.github.io/Comet/parameters/parameters_202502/use_NL_ions.html
    pub use_neutral_loss_ions: bool,
    /// Same as: https://uwpr.github.io/Comet/parameters/parameters_202502/clear_mz_range.html
    pub clear_mz_range: Option<(f64, f64)>,
}

impl Default for Configuration {
    /// Default is for high-resolution MS2 data (e.g. Orbitrap) and CID/HCD
    fn default() -> Self {
        Self {
            bin_size: 0.02,
            bin_offset: 0.0,
            use_flanking_peaks: true,
            minimum_intensity: 0.0,
            max_fragment_charge: 3,
            ions: vec![Ion::B, Ion::Y],
            use_neutral_loss_ions: false,
            clear_mz_range: None,
        }
    }
}

/// Finalized configuration for the xcorr algorithm.
/// This configuration is created from the `Configuration` struct and contains
/// a `FragmentationModel` that is built based on the selected ions.
///
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FinalizedConfiguration {
    /// Fragmentation model that is built based on the selected ions.
    pub fragmentation_model: FragmentationModel,
    /// Original configuration that was used to create this finalized configuration.
    configuration: Configuration,
}

impl Deref for FinalizedConfiguration {
    type Target = Configuration;

    fn deref(&self) -> &Self::Target {
        &self.configuration
    }
}

impl From<Configuration> for FinalizedConfiguration {
    fn from(configuration: Configuration) -> Self {
        let fragmentation_model = crate::utils::build_fragmentation_model(
            &configuration.ions,
            configuration.use_neutral_loss_ions,
        );

        Self {
            fragmentation_model,
            configuration,
        }
    }
}

impl Default for FinalizedConfiguration {
    fn default() -> Self {
        Configuration::default().into()
    }
}
