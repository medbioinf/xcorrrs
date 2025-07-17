use rustyms::model::{FragmentationModel, PrimaryIonSeries, SatelliteIonSeries};

pub enum FragmentationMethod {
    CID,
    HCD,
    ETD,
}

#[allow(clippy::from_over_into)]
impl Into<FragmentationModel> for FragmentationMethod {
    fn into(self) -> FragmentationModel {
        match self {
            FragmentationMethod::CID => FragmentationModel::cid_hcd().clone(),
            FragmentationMethod::HCD => FragmentationModel::cid_hcd().clone(),
            FragmentationMethod::ETD => FragmentationModel::etd().clone(),
        }
    }
}

#[derive(PartialEq)]
pub enum Ion {
    A,
    B,
    C,
    D,
    V,
    W,
    X,
    Y,
    Z,
}

pub struct Configuration {
    pub fragmentation_model: FragmentationModel,
    /// https://uwpr.github.io/Comet/parameters/parameters_202501/fragment_bin_tol.html
    pub bin_size: f64,
    /// https://uwpr.github.io/Comet/parameters/parameters_202501/fragment_bin_offset.html
    pub bin_offset: f64,
    /// https://uwpr.github.io/Comet/parameters/parameters_202501/fragment_bin_offset.html
    pub use_flanking_peaks: bool,
    pub minimum_intensity: f64,
    pub max_fragment_charge: usize,
}

impl Configuration {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        bin_size: f64,
        bin_offset: f64,
        fragmentation: FragmentationMethod,
        use_flanking_peaks: bool,
        ions: Vec<Ion>,
        use_neutral_loss_ions: bool,
        max_fragment_charge: usize,
        minimum_intensity: f64,
    ) -> Self {
        let mut fragmentation_model: FragmentationModel = fragmentation.into();

        if !use_neutral_loss_ions {
            let mut a_ions = fragmentation_model.a.clone();
            a_ions.neutral_losses = Vec::new();
            fragmentation_model = fragmentation_model.a(a_ions);

            let mut b_ions = fragmentation_model.b.clone();
            b_ions.neutral_losses = Vec::new();
            fragmentation_model = fragmentation_model.b(b_ions);

            let mut c_ions = fragmentation_model.c.clone();
            c_ions.neutral_losses = Vec::new();
            fragmentation_model = fragmentation_model.c(c_ions);

            let mut x_ions = fragmentation_model.x.clone();
            x_ions.neutral_losses = Vec::new();
            fragmentation_model = fragmentation_model.x(x_ions);

            let mut y_ions = fragmentation_model.y.clone();
            y_ions.neutral_losses = Vec::new();
            fragmentation_model = fragmentation_model.y(y_ions);

            let mut z_ions = fragmentation_model.z.clone();
            z_ions.neutral_losses = Vec::new();
            fragmentation_model = fragmentation_model.z(z_ions);
        }

        if !ions.contains(&Ion::A) {
            fragmentation_model = fragmentation_model.a(PrimaryIonSeries::none());
        }

        if !ions.contains(&Ion::B) {
            fragmentation_model = fragmentation_model.b(PrimaryIonSeries::none());
        }

        if !ions.contains(&Ion::C) {
            fragmentation_model = fragmentation_model.c(PrimaryIonSeries::none());
        }

        if !ions.contains(&Ion::D) {
            fragmentation_model = fragmentation_model.d(SatelliteIonSeries::default());
        }

        if !ions.contains(&Ion::V) {
            fragmentation_model = fragmentation_model.v(SatelliteIonSeries::default());
        }

        if !ions.contains(&Ion::W) {
            fragmentation_model = fragmentation_model.w(SatelliteIonSeries::default());
        }

        if !ions.contains(&Ion::X) {
            fragmentation_model = fragmentation_model.x(PrimaryIonSeries::none());
        }

        if !ions.contains(&Ion::Y) {
            fragmentation_model = fragmentation_model.y(PrimaryIonSeries::none());
        }

        if !ions.contains(&Ion::Z) {
            fragmentation_model = fragmentation_model.z(PrimaryIonSeries::none());
        }

        Self {
            bin_size,
            bin_offset,
            fragmentation_model,
            use_flanking_peaks,
            max_fragment_charge,
            minimum_intensity,
        }
    }
}

impl Default for Configuration {
    fn default() -> Self {
        Self::new(
            0.02,
            0.0,
            FragmentationMethod::CID,
            true,
            vec![Ion::B, Ion::Y],
            false,
            3,
            0.0,
        )
    }
}
