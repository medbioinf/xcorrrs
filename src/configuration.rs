use rustyms::model::FragmentationModel;

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

pub struct Configuration {
    pub fragmentation_model: FragmentationModel,
    pub bin_size: f64,
    pub theoretical_fragment_ions: bool,
    pub minimum_intensity: f64,
    pub max_fragment_charge: usize,
}

impl Configuration {
    pub fn new(
        bin_size: f64,
        fragmentation: FragmentationMethod,
        theoretical_fragment_ions: bool,
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

        Self {
            bin_size,
            fragmentation_model,
            theoretical_fragment_ions,
            max_fragment_charge,
            minimum_intensity,
        }
    }
}

impl Default for Configuration {
    fn default() -> Self {
        Self::new(0.02, FragmentationMethod::CID, true, false, 3, 0.0)
    }
}
