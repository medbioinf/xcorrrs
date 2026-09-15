use std::collections::{BTreeSet, HashMap};

use ordered_float::OrderedFloat;

#[derive(Clone, Debug)]
pub struct Precursor {
    pub charge: i32,
    /// m/z of the selected precursor ion, used to strip precursor peaks before scoring.
    pub mz: OrderedFloat<f64>,
}

#[derive(Clone, Debug)]
pub struct Spectrum {
    /// mz-values of the peaks * 10000 as i32, intensity as f64
    pub peaks: HashMap<i32, OrderedFloat<f64>>,
    /// mz-values of peaks as f64, sorted
    pub masses: BTreeSet<OrderedFloat<f64>>,
    pub precursor: Precursor,
}

pub trait Tolerance {
    fn calculate_lower_bound(&self, mass: OrderedFloat<f64>) -> OrderedFloat<f64>;
    fn calculate_upper_bound(&self, mass: OrderedFloat<f64>) -> OrderedFloat<f64>;
}

pub fn deisotope_spectrum<T: Tolerance>(
    current_spectrum: &mut Spectrum,
    tolerance_for_deisotoping: &T,
) {
    let mut peaks = current_spectrum.peaks.clone();
    let mut masses = current_spectrum.masses.clone();

    let cnt = current_spectrum.peaks.len();
    let mut new_dict_for_peaks: HashMap<i32, OrderedFloat<f64>> = HashMap::with_capacity(cnt);
    let mut temp_peaks: HashMap<i32, OrderedFloat<f64>> = HashMap::with_capacity(cnt);

    // Only looking isotoping distributions of doubly charged fragment ions
    let start = if current_spectrum.precursor.charge > 2 {
        2
    } else {
        1
    };

    for j in (1..=start).rev() {
        let mut temp_masses: BTreeSet<OrderedFloat<f64>> = BTreeSet::new();

        let mut last = OrderedFloat::<f64>(0.0);
        let mut last_intens = OrderedFloat::<f64>(0.0);
        let mut found_before = false;

        for &mass in &masses {
            // Get the peak mass as the *10000 i32 index to search in the peak list
            let r = get_mass_index(mass);
            let intensity = *peaks
                .get(&r)
                .expect("Mass index not found in peaks map during deisotoping");

            if last > OrderedFloat::<f64>(0.0) {
                let c13 = 1.0 / (j as f64);
                let lower = tolerance_for_deisotoping.calculate_lower_bound(mass);
                let upper = tolerance_for_deisotoping.calculate_upper_bound(mass);
                let mut found = false;

                for k in 1..3 {
                    let possible_isotope = last + (k as f64) * c13;
                    if possible_isotope > lower
                        && possible_isotope < upper
                        && last_intens > intensity
                    {
                        // Lies within tolerance
                        last_intens += intensity;
                        found = true;
                        found_before = true;
                        break;
                    } else if !found_before {
                        break;
                    }
                }

                if !found {
                    found_before = false;
                    let key = get_mass_index(last);
                    if let Some(existing) = temp_peaks.get_mut(&key) {
                        *existing += last_intens;
                    } else {
                        temp_peaks.insert(key, last_intens);
                        temp_masses.insert(last);
                    }

                    last = mass;
                    last_intens = intensity;
                }
            } else {
                last = mass;
                last_intens = intensity;
            }
        }

        let last_key = get_mass_index(last);
        temp_peaks.insert(last_key, last_intens);
        temp_masses.insert(last);

        // Set peaks and masses to the temp versions after deisotoping
        masses = temp_masses;
        peaks = temp_peaks;
        temp_peaks = new_dict_for_peaks;
        new_dict_for_peaks = HashMap::with_capacity(cnt);
    }

    // Update spectrum
    current_spectrum.peaks = peaks;
    current_spectrum.masses = masses;
}

pub fn get_mass_index(mass: OrderedFloat<f64>) -> i32 {
    (mass * 10000.0).into_inner() as i32
}
