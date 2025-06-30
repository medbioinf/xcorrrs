use ndarray::Array1;
use rustyms::{Element::H as Hydrogen, Fragment, MassMode};

use crate::error::Error;

/// Converts mass to charge ration (Thompson) as Dalton
///
/// # Arguments
/// * `mz` - Mass to charge ratio (Thompson)
/// * `charge` - Charge
///
pub fn mass_to_charge_to_dalton(mz: f64, charge: usize) -> f64 {
    let charge = charge as f64;
    mz * charge - Hydrogen.mass(None).unwrap().value * charge
}

/// Creates a theoretical spectrum from a list of fragments.
///
/// # Arguments
/// * `fragments` - A slice of `Fragment` objects representing the theoretical fragments.
/// * `max_charge` - The maximum charge state to consider for the fragments.
/// * `max_mz` - The maximum m/z value to consider for the fragments.
///
pub fn create_threoretical_spectrum(
    fragments: &[Fragment],
    max_charge: usize,
    max_mz: f64,
) -> Result<Array1<f64>, Error> {
    let mut mz: Vec<f64> = fragments
        .iter()
        .filter(|f| f.charge.value <= max_charge)
        .filter_map(|f| f.mz(MassMode::Monoisotopic).map(|mz| mz.value))
        .filter(|mz| *mz <= max_mz)
        .collect();

    mz.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    Ok(Array1::from(mz))
}

#[cfg(test)]
pub mod tests {
    use super::*;
    use std::path::PathBuf;

    use polars::{frame::DataFrame, prelude::*};

    #[test]
    fn test_mass_to_charge_ratio_to_dalton() {
        assert_eq!(
            mass_to_charge_to_dalton(464.888129195412, 3),
            1391.640912490542
        )
    }

    /// Reads the test data from a TSV file and processes it to include a new column for proforma peptides.
    /// TEST_NUMBER_OF_PSMS environment variable limits to the X best scored PSMs.
    ///
    pub fn read_test_data() -> DataFrame {
        // Read the Comet results from a TSV file
        let mut comet_df = CsvReadOptions::default()
            .with_has_header(true)
            .with_parse_options(
                CsvParseOptions::default()
                    .with_separator(b'\t')
                    .with_comment_prefix(Some("#")),
            )
            .try_into_reader_with_file_path(Some(PathBuf::from(
                "test_files/LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_01.tsv",
            )))
            .unwrap()
            .finish()
            .unwrap();

        // Sort by xcorr in descending order
        comet_df
            .sort_in_place(
                ["xcorr"],
                SortMultipleOptions::default().with_order_descending(true),
            )
            .unwrap();

        // Reduce the number of PSMs if the environment variable is set
        let mut comet_df = match std::env::var("TEST_NUMBER_OF_PSMS") {
            Ok(number_of_psms) => comet_df.slice(0, number_of_psms.parse::<usize>().unwrap()),
            Err(_) => comet_df,
        };

        // Create the proforma peptides column
        let modified_peptide = comet_df.column("modified_peptide").unwrap().str().unwrap();
        let profoma_peptides = modified_peptide
            .iter()
            .filter_map(|s| {
                s.map(|s| {
                    let mut proform_string = s[2..s.len() - 2]
                        .to_string()
                        .replace("[15.9949]", "[+15.9949]");
                    if proform_string.contains("C") {
                        proform_string = proform_string.replace("C", "C[+57.02146]");
                    }
                    proform_string
                })
            })
            .collect::<Vec<String>>();

        comet_df
            .with_column(
                Series::new("proforma_peptide".into(), profoma_peptides)
                    .cast(&DataType::String)
                    .unwrap(),
            )
            .unwrap();

        comet_df
    }

    /// Reads a spectrum from a Parquet file.
    ///
    /// # Arguments
    /// * `scan` - The scan number to read the spectrum for.
    ///
    pub fn get_spectrum(scan: &str) -> (Array1<f64>, Array1<f64>) {
        let spec_df = ParquetReader::new(
            std::fs::File::open(format!("test_files/spectra/scan_{scan}.parquet")).unwrap(),
        )
        .read_parallel(ParallelStrategy::None)
        .finish()
        .unwrap();

        let mz_array = spec_df["mz"]
            .f64()
            .unwrap()
            .to_ndarray()
            .unwrap()
            .to_owned();

        let intensity_array = spec_df["intensity"]
            .f64()
            .unwrap()
            .to_ndarray()
            .unwrap()
            .to_owned();

        (mz_array, intensity_array)
    }

    /// Extracts spectru, data arrays from the mzML file and saves them as Parquet files without any further metadata
    /// Important do not set TEST_NUMBER_OF_PSMS environment variable, as this will limit the spectrum extraction to a certain number of scans.
    ///
    #[test]
    #[ignore = "Spectrum extration."]
    fn spectrum_extraction() {
        let comet_df = read_test_data();

        let mut mzml_byte_reader = std::io::BufReader::new(
            std::fs::File::open("test_files/LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_01.mzML")
                .unwrap(),
        );
        let mut mzml = dihardts_omicstools::proteomics::io::mzml::reader::Reader::read_indexed(
            &mut mzml_byte_reader,
            None,
            true,
            false,
        )
        .unwrap();

        for i in 0..comet_df.height() {
            let scan = comet_df["scan"].i64().unwrap().get(i).unwrap();

            let binary_data_array_list = mzml
                .get_spectrum(&format!("controllerType=0 controllerNumber=1 scan={scan}"))
                .unwrap()
                .binary_data_array_list;

            let mz_array = binary_data_array_list
                .get_mz_array()
                .unwrap()
                .deflate_data()
                .unwrap();
            let intensity_array = binary_data_array_list
                .get_intensity_array()
                .unwrap()
                .deflate_data()
                .unwrap();

            let mut spec_frame = DataFrame::new(vec![
                Column::new("mz".into(), mz_array),
                Column::new("intensity".into(), intensity_array),
            ])
            .unwrap();

            let writer = ParquetWriter::new(
                std::fs::File::create(format!("test_files/spectra/scan_{scan}.parquet")).unwrap(),
            )
            .with_compression(ParquetCompression::Zstd(Some(
                ZstdLevel::try_new(22).unwrap(),
            )));
            writer.finish(&mut spec_frame).unwrap();
        }
    }
}
