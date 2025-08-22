use ndarray::Array1;
use rustyms::{
    system::{e, usize::Charge},
    CompoundPeptidoformIon,
    Element::H as Hydrogen,
    FragmentationModel, MassMode,
};

use crate::error::Error;

/// Da to m/z conversion
///
/// # Arguments
/// * `mass` - Mass in Dalton
/// * `charge` - Charge state
///
pub fn dalton_to_mass_to_charge(mass: f64, charge: usize) -> f64 {
    let charge = charge as f64;
    (mass + Hydrogen.mass(None).unwrap().value * charge) / charge
}

/// Creates a theoretical spectrum from a list of fragments.
///
/// # Arguments
/// * `peptide` - The peptide for which to generate the theoretical spectrum.
/// * `fragmentation_model` - The fragmentation model to use for generating fragments.
/// * `max_charge` - The maximum charge state to consider for the fragments.
///
pub fn create_threoretical_spectrum(
    peptide: &CompoundPeptidoformIon,
    fragmentation_model: &FragmentationModel,
    max_charge: usize,
) -> Result<Array1<f64>, Error> {
    let mut mz: Vec<f64> = peptide
        .generate_theoretical_fragments(Charge::new::<e>(max_charge), fragmentation_model)
        .into_iter()
        .filter_map(|f| f.mz(MassMode::Monoisotopic).map(|mz| mz.value))
        .collect();

    mz.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    Ok(Array1::from(mz))
}

#[cfg(test)]
pub mod tests {
    use crate::configuration::{Configuration, FinalizedConfiguration};

    use super::*;
    use std::path::PathBuf;

    use polars::{frame::DataFrame, prelude::*};

    /// Just a sanity check to make sure that max charge is uses as is in rustyms
    #[test]
    fn test_fragment_creation() {
        let max_charge = 6;
        let config: FinalizedConfiguration = Configuration::default().into();
        let peptide = CompoundPeptidoformIon::pro_forma("DIGSETK", None).unwrap();

        let mut charge_states = vec![false; max_charge + 1];
        charge_states[0] = true; // Charge 0 is not valid, but we use it for indexing

        peptide
            .generate_theoretical_fragments(
                Charge::new::<e>(max_charge),
                &config.fragmentation_model,
            )
            .iter()
            .for_each(|f| charge_states[f.charge.value] = true);

        assert!(charge_states.iter().all(|x| *x));
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

    /// Reads the experimental spectrum Eng provided for the peptide DIGSETK.
    ///
    pub fn get_eng_experimental_spectrum() -> (Array1<f64>, Array1<f64>) {
        let spec_df =
            ParquetReader::new(std::fs::File::open("test_files/eng/DIGSETK.parquet").unwrap())
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

    /// Reads the partial xcorr spectrum and the index provided by Eng.
    /// Return is a tuple with the index and the fast xcorr values.
    ///
    pub fn get_eng_fast_xcorr_spectrum() -> (Array1<usize>, Array1<f64>) {
        let spec_df = CsvReadOptions::default()
            .with_has_header(true)
            .with_rechunk(true)
            .with_parse_options(
                CsvParseOptions::default()
                    .with_separator(b'\t')
                    .with_comment_prefix(Some("#")),
            )
            .try_into_reader_with_file_path(Some(PathBuf::from(
                "test_files/eng/DIGSETK.process.tsv",
            )))
            .unwrap()
            .finish()
            .unwrap();

        (
            spec_df
                .column("index")
                .unwrap()
                .i64()
                .unwrap()
                .to_ndarray()
                .unwrap()
                .mapv(|x| x as usize),
            spec_df
                .column("fast_xcorr")
                .unwrap()
                .f64()
                .unwrap()
                .to_ndarray()
                .unwrap()
                .to_owned(),
        )
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
