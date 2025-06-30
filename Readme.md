# xcorrrs

Standalone Rust implementation of Comet's XCorr based on 

* > Eng JK, McCormack AL, Yates JR.   
  > An approach to correlate tandem mass spectral data of peptides with amino acid sequences in a protein database.   
  > J Am Soc Mass Spectrom.    
  > 1994;5(11):976-989.    
  > doi:10.1016/1044-0305(94)80016-2   

* > Eng JK, Fischer B, Grossmann J, Maccoss MJ.   
  > A fast SEQUEST cross correlation algorithm.   
  > J Proteome Res.    
  > 2008;7(10):4598-4602.    
  > doi:10.1021/pr800420s   

## Features

| feature | description |
| --- | --- |
| do-not-use-fast-xcorr | Enables the fast xcorr implementation. While it is functional it is not correctly implemented yet. |


## Tests
For each PSM in `test_files/LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_01.tsv` the matched peptide is converted to a Proforma String, including PSMs, and matched against the spectrum from the MS run using the xcorr reimplementation. The root means squared error is than calulated on Comet's reported xcorr and the results of the reimplementation (scaled by the highest Xcorr of both implementations). Accepted is RMSE below 0.02.

For development and debug purposes, the number of PSMs can be limited setting the environment variable: `TEST_NUMBER_OF_PSMS`

To test the fast xcorr add `--all-features` to the test command.

### Test prints
For the fast xcorr implementation prints can be generated for a visual inspection of the experimental and theoretical spectrum before and after preprocessin. The test case is skipped by default and can be executed by calling:
`cargo test --all-features -r --package xcorrrs --lib -- fast_xcorr::tests::test_prints --exact --show-output --include-ignored`

### Data
Test data is taken from [PXD028735](https://www.ebi.ac.uk/pride/archive/projects/PXD028735) ([LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_01.raw](https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/02/PXD028735/LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_01.raw)), converted to mzML using `msconvert` (3.0.24155, vendor peak picking) and identified with Comet 2025.01 using the protoeme UP000005640 (including isoforms, downloaded at 2024-05-08, compressed version at `test_files/2024-05-08_UP000005640_isoforms.fasta.gz`). Comet parameter file was mostly unchanged except for the MS and output paramters (`test_files/comet.params`).   
Due to the size of the mzML and GitHub file limits, the data array of the identified spectra were extracted and saved as parquet files without metadata to be used in tests.
