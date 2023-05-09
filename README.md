# A multiOMICs view of multiple sclerosis white matter lesions

The complex nature of the sustained inflammation of multiple sclerosis (MS) lesions in the white matter (WM) of the brain have not been resolved to date. Using several high throughput molecular profiling technologies, snRNASeq, snATACSeq and spatial transcriptomics, we seek a better understanding of the complex processes driving the compartementalized and chronic inflammation of the tissue.

The data generated for this project has been deposited to GEO and will be available on acceptance of the manuscript with a peer-reviewed journal.


## Relevant GEO records

Newly sequenced records are available under the GEO accession number GEO???. Tissue RNASeq used in a prior study has been reused. This data is available under the GEO accession number GSE138614.

This study incorporates multiple different data modalities.
- White matter brain lesion tissue spatial transcriptomics (8 slides)
- snRNASeq of microdisected white matter lesion tissue from different lesion types and several donors and controls.
- snATACSeq of microdisected white matter lesion tissue from different lesions and several donors and controls.
- previously published tissue (bulk) RNASeq of MS white matter lesions based on 10 MS patients and 5 controls (GSE138614).


## Code
The code in this repository documents most computational analyses done for this study. Most of these analyses have been performed using the [UCloud](https://escience.sdu.dk) infrastructure provided by the University of Southern Denmark. UCloud provides a virtualized computing environment, which allows users to run on high performance compute infrastructure using dockerized applications.

The code is roughly organised into directories depending on the data type. An effort has been made to integrate the different data modalities. Analyses that only depend on one of the data types are organised in the corresponding directory under the respective name. Analyses depending on one, or more data types can be found in the folder 'integrated'

In addition to the R code provided in this manuscript, [SCANet](https://pypi.org/project/scanet/) has been to identify regulatory modules in the data.

The analysis has been performed using ```R version 4.2.1 (2022-06-23)``` using a variety of packages which we provide as an explicit environment file in this repository.



## Publication

For more information please refer to our preprint:
```
Maria Elkjaer, Anne Hartebrodt, Mhaned Oubounyt et al. A single-cell multi-omics map of cell-type-specific mechanistic drivers of multiple sclerosis lesions, 08 May 2023, PREPRINT (Version 1) available at Research Square [https://doi.org/10.21203/rs.3.rs-2845466/v1]
```
