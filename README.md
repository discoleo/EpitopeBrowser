# Package EpitopeBrowser

## Description
The main purpose of this package is to facilitate the exploration of HLA epitopes using a flexible shiny app. Epitopes derived from a protein sequence can be filtered and visualized.

**Steps:**
- Import epitopes from an external csv file;
- Browse and filter the epitopes;
- Compute the population coverage: based on frequent HLA-A, HLA-B, HLA-C alleles in Europe;

**Note:**\
The epitopes can be generated using https://iedb.org and exported as a csv-file - see the example csv file in the folder /inst/examples. The application handles both HLA-I and HLA-II epitopes, although the data sets have to be processed individually.

### Examples

**Epitopes** tab: Selecting the first 2 peptides ("EYHDVRVVL" and "FLEYHDVRV") and clicking the button **Display HLA** will display the HLA alleles for the 2 peptides and also compute the total population coverage (approx 71%).

**Note:**
- Population coverage is currently compute using a haploid model (although it is easy to compute for the diploid model as well);
- Population coverage computed using the diploid model may overestimate the actual population coverage (due to linkage disequilibria);


## Authors

1. Author, Maintainer: L. Mada

