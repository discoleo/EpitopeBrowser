# Package EpitopeBrowser

## Description
The main purpose of this package is to facilitate the exploration of HLA epitopes using a flexible shiny app. Epitopes derived from a protein sequence can be filtered and visualized.

**Steps:**
- Import epitopes from an external csv file;
- Browse and filter the epitopes;
- Compute the population coverage: based on frequent HLA-A, HLA-B, HLA-C alleles in Europe;

**Note:**\
The epitopes can be generated using https://iedb.org and exported as a csv-file - see the csv file in the folder /inst/examples. The application currently handles only HLA-I epitopes.

### Examples

Selecting the first 2 peptides on the **Epitopes** page ("EYHDVRVVL" and "FLEYHDVRV") and clicking the button **Display HLA** will display the HLA alleles for the 2 peptides and also compute the total population coverage (approx 71%).


## Authors

1. Author, Maintainer: L. Mada

