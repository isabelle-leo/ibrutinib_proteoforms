#Guidelines for the R project

Anaylsis results are split by conceptual sections of the paper. These can be reproduced from PRIDE accession files - cell line data and pulldown data: PXD047187, clinical CLL data: PXD028936. (If you need a password, you're here early! Check your reviewer manuscript, or reach out to us.) 
Annotations (uniprot, interpro, and others) are possible to download from uniprot or biomart's service, but we also provide a script to fetch them directly in R, annotate_hits.R. We performed this and reproduced our results most recently on March 5, 2024.

The proteoforms eSet and NPARC results are generated in proteoform_analysis.R. Preprocessing of the CLL clinical data is performed in PSMS_to_Proteins_CLL.R. All additional scripts with the tag "cll" in the title rely on this data.

Several scripts (*analysis.R) perform functions on the main objects. The scripts read in a proteins object and their output is stored in respective subfolders of output. This can be tables, figures or modified proteins objects.

Accessory functions are stored in "functions" and sourced in analysis scripts. Some specific functions for plotting are inside the analysis scripts directly.

Meta data objects are provided in the repository for convenience, these files may not be complete. All raw data and processed data tables are provided in our publication (DOI: to be announced)

#To run the workflow

Copy raw data to ./data, e.g. target_psms.txt (Nextflow output).

Copy meta data files to ./meta.

Process proteoforms, run NPARC, run CLL aggregation, run annotate_hits

*analysis.R

Result paths are not set by default and should be be user defined
