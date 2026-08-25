## Repo for seedling manuscript: Haigh et al., "Resistance and recovery of genetically and microbially diverse Douglas-fir (*Pseudotsuga menziesii*) seedlings to repeated experimental heat events."
Contains R code used for modeling temperature and duration response, statistical analyses, figure, and table generation. 

## Contents
### data
**Chlorophyll_conc**: Chlorophyll concentration data (not used in this analysis)

**LI6800**: Amax data included in supplemental

**Chlorophyll_fluorescence**: All data used in analysis
* instrument: Data downloaded directly from Opti-Sciences fluorometer
* raw: cleaned and compiled raw fluorescence data
* processed: .csv files generated from R scripts

### results
**figures**: .png and .pdf files of figures included in the manuscript

**tables**: .pdf and .docx files of tables included in the manuscript and supplemental

### scripts
**analysis**: Scripts containing model fitting, linear mixed effects models, or emmeans comparisons

**data_cleaning**: Scripts used to process raw data and generate .csv files used for analysis

**functions**: Scripts containing functions used to model temperature and duration response as well as log-linear recovery

**plotting**: Scripts containing code to plot figures and tables included in the manuscript
