# EF-models4freu
Effect factors for freshwater eutrophication
This R script implements a model simulating the impact of nitrogen (N) emissions from human activities on freshwater fish biodiversity, as described in Zhou et al. (2023).
The model incorporates regionalized species sensitivity distributions (SSDs) relating freshwater fish species richness to N concentrations across 367 ecoregions and 48 global combinations of realms and major habitat types. It also provides effect factors (EFs) for life cycle assessment (LCA) at a 0.5° × 0.5° spatial resolution.

The script comprises the following modules:

1. compile_occurrence_records_time_JZ.R
Collects, cleans, and merges freshwater fish occurrence records, adapted from the "occ2range4fish" scripts by vbarbarossa (https://github.com/vbarbarossa/occ2range4fish.git). This script specifically extracts temporal information from occurrence data.

2. compile_pairdata_N.R
Matches species occurrence records with nitrogen concentrations simulated by IMAGE-GNM at a 0.5° × 0.5° resolution on an annual basis.

3. SSD_fish_N.R
Models the statistical relationship between nitrogen concentrations and the potentially disappeared fraction (PDF) of fish species. It also generates SSD plots for 367 ecoregions and 48 combinations of realms and major habitat types.

4. EF_fish_N.R
Calculates both average and marginal effect factors (EFs) for nitrogen emissions at a global 0.5° × 0.5° resolution.



citation: Zhou J, Mogollón JM, van Bodegom PM, Barbarossa V, Beusen AHW, Scherer L. Effects of Nitrogen Emissions on Fish Species Richness across the World's Freshwater Ecoregions. Environ Sci Technol. 2023 Jun 6;57(22):8347-8354. doi: 10.1021/acs.est.2c09333. Epub 2023 May 22. PMID: 37216582; PMCID: PMC10249400.