# Exotic species traits across levels of establishment
---

Scripts and data for reproducing the results obtained by Bernery, Marino & Bellard (2023) in the paper "Relative importance of exotic species traits in determining invasiveness across levels of establishment: Example of freshwater fish."

## Description of the Data and file structure

### 1. Traits and establishment groups of exotic fish

The data file Fish_exo_database.rds contains the taxonomy, traits, and establishment information of the 222 exotic fish (222 lines x 24 columns).
Lines refer to species.

Column description:
- Species (character): binomial name of species as in FishBase
- EhBd (numeric): Ratio between the eye position and body depth
- BlBd (numeric): Ratio of body length to body depth. Reflects the hydrodynamism of a species                 
- PFiBd (numeric): Ratio of pectoral fin position relative to body depth.
- MoBd (numeric): Ratio between the oral gape position and body depth
- Length (numeric): Maximum body length of the species (in cm)
- Freshwater (logical): TRUE is fthe species is from freshwater
- occ.introduced (numeric): number of basins were the species has at least one established exotic population
- intro.inside.native.regions (numeric): number of basins were the species has at least one established exotic population inside its native bioregion
- intro.outside.native.regions (numeric): number of basins were the species has at least one established exotic population outside its native bioregion
- Inside (binary): species has an exotic population established inside its native bioregion (1) or not (0)
- Outside (binary): species has an exotic population established outside its native bioregion (1) or not (0)
- InsideOrOutside (character): location of establishment (OnlyInside, OnlyOutside, InsideAndOutside)
- detritus (binary): species eats detritus (1) or not (0)
- nekton (binary): species eats nekton (1) or not (0)
- plants (binary): species eats plants (1) or not (0)
- zoobenthos (binary): species eats zoobenthos (1) or not (0)
- zooplankton (binary): species eats zooplankton (1) or not (0)
- Nb_diet (numeric): number of different diets consumed by the species
- RepGuild1 (character): categories of parental care (non-guarder, guarder, and bearer)
- Area.Bassins (numeric): total area of the native basins (in km²)
- MaxBio5 (numeric): 95th percentile of maximum temperature of warmest month in the native region (in degrees °C)
- MinBio6 (numeric): 5th percentile of minimal temperature of coldest month in the native region (in degrees °C)
- Amplitudetemp (numeric) : Range of temperature amplitude of the species’ native basins


### 2. How to use the scripts

The first script computes the PCoA and derives the axes of the multidimensional functional space based on the traits of all exotic fish. The second script performs the statistical analyses based on the functional space axes and creates Figure 2. The third script performs the hypervolume comparison and makes Figure 3. The last script performs the traits-by-trait analysis.

## Sharing/Access Information

Links to other publicly accessible locations of the data: https://zenodo.org/record/7982955

Trait data came from FishBase (Froese & Pauly 2019) and FISHMORPH (Brosse et al. 2021), and exotic information came from Tedesco et al. (2017) and Leroy et al. (2019). 

Brosse, S., Charpin, N., Su, G., Toussaint, A., Herrera-R, G. A., Tedesco, P. A., & Villéger, S. (2021). FISHMORPH: A global database on morphological traits of freshwater fishes. Global Ecology and Biogeography, 30(12), 2330–2336. 

Froese, R., & Pauly, D. (2019). FishBase. World Wide Web electronic publication.

Leroy, B., Dias, M. S., Giraud, E., Hugueny, B., Leprieur, F., Oberdorff, T., & Tedesco, P. A. (2019). Global biogeographical regions of freshwater fish species. Journal of Biogeography, 46(11), 2407–2419. https://doi.org/doi.org/10.1111/jbi.13674

Tedesco, P. A., Beauchard, O., Bigorne, R., Blanchet, S., Buisson, L., Conti, L., Cornu, J.-F., Dias, M. S., Grenouillet, G., Hugueny, B., & others. (2017). A global database on freshwater fish species occurrence in drainage basins. Scientific Data, 4(1), 1–6.

