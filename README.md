**Local-scale dispersal constraints promote spatial structure and arthropod diversity within a tropical sky-island**
#
The repository contains the pipeline to perform bioinformatics tools from metabarcoding data of the project *Local-scale dispersal constraints promote spatial structure and arthropod diversity within a tropical sky-island*. We focus on a single tropical sky-island, Nevado de Toluca of the Transmexican Volcanic Belt, where we sampled whole-communities of arthropods for eight orders with a comparable design at a spatial scale ranging from 50 m to 20 km, using 840 pitfall traps and whole community metabarcoding. These samples were then used to build metabarcoding libraries of arthropods using COI marker. 
#
**Primers and overhang adapter:**

**COI** expected size 418 pb

B_F 5' CCIGAYATRGCITTYCCICG 3' (Shokralla et al., 2015)

Fol-degen-R 5’ TANACYTCNGGRTGNCCRAARAAYCA 3' (Yu et al., 2012)

The overhang adapter sequence must be added to the locus‐specific primer
for the region to be targeted. The Illumina overhang adapter sequences to be
added to locus‐specific sequences are:

Forward overhang: 5’ TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG‐[locus‐
specific sequence]
Reverse overhang: 5’ GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG‐[locus‐
specific sequence]
#

## Repository organization
Contains data and scripts for the sections *Bioinformatics processing to identify OTUs at different thresholds of genetic similarity*, *Community diversity and composition* and *Similarity distance-decay and landscape connectivity* of the manuscript.

The **bin** directory cointains the scripts used in this pipeline. All scripts run using `bin` as working directory. The **data** directory is not included in this repository, but data is available at this [dryad](https://XXXXXXX).

Data comes from Cornell Institute of Biotechnology, Cornell University, USA, for sequencing on a lane of Illumina MiSeq 2x300 bp.
#

### `/bin/`

The scripts in `/bin` should be run in the order they are numbered. R functions used by some of these scripts are not numbered and have the extension `.R`. Html notebooks are provided for some of the analyses in R.

Scripts content:

* `0.install_allprograms.txt` well, not actually a script, but this fille contains many sugestions for installs packages for the first step that include the processing included quality checking, primer removal, pair merging, quality filtering, denoising, and clustering each library independently, into Arribas et al., 2020. The packages are *fastqc*, *fastx-toolkit*, *trimmomatic-0.36*, *pairfq-0.17*, *usearch-9.2*, and *usearch-10*. (EN PROCESO)
* `0.Soup_processing_steps.sh` performs a filtering haplotypes analysis of each sample. (EN PROCESO)
* `0.Steps_afterProcessing.txt`well, not actually a script. Contains the settings used to run afters steps for each library.  Here, realizar pasos dentro de las carpetas de resultados de la script anterior.  (EN PROCESO)
* `1.PrepareRastersNT.r` Script to reclassify altitude, slope, and vegetation types output rasters to desired values and to to create a flat lanscape.
* `2.PlotRaster.R` Plots each resistance among with sampling points.
* `2.Circuitscape_estimatingResDist.txt` well, not actually a script. Contains the settings used to run Circuitscape for each of the rasters.
* `3.to_get_conservative_threhold_on_original_all.r` This script gets conservative thershold on origibale table OTUs/ZOTUs by artropods order. It is important to place the correct number of columns and rows, because this can change in each group. In this step, we can join all the groups or put them separately as in this case.
* `4.to_get_Diversity_All.r`This script gets all scripts (source) of diversity using 8 arthropods order at multi-hierarchical levels.
* `5.Plot_Diversity_All.r`This script plot "Community diversity and composition", "plot global Richness by sites", "Beta diversity", "Non-Metric Multidimensional scaling (NMDS) ordinations of community similarity", and "Accumulation Curves".
* `6.to_get_BetaDiversity_DistanceDecay_IBD_IBR.r`This script gets all scripts (source) to estimate "Distance decay at the multihierchical leveles" and Isolation by resistence.
* `7.Plot_DistanceDecay_IBD_IBR`This script plot "Distance decay" and Isolation by resistance at the figure 5, figure 6 and Figure S6.
#

These scripts use the data in `genetic`, `spatial` and `meta`.

### `/genetic`

Contains genetic *data in* and *data out* for each order

Genetic *data in* corresponds to `0.Soup_processing_steps.sh` output using the subset of each order. 

Genetic *data out* corresponds to `4.to_get_Diversity_All.r`using the subset of 8 orders and `6.to_get_BetaDiversity_DistanceDecay_IBD_IBR.r` output using the subset of Diptera and Collembola order. 
#

### `/spatial`

Contais spatial data as follows:

* `/Elevation	` contains `.asc` rasters for the Nevado de Toluca. `NevTol_Alt.tif` is the original dataset, the rest are the result of reclassifying (output of `1.PrepareRastersNT` script) it for each of the altitudinal resistance surfaces. 

* `/Slope	` contains `.tif` rasters for the Nevado de Toluca slope. `NevTol_Pen.tif` is the original dataset, the rest are the result of reclassifying (output of `1.PrepareRastersNT` script) it for each of the slope resistance surfaces. 

* `/VegetationType	` contains `.tif` rasters for the Nevado de Toluca using Vegetation type. `nevado_f.tif` is the original dataset, the rest are the result of reclassifying (output of `1.PrepareRastersNT` script) it for each of the vegetation type resistance surfaces. 

* `surveyed_mountainNevadoToluca` contains sampling points (`.csv`) used in `2.PlotRasters.R`

* `IBDistanceMatrix` contains geomatrix focal point (`*.txt`) used in to run `6.to_get_BetaDiversity_DistanceDecay_IBD_IBR.r`.

* `IBResistanceFlatMatrix` flat output (`*.txt`) from Circuit scape used in to run `6.to_get_BetaDiversity_DistanceDecay_IBD_IBR.r`.

* `IBResistanceSlopeMatrix` slope output (`*.txt`) from Circuit scape used in to run `6.to_get_BetaDiversity_DistanceDecay_IBD_IBR.r`.

* `IBResistanceVegTypeMatrix` Vegetation type output (`*.txt`) from Circuit scape used in to run `6.to_get_BetaDiversity_DistanceDecay_IBD_IBR.r`.

* `IBResitanceElevationMatrix` elevation output (`*.txt`) from Circuit scape used in to run `6.to_get_BetaDiversity_DistanceDecay_IBD_IBR.r`.

* `Circuitscape` contains the focal points (`*_focalpoints.txt`) used to run Circuit scape and the output (`/out`). (EN CONTRUCCION).

### `/meta`
`ConservationForestNevadoToluca.csv` contains metadata for each of the samples sequenced in a lane Miseq. Each column names refer to:

* `ID`: sample number ID 
* `Code`: ID of the sequencing run 
 * `Label_Metabarcoding`: sample name of each library: e.g. `CON_NTO_TLC_31TCONS1`: `CON` Conservation (Treatment), `NTO`Nevado de Toluca (Mountain), `TLC`Tlacotepec (locality), `31TCONS1` ID of the sequencing sample run.    
 * `Locality`: Locality of the sampling 
* `key`: Abbreviation of the sampling location
* `Municipality`: Municipality of the sampling
* `State`: State of the sampling
* `Natural_Protected_Area`: Natural Protected Area of the sampling
* `key`: Abbreviation of the Natural Protected Area
* `Latitude`: Latitude of the sampling 
* `Longitude`: Longitude of the sampling 
* `Samplig_Altitude`: Meter above level seal of each the sampling point.
* `Treatment`: Name of the treatment 
* `Key`: Abbreviation of the treatment
* `Season`: Name of season in the sampling
* `Forest_type`: Name of the tree specie in the forest.

**END**