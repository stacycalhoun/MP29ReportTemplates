Louisiana Vegetation Subroutine (ICM-LAVegMod)
================
Madeline R. Foster-Martinez, Eric D. White, Elizabeth R. Jarrell, Denise
J. Reed, Jenneke M. Visser

<style type="text/css">
h1 { /* Header 1 */
  color: #284D36;
}
h2 { /* Header 2 */
  color: #284D36;
}
h3 { /* Header 3 */
  color: #63AE7F;
}
</style>

<!--
&#10;The various sections of this template provide examples for certain types of formatting or inserting tables/figures. Remember to remove any template text that shouldn't be included in the report. Note: text following this <!-- are comments and do not show up in the final document rendering, only in the markdown code itself, but can be removed if preferred. If you want to see how the document output will look, click "Knit" from the script editor toolbar above. Also, if you have an Rstudio version v1.4 or v2022.02 and later, there is a Visual editor that you can use to see and adjust the formatting. The tab for the Visual editor is located below the script editor toolbar in the top left. 
&#10;-->

## Contents

<!--
For the table of contents, the links appear inside parentheses and cannot have uppercase letters, spaces or punctuation. Use a hyphen (-) in place of spaces and simply remove any punctuation. e.g. '1.0 Report Summary' has the link of (#10-report-summary). The text in square brackets is the text that will appear in the table of contents.
&#10;
Sections of the table of contents can be adjusted or removed as needed. Make sure sections removed or updated from the table of contents are also removed/updated in the report structure below.
&#10;-->

1.  [INTRODUCTION](#10-introduction)
2.  [ICM-LAVEGMOD GRID AND COVERAGE
    TYPES](#20-icm--lavegmod-grid-and-coverage-types)
3.  [INITIAL CONDITIONS](#30-initial-conditions)
    1.  [STANDARD PROCESSING METHOD](#31-standard-processing-method)
    2.  [SPECIAL AREAS](#32-special-areas)
        - [CANAL SPOIL BANKS](#canal-spoil-banks)
        - [FORCED DRAINAGE POLDERS](#forced-drainage-polders)
        - [DEVELOPED AND AGRICULTURAL
          LAND](#developed-and-agricultural-land)
4.  [INPUT VEGETATION ATTRIBUTES](#40-input-vegetation-attributes)
    1.  [PROBABILITY OF MORTALITY AND
        ESTABLISHMENT](#41-probability-of-mortality-and-establishment)
    2.  [DISPERSAL CLASS](#42-dispersal-class)
    3.  [HABITAT CLASS](#43-habitat-class)
    4.  [FFIBS SCORE](#44-ffbis-score)
5.  [HABITAT CLASS FUNCTIONS](#50-habitat-class-function)
    1.  [MORTALITY PROBABILITY](#51-mortality-probability)
        - [BOTTOMLAND HARDWOOD FOREST
          CLASS](#bottomland-hardwood-forest-class)
        - [SWAMP FOREST, EMERGENT WETLAND, FLOTANT
          CLASSES](#swamp-forest-emergent-wetland-flotant-classes)
        - [BARRIER ISLAND CLASS](#barrier-island-class)
    2.  [EXPANSION LIKELIHOOD](#52-expansion-likelihood)
        - [BOTTOMLAND HARDWOOD FOREST
          CLASS](#bottomland-hardwood-forest-class)
        - [SWAMP FOREST, EMERGENT WETLAND, FLOTANT
          CLASSES](#swamp-forest-emergent-wetland-flotant-classes)
        - [BARRIER ISLAND CLASS](#barrier-island-class)
6.  [ICM-LAVEGMOD ANNUAL PROCESSES](#60-icm--lavegmod-annual-processes)
    1.  [UPDATE WATER AREA](#61-update-water-area)
        - [LAND GAIN](#land-gain)
        - [VEGETATION ESTABLISHMENT ON LAND
          GAINED](#vegetation-establishment-on-land-gained)
        - [LAND LOSS](#land-loss)
    2.  [UPDATE VEGETATION COVERAGES](#62-update-vegetation-coverages)
        - [VEGETATION MORTALITY](#vegetation-mortality)
        - [VEGETATION ESTABLISHMENT](#vegetation-establishment)
        - [FLOTANT UPDATES](#flotant-updates)
        - [FLOTANT MORTALITY](#flotant-mortality)
        - [FLOTANT ESTABLISHMENT](#flotant-establishment)
    3.  [APPLY ACUTE SALINITY STRESS](#63-apply-acute-salinity-stress)
    4.  [ASSESS COVERAGES](#64-assess-coverages)
        - [CHECK MINIMUM COVERAGE](#check-minimum-coverage)
        - [CHECK TOTAL SUM](#check-total-sum)
        - [CALCULATE FFIBS SCORE](#calculate-ffibs-score)
        - [CALCULATE PERCENT VEGETATED
          LAND](#calculate-percent-vegetated-land)
        - [PREPARE OUTPUT](#prepare-output)
7.  [REFERENCES](#70-references)

[SUPPLEMENTAL MATERIALS](#supplemental-materials)

<!-- 
&#10;Sections below can be adjusted or removed as needed. Make sure section names match those in the table of contents so links will work.The number of # indicates the level of heading. There needs to be a space between the last # and the heading text for the output to show a heading.
&#10;-->

## 1.0 Introduction

Wetland vegetation is a critical driver of coastal processes. Vegetation
impacts hydrodynamics, altering flow paths and where sediment deposition
occurs; accretes organic matter, increasing the elevation of the
wetlands; and creates niche habitat, supporting fauna. Vegetation
species shift as environmental conditions change, and ICM-LAVegMod
predicts these vegetation shifts on an annual basis based on predictions
of changing conditions from other ICM subroutines (e.g., ICM-Hydro and
ICM-Morph). ICM-LAVegMod models the species distribution of 41
vegetation species that cover the full range of habitat types along the
Louisiana coast, from bottomland hardwood and swamp forest to saline
marsh. This report describes the rules and processes that govern how
species coverage expands and contracts in each model year in more
detail.

ICM-LAVegMod requires several input files to run, including:
spreadsheets containing the probability of establishment and mortality
of each vegetation species, attributes describing each vegetation
species, the initial vegetation conditions at the start of the model
run, and numerous environmental condition files representing the
hydrodynamics and landscape character across the model domain. There are
four main processes that occur each model year. First, the percentage of
water and land available for wetland vegetation establishment is updated
based on input from ICM-Morph. Second, the vegetation coverages are
updated based on the environmental conditions from ICM-Hydro. This
process occurs in two steps: vegetation mortality is assessed, removing
vegetation coverage and creating bareground for new vegetation
establishment; then vegetation establishment is assessed. For a
vegetation species to establish in a grid cell, it must have suitable
environmental conditions and be present in the dispersal zone of the
cell. The size of the dispersal zone depends on the vegetation species.
Third, acute salinity stress is applied, if it occurred in the cell
during the model year. Fourth, the coverages within each cell are
assessed and prepared for output.

ICM-LAVegMod outputs are used by ICM-Habitat Suitability Index (ICM-HSI)
to determine changes in availability and quality of habitat for
wildlife, fish, and shellfish species; by ICM-Morph to determine the
amount of organic accretion and to track the change in coverage type
(i.e., flotant, unvegetated, vegetated); and by ADCIRC+SWAN to determine
drag coefficient across the model domain. ICM-LAVegMod interacts
particularly closely with ICM-Morph. Within each ICM-LAVegMod grid cell,
ICM-Morph tracks the location and amount of land and open water, while
ICM-LAVegMod tracks the distribution (i.e., percent coverage) of each
vegetation species on the land portion of the cell.

The current ICM-LAVegMod subroutine is the product of interactive
development over many years. ICM-LAVegMod was part of the suite of
predictive models used to support the 2012 Coastal Master Plan (e.g.,
Visser et al., 2013) and was fully integrated into the ICM for the 2017
Coastal Master Plan (e.g., Visser & Duke-Sylvester, 2017). For each
master plan iteration, revisions and updates have been tested and
incorporated into the model. For comparison of the current subroutine
described here to previous models, refer to documentation from the
appropriate plans. For example, similar information for the version of
ICM-LAVegMod used for the 2017 Coastal Master Plan is provided in
Section 3.4.1 of [Attachment C3-22 to the 2017 Coastal Master
Plan](https://coastal.la.gov/wp-content/uploads/2017/04/Attachment-C3-22_FINAL_03.07.2017.pdf).

![](https://github.com/stacycalhoun/MP29ReportTemplates/blob/main/ICM-LAVegMod/C8_ICM-LAVegMod_Fig1.png)
<caption>

Figure 1. Spatial Resolution for ICM subroutine separate, overlapping
grids in the area around Marsh Island in Vermilion Bay.
</caption>

## 2.0 ICM-LAVegMod Grid and Coverage Types

The ICM-LAVegMod grid contains 480 m x 480 m orthogonal cells. The size
was set to align with the ICM-Morph grid, which has 30 m x 30 m
orthogonal pixels (Figure 1). The ICM-LAVegMod domain extends across all
ecoregions[^1] south of Interstate-10 (Figure 2) and contains 173,898
grid cells.

![](https://github.com/stacycalhoun/MP29ReportTemplates/blob/main/ICM-LAVegMod/C8_ICM-LAVegMod_Fig2.png)
<caption>

Figure 2. ICM-LAVegMod grid domain. White areas are outside of the ICM
ecoregions. Gray areas are NOTMOD. Shades of green are vegetation, and
blue is water.
</caption>

Each grid cell can contain 46 different coverage types: 41 vegetation
species; water; unvegetated wetland that was created this model year,
called *new bareground*; unvegetated wetland that was created the
previous model year, called *old bareground*; dead thick-mat flotant,
called *flotant bareground*; and upland or developed land, which is not
modeled and therefore called *NOTMOD*. The percent of each coverage type
within each grid cell is tracked over time at an annual timestep. Since
*NOTMOD* portions of each cell do not change over the course of the
model runs, if a grid cell is more than 95% *NOTMOD*, it is removed from
the grid and not assigned a grid cell ID value (example shown in Figure
1). Exceptions to this rule are described in Special Areas below.

![](https://github.com/stacycalhoun/MP29ReportTemplates/blob/main/ICM-LAVegMod/C8_ICM-LAVegMod_Fig3.png)
<caption>

Figure 3. ICM-LAVegMod grid in the area surrounding Cote Blanche. Color
legend is the same as in Figure 2. All grid cells (shown by white lines)
containing 95% or more of upland or developed land (NOTMOD) are removed
from the grid.
</caption>

## 3.0 Initial Conditions

Every cell within ICM-LAVegMod grid has initial coverage values of
vegetation and land types that are derived from the coastal Louisiana
landscape. The initial conditions for the 2023 Coastal Master Plan come
from the 2018 land use land cover (LULC) data collected by USGS
(Baustian et al., 2020). This data is processed satellite imagery and
has a 10 m resolution, which is higher resolution than the ICM-LAVegMod
grid with 480 m x 480 m cells. The USGS dataset contains more coverage
types than the 46 used in ICM-LAVegMod. The sections below outline the
standard method applied to produce the ICM-LAVegMod initial conditions
from the LULC dataset and special areas for which exceptions to this
method were applied. The full ICM-LAVegMod domain was evaluated,
producing an asc+ file of initial conditions. The format for this file
is a comma-separated table where each row is the cell ID and every
column is a coverage type. The initial conditions are generated once and
then used to initiate the model for all subsequent runs, including both
future with action (FWA) and future without action (FWOA).

### 3.1 Standard Processing Method

The following steps outline the standard processing method for
generating the ICM-LAVegMod initial conditions from the USGS LULC
dataset.

1.  *Reclassify the USGS LULC dataset to ICM-LAVegMod coverages*: The
    USGS LULC dataset contains more categories than are used in
    ICM-LAVegMod. The 10 m LULC cells were classified according to the
    conversions shown in Table S1 (in Supplemental Materials). For
    example, cells classified as water, palustrine aquatic bed, and
    estuarine aquatic bed in the USGS LULC dataset are all reclassified
    as water in the ICM-LAVegMod initial conditions. As shown under
    *NOTMOD* Table S1. ICM-LAVegMod does not model upland or developed
    land types. There are also 11 coverage types used by ICM-LAVegMod
    that are not included in the USGS LULC dataset (e.g., *new
    bareground* and *Sagittaria latifolia*). These coverages are set to
    zero in the initial conditions and are then allowed to increase
    based on modeled conditions as simulations run.

2.  *Tabulate the coverages within each ICM-LAVegMod cell*: The area of
    each ICM-LAVegMod coverage type within a cell was summed and
    tabulated, and then the total area for each coverage type was
    converted to a decimal percentage. For example, the ICM-LAVegMod
    cell shown in Figure 4. Figure 4 contains 1,443 LULC cells of water
    (144,300 m<sup>2</sup>), 854 LULC cells of *Spartina patens* (85,400
    m<sup>2</sup>), 5 LULC cells of *Spartina alterniflora* (500
    m<sup>2</sup>), and 2 LULC cells of *Schoenoplectus californicus*
    (200 m2). The initial conditions for this cell in ICM-LAVegMod are
    therefore 0.626 water, 0.371 *Spartina patens*, 0.002 *Spartina
    alterniflora*, 0.001 *Schoenoplectus californicus*, and 0 for all
    other coverage types.

![](https://github.com/stacycalhoun/MP29ReportTemplates/blob/main/ICM-LAVegMod/C8_ICM-LAVegMod_Fig4.png)
<caption>

Figure 4. Left: Example of an ICM-LAVegMod grid cell and the USGS LULC
dataset (colored squares). Right: Result of tabulating the coverage
types within the ICM-LAVegMod cell shown on the right. The percentages
shown in the pie chart are tabulated and used as the initial conditions.
</caption>

3.  *Determine cells to remove*: If an ICM-LAVegMod cell contained more
    than 95% of *NOTMOD*, meaning coverage categories that are not
    modeled in ICM-LAVegMod, then that cell was not included in model
    runs. The *NOTMOD* coverages in each cell do not change over the
    course of the simulation; therefore, removing these cells reduces
    memory requirements without altering the model results (see example
    shown in Figure 3).

4.  *Adjust for compatibility with ICM-Morph*: ICM-Morph operates on a
    30 m raster structure, and each pixel is initially assigned one of
    five coverage types: water, vegetated land, flotant, bareground, or
    NOTMOD. This assignment is done by resampling the 10 m USGS LULC
    dataset using a nearest neighbor algorithm in ArcGIS. Since the
    methods for creating the initial coverages differ between
    ICM-LAVegMod and ICM-Morph, the results also differ. To ensure
    compatibility, the percent coverages for ICM-LAVegMod are adjusted
    to match those of ICM-Morph. For example, the coverage of vegetated
    wetland in the ICM-Morph pixels within one ICM-LAVedMod cell matches
    the sum of all vegetated land in ICM-LAVegMod. These adjustments
    were minor (\< 5% within the cell). When vegetated land needed to be
    reduced, all species were reduced by the same proportion.

### 3.2 Special Areas

The following areas are special categories that require additional or
alternative processing steps to create the initial conditions. The
underlying data is the same USGS LULC dataset described above.

#### Canal Spoil Banks

When access canals were historically dredged through wetlands, the
dredged materials or “spoils” were commonly piled on the sides of the
canal, creating low linear levees or “banks.” Louisiana wetlands contain
an estimated 33,000 km of these dredged material levees (Turner &
McClenachan, 2018). Since these areas are elevated relative to
surrounding wetlands, they often contain upland vegetation species, but
over time these areas subside to elevations that can support wetland
species. Within the USGS LULC dataset, dredged material levees are often
classified as upland species that are not modeled in ICM-LAVegMod.
Classifying these areas as *NOTMOD* removes them from the ICM-LAVegMod
grid for modeling, even in future years when elevations could become
suitable for wetland species. To enable their inclusion, a different
classification conversion was used within the dredged material levee
areas. Pate (2014) identified and created polygons containing dredged
material levees along the Louisiana coast. Within these polygons, all
LULC coverages first assigned as *NOTMOD* (Table S1 in Supplemental
Materials, below) were reclassified to *old bareground* for the initial
conditions. This reclassification allows surrounding vegetation to
establish in these areas if conditions are right during model
simulations.

#### Forced Drainage Polders

Many communities and agricultural areas across the Louisiana coast are
drained with pumps and surrounded by low-lying levees that prevent tidal
flooding. These ‘polders’ were identified during the development of
input files for ICM-Morph (described below) and were subsequently
removed from the ICM-LAVegMod domain. These areas are cutoff from
natural hydrologic forcing and thus the processes driving vegetation and
morphology dynamics are not expected to align with patterns observed in
surrounding wetland areas. All coverages from the UGSS LULC dataset
within these areas were converted to *NOTMOD*.

#### Developed and Agricultural Land

Agricultural land, open fields, and other predominantly developed areas
in the upper-basins of the domain were often classified as bareland in
the USGS LULC dataset. In the standard method, these areas would become
*old bareground* in ICM-LAVegMod and thus be available for wetland
vegetation establishment in model simulations. While it is possible that
subsidence and sea-level rise might bring them under tidal influence,
and thus eligible for consideration in the ICM, it is likely that such
inundation in the future will be prevented by drainage or small levees
(i.e., they will become like polders). Thus, the decision was made to
not consider them in ICM-LAVegMod, and bareland LULC cells were
reclassified as *NOTMOD* in the areas listed below. This conversion left
a lot of grid cells with only *NOTMOD* and small coverages of water.
Cells with greater than 90% *NOTMOD* and greater than 95% of total
*NOTMOD* plus water were removed from the ICM-LAVegMod grid. This
classification and cell removal criteria was only applied to cells
meeting the conditions in the following areas: the entire Verret
Ecoregion, the entire Upper Verret Ecoregion, the entire Upper Barataria
Ecoregion, the portion of the Eastern Terrebonne Ecoregion that is north
of the Gulf Intracoastal Waterway (GIWW), the portion of the Western
Terrebonne Ecoregion that is north of the GIWW, and the portion of the
Maurepas Ecoregion south of Interstate-10[^2].

## 4.0 Input Vegetation Attributes

ICM-LAVegMod models the change in distribution for 41 vegetation species
in response to changing environmental conditions. Each species has
attributes used in the model that do not change over the course of the
simulation. These attributes are probability of mortality and
establishment tables, dispersal class, habitat class, and FFIBS score,
which are each described in the sections below. Table S2 (in
Supplemental Materials, below) contains the assigned attributes for each
species.

### 4.1 Probability of Mortality and Establishment

Each vegetation species has a relationship between environmental
conditions and the probability of mortality and establishment. These
relationships are given in input tables and were derived from vegetation
data from the Louisiana Coastwide Reference Monitoring System (CRMS).
These tables were updated for the 2023 Coastal Master Plan and can be
found in Baustian et al. (2020). If the given environmental condition
from ICM-Hydro falls between values in the table, the mortality and
establishment probabilities are linearly interpolated.

### 4.2 Dispersal Class

The dispersal class defines how far from its location a species can
establish. In the version of ICM-LAVegMod used to support the 2017
Coastal Master Plan modeling, species could establish within a cell if
they were present in one of the eight surrounding cells (i.e., Moore’s
neighborhood). For the 2023 Coastal Master Plan, ICM-LAVegMod has been
updated such that there are now three classes of dispersal: low, medium,
and high (Baustian et al., 2020). Low dispersal class species have the
same restrictions as for the 2017 Coastal Master Plan modeling and can
establish from the eight surrounding cells (light blue bounding box in
Figure 5). Medium dispersal class species can establish from two away
(darker blue bounding box in Figure 5). High dispersal class species do
not have to be present in the surrounding area to establish. These
species are considered to be “weedy” due to their high dispersal rates.

![](https://github.com/stacycalhoun/MP29ReportTemplates/blob/main/ICM-LAVegMod/C8_ICM-LAVegMod_Fig5.png)
<caption>

Figure 5. The extent of different dispersal classes: Low dispersal class
species can establish in the white cell from the cells bound by the
light blue border, and medium dispersal class species can establish in
the white cell from the cells bound by the darker blue border.
Background image is from the USGS LULC dataset. Black lines show the
ICM-LAVegMod grid.
</caption>

### 4.3 Habitat Class

Habitat classes are used to group vegetation species that are often
found in areas with similar environmental conditions and with similar
behaviors within the code. There are five habitat classes within the
ICM-LAVegMod code that correspond to generalized habitat types:
bottomland hardwood forest, swamp forest, emergent wetland, flotant, and
barrier island. The emergent wetland class includes species from
freshwater, intermediate, brackish, and saline marshes (often referred
to as FIBS or habitat type). The functions of the different habitat
classes are described in Section 5.0.

### 4.4 FFIBS Score

The FFIBS score is a numeric value that indicates the salinity regime of
the vegetation species. FFIBS stands for **F**orested, **F**resh,
**I**ntermediate, **B**rackish, and **S**aline. The values range from 0
for forested to 24 for saline. Species that exist in environments not
strictly in one category have scores that are between the different
cutoffs. For example, *Sagittaria lancifolia* and *Schoenoplectus
californicus* are both intermediate species, but *Schoenoplectus
californicus* has a FFIBS score of 2.75 because it can tolerate higher
salinity than *Sagittaria lancifolia*, which has a FFIBS score of 1.5.
Additional information on FFIBS score can be found in Baustian et
al. (2020) and Visser et al. (2002). Coverages that do not have a FFIBS
score are assigned a value of -9999. Calculation of FFIBS score for a
grid cell in ICM-LAVegMod is described in Section 2.6. The FFIBS score
of each cell is an input for ICM-Morph and is used to calculate the
organic matter accumulation rate (OMAR), which is described in ICM-Morph
Model Subroutines: Organic Matter Accumulation, below. Flotant species
do not have a FFIBS score and do not influence the OMAR on land.

## 5.0 Habitat Class Functions

There are five habitat classes within the ICM-LAVegMod code that
correspond to generalized habitat types: bottomland hardwood forest,
swamp forest, emergent wetland, flotant, and barrier island. The
emergent wetland class includes species from freshwater, intermediate,
brackish, and saline marshes. Each species and their assigned habitat
class are listed in Table 1. Using classes within Python allows for more
flexibility in the code. One function can be called multiple times but
will perform a different action depending on the class. The following
sections describe the two functions key to ICM-LAVegMod, mortality
probability and expansion likelihood, and how they differ for different
habitat classes.

<caption>

Table 1. Vegetation species included in ICM-LAVegMod, along with
corresponding habitat classification
</caption>

| **Habitat Class** | **Scientific Name** |
|:--:|----|
| Bottomland Hardwood Forest | *Quercus laurifolia, Quercus lyrate, Quercus nigra, Quercus texana, Quercus virginiana, Ulmus americana* |
| Swamp Forest | *Nyssa aquatica, Salix nigra, Taxodium distichum* |
| Flotant | *Eleocharis baldwinii, Panicum hemitomon* |
| Emergent Wetland | *Colocasia esculenta, Morella cerifera, Panicum hemitomon, Sagitttaria latifolia, Zizaniopsis miliacea, Cladium mariscus, Eleocharis cellulose, Iva frutescens, Typha domingensis, Paspalum vaginatum, Phragmites australis, Polygonum punctatum, Sagittaria lancifolia, Schoenoplectus californicus, Schoenoplectus americanus, Juncus roemerianus, Spartina alterniflora, Distichlis spicata, Avicennia germinans, Spartina patens, Spartina cynusuroides, Schoenoplectus robustus* |
| Barrier Island | *Uniola paniculate, Strophostyles helvola, Sporobolus virginicus, Spartina patens, Solidago sempervirens, Panicum amarum, Distichlis spicata, Baccharis halimifolia* |

There are three species that are in both the emergent wetland and
barrier island habitat classes: *Baccharis halimifolia*, *Distichlis
spicata*, and *Spartina patens*. These are kept distinct within the code
and follow the rules of the habitat class depending on location (i.e.,
when in a barrier island area, they follow the barrier island class).
Note, the results from the barrier island class have limited utility due
to the large grid size. The functions of this class depend on the
elevation of the grid cell, and the complex topography of the barrier
islands is smoothed out at this scale. A possible future improvement
would be to apply the barrier island class to the more finely resolved
grid used within the ICM-Barrier Island subroutine. While the barrier
island class species are included in ICM-LAVegMod following the
functions described below, the results are not used by other ICM
subroutines due to this limited utility.

### 5.1 Mortality Probability

The mortality probability determines the decrease in coverage of each
species. A value of 100% completely removes the species coverage from an
ICM-LAVegMod cell, and conversely, a value of 0% leaves the species
coverage unchanged. The mortality probability is a function of one or
more variables depending on the habitat class. The relationships between
environmental variable(s) and the mortality probability are given in the
mortality tables, an input to ICM-LAVegMod. The specific criteria for
calculating the mortality probability is given below for each habitat
class.

#### Bottomland Hardwood Forest Class

If the mean annual salinity is greater than 1.0 ppt, the mortality
probability is set to 100%. If the annual salinity is less than 1.0 ppt,
the mortality probability is a function of the elevation of the cell
(i.e., the height above mean annual water level).

#### Swamp Forest, Emergent Wetland, Flotant Classes

The mortality probability is a function of both mean annual salinity and
water level variability.

#### Barrier Island Class

The mortality probability is a function of the elevation of the cell.
Barrier island species can only exist in the designated barrier island
area. If the ICM-LAVegMod cell being evaluated is not in the barrier
island area, the mortality probability is set to 100%.

### 5.2 Expansion Likelihood

The expansion likelihood determines how the coverage of each species
expands. A value of 100% indicates the species has a high likelihood of
expansion, whereas a value of 0% means the species cannot expand. The
expansion likelihood is dependent upon the establishment probability,
which is a function of one or more environmental variables depending on
the habitat class. The relationships between environmental variable(s)
and the establishment probability are given in the establishment tables.

For most habitat classes, the expansion likelihood is also dependent on
the coverage of the species in the surrounding area, called the
dispersal coverage. The dispersal coverage is calculated as:

$$D_{i} = \frac{C_{T,i}}{A_{surrounding}} $$

<div align="right">

<caption>

Eq. 1
</caption>

</div>

Where $D_i$ is the dispersal coverage for the $i^{th}$ vegetation
species; $C_{T,i}$ is the total coverage for the $i^{th}$ vegetation
species in $A_{surrounding}$; and $A_{surrounding}$ is the total area of
the surrounding cells. If the species is in the low dispersal class,
$A_{surrounding}$ is the area of nine cells, and if the species is in
the medium or high dispersal classes, $A_{surrounding}$ is the area of
25 cells. The specific criteria for how establishment probability is
calculated for each habitat class described below.

### 5.3 Bottomland Hardwood Forest Class

For expansion to occur, the tree establishment condition must have been
met, meaning there must have been a 14 day period with no inundation
followed by a 14 day period with water depth less than 14 cm in the
current model year and the mean annual salinity must be less than 1.0
ppt. ICM-Hydro determines if these criteria are met for each
ICM-LAVegMod cell (from compartment-level water surface elevation and
grid-level average land surface elevation), and this information is
passed to ICM-LAVegMod each simulation year. If these criteria are not
met, then the establishment probability is set to 0% for bottomland
hardwood species. If these criteria are met in the cell being evaluated,
then expansion can occur, and the establishment probability is
determined. The establishment probability is a function of elevation of
the cell (i.e., the height above water). The establishment probability
is then multiplied by the dispersal coverage of the species to produce
the expansion likelihood.

### 5.4 Swamp Forest, Emergent Wetland, Flotant Classes

The establishment probability is a function of both mean annual salinity
and water level variability. The expansion likelihood is the product of
the establishment probability and the dispersal coverage.

### 5.5 Barrier Island Class

The establishment probability is a function of the elevation of the
cell. The dispersal coverage is not considered, meaning the
establishment probability is the same as the expansion likelihood. If
the cell is not in a set barrier area, the expansion likelihood is set
to 0%.

## 6.0 ICM LaVegMod Annual Process

ICM-LAVegMod updates vegetation coverages each year of the model
simulation. The following sections describe each processing step for
this update in the order it occurs. Each process applies to every
ICM-LAVegMod cell, but the processes are explained on an individual cell
basis. Inputs from ICM-Hydro and ICM-Morph are required. From ICM-Hydro,
there are five inputs: the annual mean salinity, water level
variability, height above annual mean water level, acute salinity
stress, and tree establishment condition. The annual mean salinity is
the annual average of the daily mean salinities calculated within the
open water area of each respective ICM-Hydro compartment. The water
level variability is the standard deviation of the hourly water
level[^3]. The acute salinity stress input indicates which ICM-LAVegMod
cells experienced salinity of greater than or equal to 5.5 ppt for
consecutive 2 weeks[^4]. Similarly, the tree establishment condition
indicates which ICM-LAVegMod cells met the conditions required for tree
establishment, described in Section 5.0. All of these inputs are at the
scale of the ICM-Hydro compartments. There is one input from ICM-Morph,
the percent water coverage in each ICM-LAVegMod cell. Due to the order
the subroutines of the ICM are run, the ICM-Hydro inputs are for the
same model year, and the ICM-Morph input is from the previous model
year. These annual inputs are in addition to the inputs attributes that
do not change throughout the simulation, described in Input Vegetation
Attributes section 2.4, above. The annual processes described here are
shown in the left side of Figure 6; the right side of Figure 6 shows an
example of how these processes govern the coverage changes of one cell.

![](https://github.com/stacycalhoun/MP29ReportTemplates/blob/main/ICM-LAVegMod/C8_ICM-LAVegMod_Fig6.png)
<caption>

Figure 6. Left: The four main processes performed each model year in
ICM-LAVegMod with details on how each process is governed. Right: An
example scenario in one grid cell. Vegetation abbreviations are as
follows: COES = Colocasia esculenta; SALA = Sagittaria lancifolia; PAHE2
= Panicum hemitomon.
</caption>

### 6.1 Update Water Area

Changes in water area that occurred in ICM-Morph are assessed and
incorporated into the ICM-LAVegMod coverages. ICM-Morph governs all
elevation changes and determines what area is eligible for vegetation
establishment, referred to as land, and what is open water. The
percentage of water in each ICM-LAVegMod cell is compared to the
percentage water in ICM-Morph, calculated at the end of the previous
model year. If the percentage is the same, land was not lost or gained,
and no further action is needed. If the percentage has increased or
decreased, land has been lost or gained,

#### Land Gain

If the percentage of water from ICM-Morph is less than the prior model
year, land was gained. The difference between the two percentage values,
$\Delta water$, is subtracted from the water coverage and added to the
*new bareground* coverage.

#### Vegetation Establishment on Land Gained

Vegetation species with high rates of dispersal are often the first to
establish on newly formed land. This step mimics that process by
allowing vegetation in the high dispersal class (i.e., “weedy” species)
to establish before other species on new land created during the
previous model year in ICM-Morph.

If new land was gained, vegetation species in the high dispersal class
(i.e., “weedy” species) are allowed to establish first. Establishment
probability is based on the mean annual salinity and water level
variability in the ICM-LAVegMod cell. If no species can establish due to
unfavorable conditions, then the area of the new land gained is added to
*old bareground*. This *old bareground* will be available for vegetation
establishment from all species during the Vegetation Establishment step.

#### Land Loss

If the percentage of water from ICM-Morph is greater than the year
prior, land was lost. The difference between the two percentage values,
$\Delta water$, is added to the *water coverage*. If any *old
bareground* exists in the ICM-LAVegMod cell, that coverage is lost
first. If there was no *old bareground* coverage or if it was less than
the amount of land lost, the vegetation coverages are all
proportionately decreased following Eq. 2. Flotant and *NOTMOD*
coverages, which are not altered by ICM-Morph, are not decreased and
remain the same.

$$C_{1,i} = C_{0,i} \left(\frac{A_\text{land, ICM-LAVegMod}}{A_\text{land, ICM-Morph}}\right)$$

<div align="right">

<caption>

Eq. 2
</caption>

</div>

$C$ is the coverage for the $i^{th}$ vegetation species; a subscript of
$1$ indicates the resulting value at the end of the step; a subscript of
$0$ indicates the value at start of the step;
$A_\text{land, ICM-LAVegMod}$ is the total vegetated land within the
cell from ICM-LAVegMod in the previous model year; and
$A_\text{land, ICM-Morph}$ is the total vegetated land within the cell
after updates from ICM-Morph in the previous model year.

### 6.2 Update Vegetation Coverages

#### Vegetation Mortality

This step assesses the amount of vegetation mortality that occurred
within each cell due to the environmental conditions in the current
model year. The mortality probability of each vegetation species is
calculated with the exception of flotant species, which are handled
separately (see Flotant Updates). The reduction in coverage of each
species is calculated based on the mortality tables, and the specific
criteria for each class is described in Section 5.0. The mortality
probability is directly applied to each species coverage. For example,
if the mortality probability is 20% based on the conditions associated
with an ICM-LAVegMod cell, the coverage of that species is reduced by
20% in that cell (Eq. 3). If a species coverage is reduced, the area it
no longer occupies is added to *new bareground* coverage (Eq. 4)

$$C_{1,i} = C_{0,i}*(1-M_i)$$

<div align="right">

<caption>

Eq. 3
</caption>

</div>

$$BG_{new} = \sum_{i=1}^n C_{0,i} * M_i $$

<div align="right">

<caption>

Eq. 4
</caption>

</div>

Where $C_i$ is the coverage for the $i^{th}$ vegetation species, $M_i$
is the mortality probability for the $i^{th}$ vegetation type;
$BG_{new}$ is *new bareground*; $n$ is the total number of vegetation
species; a subscript of $1$ indicates the resulting value at the end of
the step; and a subscript of $0$ indicates the value at start of the
step. At the end of this step, the *new* and *old bareground* coverages
are summed, and this area becomes the total area of unoccupied land,
which is available for vegetation establishment.

#### Vegetation Establishment

The total amount of unoccupied land (i.e., the sum of *new* and *old
bareground*) updated in the previous step is available for vegetation
establishment. In coastal wetlands, species must compete to establish on
the unoccupied land. Here, the establishment probability and the
dispersal coverage of the species in the surrounding ICM-LAVegMod cells
are used to calculate an expansion likelihood. Species with a higher
establishment probability and/or a greater coverage in the surrounding
areas have a larger expansion likelihood. The expansion likelihoods of
all species are then normalized to model the process of natural
competition. Two species with equal expansion likelihoods will occupy
equal portions of the available land.

In this step, the expansion likelihood of each vegetation species is
calculated as follows:

$$L_i = E_i * D_i$$

<div align="right">

<caption>

Eq. 5
</caption>

</div>

Where $L_i$ is the expansion likelihood for the $i^{th}$ vegetation
species; $E_i$ is the establishment probability for the $i^{th}$
vegetation type; and $D_i$ is dispersal coverage for the $i^{th}$
vegetation species. The establishment probability is calculated based on
the environmental conditions and comes from the establishment tables
(Baustian et al., 2020). The specific criteria of how expansion
likelihoods are calculated for each class is described in Section 5.0.
The resulting expansion likelihoods are normalized to divide the
unoccupied land accordingly, and the change in coverage for each species
follows:

$$C_{1, i} = C_{0, i} + \left(\frac{L_{i}}{\sum_{i = 1}^n L_i}\right)A_{unoccupied}$$

<div align="right">

<caption>

Eq. 6
</caption>

</div>

Where $C_i$ is the coverage for the $i^{th}$ vegetation species, $L_i$
is the expansion likelihood for the $i^{th}$ vegetation species;
$A_{unoccupied}$ is the total unoccupied area; $n$ is the total number
of vegetation species; a subscript of $1$ indicates the resulting value
at the end of the step; and a subscript of $0$ indicates the value at
start of the step.

If only one species has a non-zero expansion likelihood, it will occupy
all of the unoccupied land, and old and *new bareground* are reset to
0%. This step differs from the Vegetation Establishment on Land Gained
step in two ways: both *old* and *new bareground* are available for
vegetation establishment, as opposed to only *new bareground*, and all
vegetation species, as opposed to only the high dispersal class species,
have the opportunity to establish.

If all species have expansion likelihoods of zero, high dispersal
species are given the opportunity to establish without the restriction
of having dispersal coverage in the area. This modified form of
expansion likelihood is called spread likelihood and is equal to the
establishment probability (i.e., $E_i = S_i$). If there are more than
one non-zero spread likelihoods, the values are normalized to determine
the apportionment of the unoccupied area:

$$C_{1, i} = C_{0, i} + \left(\frac{S_{i}}{\sum_{i = 1}^m S_i}\right)A_{unoccupied}$$

<div align="right">

<caption>

Eq. 7
</caption>

</div>

Where $C_i$ is the coverage for the $i^{th}$ vegetation species, $S_i$
is the spread likelihood for the $i^{th}$ vegetation species;
$A_{unoccupied}$ is the total unoccupied area; $m$ is the total number
of vegetation species in the high dispersal class; a subscript of $1$
indicates the resulting value at the end of the step; and a subscript of
$0$ indicates the value at start of the step.

If both the expansion and spread likelihoods for all species are zero,
the unoccupied land remains unoccupied, and the coverage of *new
bareground* and *old bareground* are unchanged.

#### Flotant Updates

This step assesses the amount of flotant mortality that occurred within
each ICM-LAVegMod cell due to the environmental conditions in the
current model year. The mortality probability, or the percent reduction
in coverage, is calculated based on the mortality tables, which are a
function of annual mean salinity and water level variability. The
coverage lost by each species is tracked separately but is summed to
produce unoccupied flotant area.

#### Flotant Establishment

Flotant establishment can occur on any of the unoccupied flotant area
created in the previous step, in addition to *bareground flotant*
created in the previous model year. *Bareground flotant* is the area
occupied by thick-mat flotant that died in the previous model year. The
expansion likelihood for each species is calculated in the same way as
described in the Vegetation Establishment step above. If one or both of
the flotant species has a non-zero expansion likelihood, then all of the
unoccupied flotant and *bareground flotant* become established with
flotant vegetation.

If both flotant species have expansion likelihoods of zero, then the
following occurs: the *bareground flotant* and area of dead thin-mat
flotant are summed to become *dead flotant* coverage, and the area of
dead thick-mat flotant becomes the *bareground flotant* for the next
model year. The *dead flotant* coverage in each cell is passed to
ICM-Morph at the end of the ICM-LAVegMod run, where it is converted to
open water with a depth of 1 m. Within ICM-LAVegMod, the *dead flotant*
coverage is added to the water coverage. The *dead flotant* coverage is
reset to 0% at the start of every model year.

### 6.3 Apply Acute Salinity Stress

Freshwater wetland and flotant vegetation can experience mortality if
exposed to elevated salinities for a short period of time (i.e., on the
order of weeks). This impact of salinity may not be captured in the
annual mean salinity values. Since the previous vegetation mortality
steps depend on mean annual salinity, this acute salinity stress must be
assessed separately.

If the two-week average salinity in the ICM-LAVegMod cell exceeded 5.5
ppt at any time over the course of the ICM-Hydro model year, the cell is
flagged for acute salinity stress. The coverage of all freshwater marsh
species within the cell is converted to *new bareground*. The coverage
of thin-mat flotant is added to the *dead flotant* coverage and
converted to open water, and the coverage of thick-mat flotant is added
to the *bareground flotant*. Acute salinity stress is not applied to
bottomland hardwood or swamp forest species because they already have
low tolerance for salinity.

### 6.4 Assess Coverages

Once the previous steps are performed, the vegetation coverage updates
for the model year are complete. The following processing steps ensure
the model performed without error and prepare outputs for use in the
other subroutines of the ICM.

#### Check Minimum Coverage

As described in the mortality and establishment steps above, coverages
are changed by adding and subtracting percentages of the existing
coverage. The probability of mortality and tables, which are based on
CRMS data, have some overlap. This overlap means a species can be
completely removed in the mortality step and still have a non-zero
expansion likelihood, which can lead to small fractions of vegetation
coverage (e.g., $1 e^{-16}$) remaining. This issue can be problematic
when salinity regimes shift; due to the dependence on the presence of
vegetation in the expansion likelihood, a small presence of a species
can lead to unrealistic establishment.

To avoid this issue, fractions of vegetation below the minimum coverage
threshold of 1 m<sup>2</sup> are removed. This minimum threshold value
is $4.34𝑒^{−4}$ % of the 480 m x 480 m cells. If the coverage of any
species is less than 1 m<sup>2</sup>, the coverage is set to 0%.

#### Check Total Sum

All 46 coverages within each ICM-LAVegMod cell are summed, including the
non-vegetation type coverages (e.g., *NOTMOD* coverage). If the total
equals 100% $\pm$ 0.5%, it indicates the model performed as expected. If
not, a warning is raised alerting the user that an error occurred.

#### Calculate FFIBS Score

The FFIBS score of each cell is an indicator of the salinity regime of
the cell. FFIBS stands for *F*orested, *F*resh, *I*ntermediate,
*B*rackish, and *S*aline. For analyses using previous versions of
ICM-LAVegMod, the habitat type for each cell was determined by the
dominant vegetation species. This calculation has been improved for the
2023 Coastal Master Plan modeling by accounting for all vegetation
species within the cell. The FFIBS value of each species is weighted by
the area it occupies, providing a more accurate representation of the
salinity regime in the ICM-LAVegMod cell (Baustian et al., 2020). The
weighted FFIBS values are passed to ICM-Morph, where they are used to
calculate the OMAR in each cell, as described below.

The weighted FFIBS score of each cell is calculated as follows:

$$FFBIS_{weighted} = \frac{\sum_{i = 1}^n\left(FFIBS_i * C_i\right)}{\sum_{i=1}^n C_i}$$

<div align="right">

<caption>

Eq. 8
</caption>

</div>

Where $FFIBS_{weighted}$ is the weighted FFIBS score for the cell;
$FFIBS_i$ = the FFIBS score for the $i^{th}$ vegetation species; $C_i$
is the coverage for the $i^{th}$ vegetation species; and $n$ is the
total number of vegetation species. The flotant and barrier island
species are not included in this calculation. If no coverages have a
FFIBS score (e.g., the cell is 100% water or *old bareground*), then a
score of -9999 is assigned to the cell. The FFIBS score for each species
is listed in Table S2 in Supplemental Materials, below.

#### Calculate Percent Vegetated Land

The percent of each habitat type out of the total vegetation land in
each ICM-LAVegMod cell is calculated. These values are used as inputs to
the ICM-HSI subroutine. The habitat types include: bottomland hardwood
forest, swamp forest, fresh marsh, intermediate marsh, brackish marsh,
and saline marsh. Saline marsh includes barrier island species. The
percentages are calculated based on the total vegetated land in the
cell, not the total area. For example, if a cell only contains water and
swamp forest species, the percent land of swamp forest species would be
100%, and all other habitat types would be 0%.

#### Prepare Output

ICM-LAVegMod output for the model year is an asc+ file. The top portion
of this file is a table indicating the column and row location of each
cell ID in the ICM-LAVegMod grid. The bottom portion is a
comma-separated table where every row is a cell and each column is an
output value. There are 54 output values for each cell: the percent
coverage for each of the 41 vegetation species; the percent coverage for
the five non-vegetation coverage types, which are *water*, *old
bareground*, *new bareground*, *flotant bareground*, and *NOTMOD*; the
percent coverage of *dead flotant*; the percent of vegetated land that
is bottomland hardwood forest (pL_BF), swamp forest (pL_SF), fresh marsh
(pL_FM), intermediate marsh (pL_IM), brackish marsh (pL_BM), and saline
marsh (pL_SM); and the cell weighted FFIBS score. These outputs can then
be further processed to analyze and map different vegetation
distributions.

## 7.0 References

References can be copied and pasted here.

## Supplemental Materials

<caption>

Table S1. List of Coverage types (vegetation and land) within
ICM-LAVegMod and the corresponding coverages from the USGS LULC dataset
</caption>

| **Full or Scientific Name of ICM-LAVegMod Coverage Type** | **ICM-LAVegMod Symbol** | **USGS LULC Dataset Classes** |
|----|----|----|
| Water | WATER | Water; Palustrine Aquatic Bed; Estuarine Aquatic Bed |
| Not Modeled | NOTMOD | Developed, High Intensity; Developed, Medium Intensity; Developed, Low Intensity; Developed, Open Space; Cultivated Crops; Pasture/Hay; Grassland/Herbaceous; Upland - Mixed Deciduous Forest; Upland - Mixed Evergreen Forest; Upland Mixed Forest; Upland Scrub/Shrub |
| Old Bareground | BAREGRND_OLD | Unconsolidated Shore; Bare Land |
| New Bareground | BAREGRND_NEW | None |
| *Quercus laurifolia* | QULA3 | Cottonwood - willow mixing - bottomland hardwood sites - occasional flooding |
| *Quercus lyrate* | QULE | Lower site bottomland hardwoods such as overcup oak and water hickory |
| *Quercus nigra* | QUNI | Lower site bottomland hardwoods such as water oak - lower site ash |
| Quercus texana | QUTE | Sweetgum/nutall/willow oak - bottomland hardwoods seasonal flooding |
| *Quercus virginiana* | QUVI | Bottomland hardwoods/longleaf/slash pine mix infrequent flooding; Bottomland hardwoods/loblolly pine mix infrequent flooding; Sycamore/pecan/american elm - infrequently flooding; Sweetgum/yellow poplar; Swamp chestnut oak/cherrybark oak - bottomland hardwoods - infrequently flooding; River birch / sycamore - bottomland hardwood sites - infrequent flooding; Live oak / bottomland hardwoods mix |
| *Ulmus americana* | ULAM | Higher site bottomland hardwoods such as sugarberry/elm/greenash |
| *Nyssa aquatica* | NYAQ2 | Swamp tupelo dominant - Sweetbay mixing; Sweetbay dominant - swamp tupelo mixing; Tupelo dominant - cypress co-dom - low sites - freq. flooded |
| *Salix nigra* | SANI | Willow - low sites - wax myrtle mixing |
| *Taxodium distichum* | TADI2 | Red maple lowland; Cypress dominant - tupelo mixing - low sites - frequently flooded |
| *Eleocharis baldwinii* | ELBA2_Flt | ELBA2_FLT |
| *Panicum hemitomon* | PAHE2_Flt | PAHE2_Flt |
| Flotant Bareground | Bareground_flt | None |
| *Colocasia esculenta* | COES | COES |
| *Morella cerifera* | MOCE2 | MOCE2 |
| *Panicum hemitomon* | PAHE2 | PAHE2 |
| *Sagittaria latifolia* | SALA2 | None |
| *Zizaniopsis miliacea* | ZIMI | ZIMI |
| *Cladium mariscus* | CLMA10 | CLMA10 |
| *Eleocharis cellulose* | ELCE | ELCE |
| *Iva frutescens* | IVFR | Deciduous shrub scrub species - mixed; IVFR |
| *Paspalum vaginatum* | PAVA | PAVA |
| *Phragmites australis* | PHAU7 | PHAU7 |
| *Polygonum punctatum* | POPU5 | POPU5 |
| *Sagittaria lancifolia* | SALA | SALA |
| *Schoenoplectus californicus* | SCCA11 | SCCA11 |
| *Typha domingensis* | TYDO | TYDO |
| *Schoenoplectus americanus* | SCAM6 | SCAM6 |
| *Schoenoplectus robustus* | SCRO5 | SCRO5 |
| *Spartina cynusuroides* | SPCY | SPCY |
| *Spartina patens* | SPPA | SPPA |
| *Avicennia germinans* | AVGE | AVGE |
| *Distichlis spicata* | DISP | DISP |
| *Juncus roemerianus* | JURO | JURO |
| *Spartina alterniflora* | SPAL | SPAL |
| *Baccharis halimifolia* | BAHABI | None |
| *Distichlis spicata* | DISPBI | None |
| *Panicum amarum* | PAAM2 | None |
| *Solidago sempervirens* | SOSE | None |
| *Spartina patens* | SPPABI | None |
| *Sporobolus virginicus* | SPVI3 | None |
| *Strophostyles helvola* | STHE9 | None |
| *Uniola paniculate* | UNPA | None |

<caption>

Table S2. Vegetation species modeled in ICM-LAVegMod and their
attributes
</caption>

| **Scientific Name** | **Dispersal Class** | **Habitat Type** | **FFIBS score** |
|----|----|----|----|
| *Quercus laurifolia* | Low | Fresh | 0 |
| *Quercus lyrata* | Low | Fresh | 0 |
| *Quercus nigra* | Low | Fresh | 0 |
| *Quercus texana* | Low | Fresh | 0 |
| *Quercus virginiana* | Low | Fresh | 0 |
| *Ulmus americana* | Low | Fresh | 0 |
| *Nyssa aquatica* | Low | Fresh | 0 |
| *Salix nigra* | High | Fresh | 0 |
| *Taxodium distichum* | Low | Fresh | 0 |
| *Eleocharis baldwinii* | Medium | Fresh | -9999 |
| *Panicum hemitomon* | Low | Fresh | -9999 |
| *Colocasia esculenta* | High | Fresh | 0.25 |
| *Morella cerifera* | Medium | Fresh | 0.25 |
| *Panicum hemitomon* | Low | Fresh | 0.25 |
| *Sagitttaria latifolia* | High | Fresh | 0.25 |
| *Zizaniopsis miliacea* | High | Fresh | 0.25 |
| *Cladium mariscus* | Medium | Intermediate | 1.5 |
| *Eleocharis cellulosa* | Medium | Intermediate | 1.5 |
| *Iva frutescens* | Medium | Intermediate | 2.75 |
| *Paspalum vaginatum* | Medium | Intermediate | 2.75 |
| *Phragmites australis* | High | Intermediate | 2.75 |
| *Polygonum punctatum* | High | Intermediate | 1.5 |
| *Sagittaria lancifolia* | Medium | Intermediate | 1.5 |
| *Schoenoplectus californicus* | Medium | Intermediate | 2.75 |
| *Typha domingensis* | High | Intermediate | 2.75 |
| *Schoenoplectus americanus* | Medium | Brackish | 7.15 |
| *Schoenoplectus robustus* | Medium | Brackish | 11.5 |
| *Spartina cynusuroides* | Medium | Brackish | 11.5 |
| *Spartina patens* | Medium | Brackish | 7.15 |
| *Avicennia germinans* | High | Saline | 24 |
| *Distichlis spicata* | Medium | Saline | 17.5 |
| *Juncus roemerianus* | Medium | Saline | 17.5 |
| *Spartina alterniflora* | Low | Saline | 24 |
| *Uniola paniculate* | NA | NA | -9999 |
| *Strophostyles helvola* | NA | NA | -9999 |
| *Sporobolus virginicus* | NA | NA | -9999 |
| *Spartina patens* | NA | NA | -9999 |
| *Solidago sempervirens* | NA | NA | -9999 |
| *Panicum amarum* | NA | NA | -9999 |
| *Distichlis spicata* | NA | NA | -9999 |
| *Baccharis halimifolia* | NA | NA | -9999 |

[^1]: For a description and map of ecoregions used for the 2023 Coastal
    Master Plan, please refer to the Ecoregion and Regional Boundaries
    section of [Attachment C2: 50-Year FWOA Model Output, Regional
    Summaries –
    ICM](https://coastal.la.gov/wp-content/uploads/2024/02/C2_50YearFWOAModelOutputRegionalSummaries-ICM_v3.pdf#page=19).

[^2]: An ecoregion map is available in Figure 3 of Attachment C2:
    50-Year FWOA Model Output, Regional Summaries – ICM

[^3]: ICM-LAVegMod statistical relationships rely on standard deviation
    of hourly water level during the year; however, due to memory
    limitations, hourly water levels are not saved for every ICM-Hydro
    compartment for the entire year, however daily tidal range is saved.
    Using observed CRMS water level data, a linear correlation was fit
    to predict standard deviation of hourly water level from daily tidal
    range values. See [code on
    Github](https://github.com/CPRA-MP/ICM_Hydro/blob/e55c22d208b9322fc2cff1609ac87dfc46661fce/2D_ICM_summaries.f#L225)
    for relationship used.

[^4]: This acute salinity threshold value is set in ICM.py – a value of
    5.5 was used for 2023 modeling.
