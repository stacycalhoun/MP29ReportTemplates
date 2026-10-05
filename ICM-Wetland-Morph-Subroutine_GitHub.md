Wetland Morphology Subroutine
================
Madeline R. Foster-Martinez, Eric D. White, Elizabeth R. Jarrell, Denise
J. Reed, Jenneke M. Visser

## Contents

1.  [INTRODUCTION](#10-introduction)
2.  [SOFTWARE ENVIRONMENT](#20-software-environment)
    1.  [FILE FORMATS](#21-file-formats)
3.  [INITIAL CONDITIONS AND INPUT
    FILES](#30-initial-conditions-and-input-files)
    1.  [TOPOBATHYMETRIC DIGITAL ELEVATION MODEL
        (DEM)](#31-topobathymetric-digital-elevation-model-dem)
    2.  [LANDSCAPE
        COMPOSITION/LANDTYPE](#32-landscape-compositionlandtype)
    3.  [POLDERS](#33-polders)
    4.  [DEEP SUBSIDENCE](#34-deep-subsidence)
    5.  [SHALLOW SUBSIDENCE](#35-shallow-subsidence)
    6.  [MARSH EDGE EROSION RATES](#36-marsh-edge-erosion-rates)
    7.  [ORGANIC MATTER ACCUMULATION
        RATES](#37-organic-matter-accumulation-rates)
    8.  [ACTIVE DELTAIC ZONES](#38-active-deltaic-zones)
    9.  [MAINTAINED NAVIGATION CHANNELS AND
        RIVERS](#39-maintained-navigation-channels-and-rivers)
    10. [MODEL GRID LOOKUP MAPS](#310-model-grid-lookup-maps)
    11. [INPUT PARAMETERS CONTROL
        FILE](#311-input-parameters-control-file)
4.  [ICM-MORPH MODEL SUBROUTINES](#40-icm-morph-model-subroutines)
    1.  [MAIN CONTROL PROGRAM AND PARAMETER
        SETTINGS](#41-main-control-program-and-parameter-settings)
    2.  [PRE-PROCESSING INPUT FILES](#42-pre-processing-input-files)
    3.  [INUNDATION DEPTHS](#43-inundation-depths)
    4.  [MARSH EDGE DELINEATION](#44-marsh-edge-delineation)
    5.  [MINERAL SEDIMENT DEPOSITION](#45-mineral-sediment-deposition)
        - [EROSION OF WATER BODY
          BOTTOMS](#erosion-of-water-body-bottoms)
    6.  [ORGANIC MATTER ACCUMULATION](#46-organic-matter-accumulation)
    7.  [FLOTANT MARSH LOSS](#47-flotant-marsh-loss)
    8.  [IMPLEMENT SHORELINE PROTECTION
        PROJECTS](#48-implement-shoreline-protection-projects)
    9.  [MARSH EDGE EROSION](#49-marsh-edge-erosion)
    10. [MAP SPATIAL DISTRIBUTION OF
        BAREGROUND](#410-map-spatial-distribution-of-bareground)
    11. [CHRONIC INUNDATION STRESS
        THRESHOLDS](#411-chronic-inundation-stress-thresholds)
    12. [UPDATE ELEVATION](#412-update-elevation)
    13. [CLASSIFY NEW SUBAERIAL LAND](#413-classify-new-subaerial-land)
    14. [UPDATE LANDTYPE
        CLASSIFICATION](#414-update-landtype-classification)
    15. [CALCULATE INUNDATION FOR USE IN HSI
        EQUATIONS](#415-calculate-inundation-for-use-in-hsi-equations)
    16. [IMPLEMENT MARSH CREATION AND LANDBRIDGE
        PROJECTS](#416-implement-marsh-creation-and-landbridge-projects)
    17. [IMPLEMENT RIDGE RESTORATION AND LEVEE
        PROJECTS](#417-implement-ridge-restoration-and-levee-projects)
    18. [POST-PROCESSING/SUMMARY
        OUTPUTS](#418-post-processingsummary-outputs)
5.  [MODELING SUBMERGED AQUATIC VEGETATION
    (SAV)](#50-modeling-submerged-aquatic-vegetation-sav)
6.  [REFERENCES](#60-references)

## 1.0 Introduction

ICM-Morph is a relative elevation model of wetland morphology developed
specifically for use in coastal Louisiana. ICM-Morph originated in the
Wetland Morphology model developed for the 2012 Coastal Master Plan and
was subsequently integrated with other models to form the ICM for the
2017 Coastal Master Plan (White et al., 2017).

ICM-Morph updates the elevation of wetlands as a function of subsidence,
erosion and deposition of water bottoms, mineral sediment deposition and
organic matter accretion of vegetated wetlands, marsh edge erosion.
Changes to inundation patterns are then modeled as a result of relative
sea level rise (as modeled by the elevation change model and assumed
rates of eustatic sea level rise). Inundation stress of vegetated marsh
areas may lead to chronic stress and ultimately loss of vegetation and
eventual loss of land to open water. These dynamics will then result in
changing coastal hydrology (as modeled by ICM-Hydro), which in turn
impact vegetation coverage (as modeled by ICM-LAVegMod); vegetation
dynamics subsequently impact organic accretion rates, and the hydrologic
changes impact mineral sediment deposition. Thus, the parsimonious
relative elevation change model is the foundation for the morphological
aspects of the ICM and all other modeling tools used for the development
of the Louisiana Coastal Master Plan.

## 2.0 Software Environment

Previous versions of the wetland morphology models used for the 2012 and
2017 master plans were largely built using proprietary geoprocessing
software tools in the ESRI ArcGIS platform. In order to make the ICM
completely platform-independent, the model algorithms were converted to
Fortran, a compiled software language that is frequently used in high
performance computing (HPC) systems. Both the ICM-Hydro and ICM-BI
subroutines were already coded in Fortran; therefore, converting the
ICM-Morph code to Fortran required no additional compilers or expertise
that was not already needed for other ICM subroutines.

Once converted to Fortran, there were two distinct advantages as
compared to the old ESRI-dependent versions of the model. First, Fortran
programs are able to be compiled in both Windows and Linux computing
environments; ICM-Morph has been tested, and shown to work, in both of
these environments. Production runs for the 2023 plan were performed on
Bridges-2, a Linux-based HPC hosted by the Pittsburgh Supercomputing
Center. The second advantage of the Fortran ICM-Morph is the ability to
utilize large memory arrays and binary data files. In previous
ESRI-dependent versions of ICM-Morph, every calculation step required a
raster file to be written to disc. While many of these temporary
calculation rasters were not permanently stored, the reading and writing
of files to disc is a potential performance bottleneck. Even more so,
when moving to HPC environments, which are generally optimized for
processors to utilize large memory arrays, and not around file
read/write speeds.

### 2.1 File Formats

The file types used by ICM-Morph generally consist of three different
types. First, raster based files are initially read into the program as
ASCI XYZ text files where the first two columns are the X and Y
coordinates (in UTM Zone 15N coordinates) and the third column (Z)
contains the value for that pixel for whatever dataset is being
represented in the raster. This could be elevation, landscape
composition, subsidence rates, etc. A full accounting of the file types
used are provided in the following section. Once these files are read
into the program in XYZ format, ICM-Morph will save all outputs (and
convert the input files at the end of the first year) into binary arrays
that are an order—of-magnitude smaller in size. In addition to being
substantially smaller in size, the reading and writing to text file in
XYZ format is very time consuming; by using binary array files, the
model run time was reduced from over four hours per simulated year to
less than ten minutes. The third file format utilized by ICM-Morph is
comma separated values (CSV) files. These files are used for storing a
large amount of input and output data that does not need to be mapped,
such as lookup tables and summary data that may be based on the
ICM-Hydro compartment polygons, or the ICM-LAVegMod grid cells.

While the binary arrays that are used by ICM-Morph drastically improved
model performance, with respect to memory usage and speed, there is one
distinct downside; these files are not human- nor universally
machine-readable. They can only be read by another Fortan program that
is being run on a computer processor that has the same “endianness”
settings as the computer used to originally generate the binary array.
Since all of the master plan simulations are conducted on Bridges-2,
this did not pose any issues with respect to the 2023 modeling efforts.

However, since there is a desire to use these data generated by
ICM-Morph in perpetuity, new post-processing programs were developed
that are able to be run, on Bridges-2, after the completion of the
entire simulation. These post-processors are run on Bridges-2 and read
in the binary arrays and convert them first back into ASCI XYZ text
files, then converts them into TIF raster files that are readable on any
computer that is equipped with GIS software. The TIFs are significantly
smaller than the ASCI XYZ text-based rasters, however there is
considerable computational time that is required to complete these file
conversions. The model run times are not impacted, however, because all
raster file conversions are now down outside of the ICM and can be done
in parallel, not impacting the overall ICM simulation run times for
production runs.

All of the raster post-processing Fortran source code files are
available on the CPRA Master Plan GitHub site[^1].

## 3.0 Initial Conditions and Input Files

A number of input files are required to be prepared prior to running the
ICM and ICM-Morph. As discussed in the previous section, many of these
input files are raster datasets that must be pre-processed into ASCI XYZ
format for initial use in ICM-Morph. While this pre-processing can be
cumbersome, and care must be taken to ensure that all raster files are
identical in resolution and extent, using these GIS-independent file
formats result in simple Fortran compiler settings with no external
dependencies such as NetCDF or GDAL libraries during model runtime
(GIS-readable rasters were generated as a post-processing step upon
completion of each simulation, as described in the previous section).
The following datasets are required to run ICM-Morph:

### 3.1 Topobathymetric Digital Elevation Model (DEM)

One of the primary datasets needed for ICM-Morph is a digital elevation
model (DEM) for the water and wetland areas of coastal Louisiana. The
preparation of a combined topographic and bathymetric (topobathymetric)
DEM is a complicated and time consuming process, which is fully
documented in [Attachment B1: Landscape Input
Data](https://coastal.la.gov/wp-content/uploads/2023/07/B1_LandscapeInputData_Jul2023_v5.pdf).
Once the initial conditions DEM was developed for the 2023 Coastal
Master Plan, minor edits were made to incorporate projects on the
landscape that were assumed to be in the FWOA simulations, but were not
present on the landscape when the elevation was collected. This updated
initial conditions DEM was then converted into the ASCI XYZ format and
read into ICM-Morph as the initial conditions elevation at the start of
the model run. At the very end of the model year, ICM-Morph saves the
DEM of the updated landscape to the model run folder, so that it can be
used for the starting elevation during the next model year.

### 3.2 Landscape Composition/Landtype

Similar to the DEM file, another important dataset for ICM-Morph is the
landscape composition. This refers to the land/water status of each
ICM-Morph pixel. In fact, ICM-Morph does not just differentiate between
land and water, but also further defines the land into: vegetated
wetland, unvegetated wetland, developed land/upland/fastlands, and
flotant marsh (Table 1). Like the DEM, the methodology for developing
the initial conditions version of this file is provided in [Attachment
B1: Landscape Input
Data](http://coastal.la.gov/wp-content/uploads/2023/07/B1_LandscapeInputData_Jul2023_v5.pdf).
This is also a file that is updated at the end of each model year by
ICM-Morph so that the changing landscape will be used as starting
conditions for each subsequent model year.

<caption>

Table 1. Values used in ICM-Morph for landtype classifications.
</caption>

| **Value** | **Landtype**                   |
|-----------|--------------------------------|
| 1         | Vegetated wetland              |
| 2         | Open water                     |
| 3         | Unvegetated wetland/bareground |
| 4         | Developed/upland/fastland      |
| 5         | Flotant marsh                  |

### 3.3 Polders

Leveed areas across the coastal zone were identified and delineated.
These poldered areas were mapped and are provided as an input file. All
ICM-Morph pixels that are located within a poldered area are treated as
developed/upland (landtype=4) and will not have any land change
processes applied to them (with the exception of elevation change due to
deep subsidence, as described below).

### 3.4 Deep Subsidence

Subsidence in ICM-Morph is treated as two separate rates. The first, is
a spatially heterogeneous rate that represents deep subsidence. This
rate map was derived from GPS benchmarks throughout coastal Louisiana
and remains unchanged in time; each year the same rates of deep
subsidence are applied, although the rates do vary spatially. The
methodology for deriving this input file is provided in [Attachment B3:
Determining Subsidence Rates for use in Predictive
Modeling](https://coastal.la.gov/wp-content/uploads/2023/08/B3_DeterminingSubsidenceRates_Mar2021_v3.pdf).

### 3.5 Shallow Subsidence

The second subsidence rate used in ICM-Morph is the shallow subsidence.
This rate is meant to represent processes in the coastal wetland areas
near surface – and the rates were derived from CRMS observations. These
rates are provided as an input CSV file and include spatial averages for
of shallow subsidence rates for each ecoregion (lower, upper, and median
quartiles per ecoregion). Again, the methodology for deriving this input
file is provided in [Attachment B3: Determining Subsidence Rates for use
in Predictive
Modeling](https://coastal.la.gov/wp-content/uploads/2023/08/B3_DeterminingSubsidenceRates_Mar2021_v3.pdf).

### 3.6 Marsh Edge Erosion Rates

Edge erosion is not dynamically modeled within ICM-Morph, but instead is
represented by projecting historic rates of edge erosion into the
future. Historic rates of marsh edge erosion were calculated from high
resolution aerial imagery; the methodology is provided in [Attachment
B1: Landscape Input
Data](https://coastal.la.gov/wp-content/uploads/2023/07/B1_LandscapeInputData_Jul2023_v5.pdf).
Once developed for initial conditions, this map of historic edge erosion
was updated to represent newly built shoreline protection features. Edge
erosion rates were set to zero for all areas assigned to be within the
influence area of a shoreline protection project. The influence area was
defined as within a 200 m buffer of shoreline protection feature’s
centerline.

### 3.7 Organic Matter Accumulation Rates

As discussed in later sections, the organic matter accumulation rates
used by ICM-Morph are provided as quartile ranges for each FFIBS
category. The rates are further divided into Chenier Plain and Delta
Plain regions. This information is passed into ICM-Morph via a CSV table
that provides OMAR values for each FFIBS category for each ecoregion.
Ecoregions were assigned to be either in the Delta Plain, or the Chenier
Plain, and assigned appropriate OMAR values from these assignments. This
allows for future changes where the OMAR values may be provided at a
more granular (e.g., ecoregion) scale. These regional OMAR values
(provided in Table 2, below) were updated based on soil surveys
conducted in 2018 and discussed in detail in [Attachment D2:
ICM-Wetlands, Vegetation & Soils Model
Improvements](https://coastal.la.gov/wp-content/uploads/2023/08/D2_2023ICM-Wetlands-Veg-Soils-Model-Improvements_Jun2020_v2.pdf)
(Baustian et al., 2020).

### 3.8 Active Deltaic Zones

In addition to the categorically defined OMAR values for the Chenier
Plain and Delta Plain, a third region was also differentiated in order
to better parameterize marsh accretionary processes within active,
river-connected, delta splays. These active deltaic OMAR rates are
applied only to fresh and/or intermediate marshes (e.g.,
river-connected) located within predefined footprints of active delta
splays. These predefined footprints are based on the ICM-Hydro
compartments. An input CSV file was prepared that flags each ICM-Hydro
compartment as being an active delta, or not. This input file was
initially set such that only ICM-Hydro compartments within the Bird’s
Foot Delta, the Atchafalaya/Wax Lake Deltas, and existing river
diversions (e.g., Davis Pond and Caernarvon) were treated as active
deltaic zones. However, for any simulation in which a newly proposed
diversion was implemented in the model (e.g., the Mid-Barataria and
Mid-Breton Sediment Diversions in FWOA simulations) this list of active
deltaic compartments was updated so that newly built land within the
diversion outfall areas would be treated, with respect to OMAR values,
as an active river-connected delta splay.

### 3.9 Maintained Navigation Channels and Rivers

Several water bodies in the Louisiana coastal zone receive regular
dredging in order to maintain proper channel depths for navigational
purposes. Therefore, all federal navigational channels in the coastal
zone have a separate DEM file that contains the initial bathymetric
elevation of the channel. This separate DEM is used, as described below,
to ensure that the elevation within these channels are maintained
throughout the entire simulation.

### 3.10 Model Grid Lookup Maps

As described above, several input files and datasets used by ICM-Morph
are provided via CSV tables that are based on ecoregions, ICM-Hydro
compartments, and ICM-LAVegMod grid cells. In order to properly link to
each of these lookup tables, three additional model grid map rasters are
necessary that map each ICM-Morph pixel to the corresponding
ICM-LAVegMod grid cell, each ICM-Hydro compartment, and each ecoregion.
These lookup maps all have identical resolution, extent, and NoData
settings as the DEM and landtype raster files described above.

### 3.11 Input Parameters Control File

The majority of information that is needed to run an instance of
ICM-Morph is passed into the compiled code via a CSV control file titled
*input_params.csv*. This file is generated programmatically by

**ICM.py** and is updated for every year of the simulation. Some of the
parameters included are constant for every model year, whereas other
parameters may change from one year to the next. An example of the
latter would be the files that are used to implement a master plan
project on the landscape during a specific year.

This file also contains the directory paths for the many input/output
files generated for, and by, ICM-Morph. The variables that are passed
into ICM-Morph via this *input_params* file are listed in the **set_io**
subroutine - a description of each variable is also provided in
**ICM.py**[^2].

## 4.0 ICM-Morph Model Subroutines

The following section steps through the model subroutines that are
contained in ICM-Morph. At the start of each section, two lists will be
provided; the left column contains the subroutine name of the ICM-Morph
code; the right column contains the name of the corresponding Fortran
source code file available on the CPRA Github site[^3]. Throughout this
section, as various subroutines (either internal to ICM-Morph, or other
ICM components) are referenced, the subroutine name will be
**emphasized** so that the cross-subroutine connections within the model
are evident.

### 4.1 Main Control Program and Parameter Settings

|                                  |                                  |
|----------------------------------|----------------------------------|
| Subroutine: **main**             | Source code: WM_main.f90         |
| Subroutine: **set_io**           | Source code: WM_set_io.f90       |
| Subroutine: **params_alloc_io**  | Source code: WM_params_alloc.f90 |
| Subroutine: **params_alloc**     | Source code: WM_params_alloc.f90 |
| Subroutine: **dem_params_alloc** | Source code: WM_params_alloc.f90 |
| Module: **params**               | Source code: WM_params.f90       |

The primary control program for ICM-Morph is housed within the **main**
subroutine which steps through all subsequent subroutines of ICM-Morph.
The first four subroutines that are called are **set_io**,
**params_alloc_io**, **params_alloc**, **dem_params_alloc**, and
**params**. These four subroutines are used to read input variables and
allocate variables and memory arrays that are required to run ICM-Morph.
Processing messages and model performance times are written to console
and log files throughout this subroutine.

### 4.2 Pre-processing Input Files

|                               |                                   |
|-------------------------------|-----------------------------------|
| Subroutine: **preprocessing** | Source code: WM_preprocessing.f90 |

Once input variables and file paths are initialized, the
**preprocessing** subroutine programmatically steps through all input
files and reads them into allocated memory arrays. This subroutine is
passed a binary flag that indicates whether the binary files
representing raster arrays (as discussed above) are used, or if the ASCI
XYZ text-based raster files are utilized. If the ASCI XYZ rasters are
used, this subroutine takes over 14 minutes to run on Bridges-2 for the
2023 Coastal Master Plan, whereas if binary files are used the
subroutine finishes in 42 seconds.

### 4.3 Inundation Depths

|                                   |                                       |
|-----------------------------------|---------------------------------------|
| Subroutine: **inundation_depths** | Source code: WM_inundation_depths.f90 |

The **inundation_depths** subroutine is used to calculate, for every
pixel within the ICM-Morph domain, a water depth. This subroutine uses
the initial topobathy DEM raster from the start of the model year to set
a bottom elevation for any given pixel. The water surface elevation,
however, can be calculated for one of fourteen different time periods;
which is set by the variable *tp*, which is passed into the subroutine
when it is called from **main**. This variable must have a value set to
any integer from 1 to 14. If *tp* is a value less than or equal to a
value of 12, then the *tp* represents the elapsed month of the year. If
*tp* is set to 13, it represents the annual average; if the value is 14,
it represents the annual average for the previous year. For example if
**inundation_depths** is called and is passed (*tp*=3), then the
subroutine will calculate the average inundation depth for the month of
March during the current simulation year. The twelve monthly (and annual
mean) water surface elevations are read into the model for each
ICM-Hydro compartment in a CSV file that is processed during
**preprocessing**. Both the current year and previous years’ water
surface elevations are read into the model. The previous year’s average
inundation depth is required since chronic inundation stress thresholds
used to define inundation loss require inundation depth data for two
consecutive years. This is discussed later during the section describing
the **inundation_thresholds** subroutine.

In addition to providing arrays of inundation depth for each pixel, this
subroutine also tabulates the number of pixels within each ICM-Hydro
compartment that were wet during each respective month. This wetted area
tabulation per compartment is utilized later during the
**mineral_deposition** subroutine.

### 4.4 Marsh Edge Delineation

|                                  |                                      |
|----------------------------------|--------------------------------------|
| Subroutine: **edge_delineation** | Source code: WM_edge_delineation.f90 |

*Marsh edge*, as defined in ICM-Morph refers to a 30-m pixel that is
classified as either vegetated wetland or bareground and is adjacent to
at least one pixel classified as open water. The entire landtype raster
is looped over using a 3 x 3 moving window; pixels that are classified
as either upland or flotant marsh are excluded from analysis since
neither of these two landtypes has any edge processes applied to them
and are therefore not included in the definition of marsh edge. The
identified marsh edge pixels are saved to a new raster with an assigned
integer value of 1; all other pixels within the domain that are not
marsh edge receive a value of 0.

### 4.5 Mineral Sediment Deposition

|                                    |                                        |
|------------------------------------|----------------------------------------|
| Subroutine: **mineral_deposition** | Source code: WM_mineral_deposition.f90 |

Mineral sediment deposition on the marsh surface and on water body
bottoms is calculated at a monthly timestep[^4]. This
**mineral_deposition** subroutine utilizes two datasets; the monthly
inundation depth data calculated during **inundation_depths** and mass
per unit area sediment depositional loads calculated monthly by
ICM-Hydro and read into ICM-Morph during **preprocessing**. The monthly
sediment depositional mass loading data is tabulated in ICM-Hydro to be
total mass of mineral sediment deposition over three different zones of
each ICM-Hydro compartment: the open water area, the marsh edge area
(same definition as used in **edge_delineation**), and the marsh
interior (which is set equal to the total marsh area within the
ICM-Hydro compartment less the area of marsh edge). ICM-Hydro uses
average elevations for each of these three zones to calculate sediment
deposition, it does not account for the topographical details at each
pixel in the same way that ICM-Morph does. Therefore, the first step in
this **mineral_deposition** subroutine is to correct the mineral
sediment loading rates so that the area over which the monthly
deposition occurs is only those areas within each ICM-Hydro compartment
that was actually inundated at some point during the month for which
sediment deposition is being calculated for. This is done by conserving
the mass of sediment deposited during the month for a given ICM-Hydro
compartment, but updating the area to include only those areas wetted
during the month. For example, if a given compartment were calculated to
deposit 100 g/m<sup>2</sup> of mineral sediment on the marsh interior
area, but only 50% of the marsh interior were at elevations low enough
to have been inundated during the respective month – then the areal
loading of sediment would be pro-rated by this 50% factor. This would
result in any marsh interior pixel inundated during the month to receive
200 g/m<sup>2</sup> and all dry marsh interior pixels would receive no
mineral sediment deposition.

Flotant marsh areas are treated as having no inundation depth, since
they will float on top of the water column; therefore, no mineral
sediment deposition is modeled to occur on flotant marsh areas within
ICM-Morph.

Once the mineral deposition mass per inundated area is determined, the
mass loading is converted to a vertical accretion due to mineral
deposition, in units of centimeters. This conversion to accretion depth
is done by dividing the mineral mass loading per area (units of
g/cm<sup>2</sup>) by the mineral sediment bulk density (units of
g/cm<sup>3</sup>). Open water areas use a bed sediment bulk density
value, and marsh areas use the mineral soil self-packing density as
determined via an ideal mixing model (described in [Attachment D2:
ICM-Wetlands, Vegetation & Soils Model
Improvements](https://coastal.la.gov/wp-content/uploads/2023/08/D2_2023ICM-Wetlands-Veg-Soils-Model-Improvements_Jun2020_v2.pdf)).
The values for both of these density variables are passed into ICM-Morph
via the input parameters control file. For 2023 modeling, the open water
bed bulk density and mineral self-packing density values were 0.835
g/cm<sup>3</sup> and 2.106 g/cm<sup>3</sup>, respectively.

The final step of the **mineral_deposition** subroutine is to apply
low-pass filter on the total annual mineral accretion to ensure that any
numerical instabilities that may have resulted in runaway accretionary
rates are kept in check. The maximum allowable mineral accretionary
rates are set in the input parameters control file and are set to 50 cm
and 10 cm of mineral sediment accretion in open water bodies and marsh
surfaces, respectively.

#### Erosion of Water Body Bottoms

In addition to mineral sediment accretion, the **mineral_deposition**
subroutine is also where erosion of water bottoms is calculated. For
open water areas, ICM-Hydro will report out negative mass depositional
rates if a water body is erosive, as opposed to depositional. Therefore,
for open water areas, the calculated accretion depths may be negative,
indicating net erosion for a given month. While the erosional
mass-to-erosional depth is calculated in the same manner as deposition
in water bodies (the same bulk density value of 0.835 g/cm<sup>3</sup>),
there is a separate erosion depth threshold applied for the total
cumulative erosion of water bodies. Like the accretionary limits, the
maximum erosion threshold value is passed into ICM-Morph via the input
parameters control file. For 2023, a value of 50 cm was used, resulting
in a maximum scour of water bodies of 50 cm in any given model year.

Unlike water body bottoms, erosion of sediments from the marsh interior
surface are not modeled in ICM-Morph; only marsh edge erosional
processes are incorporated. Marsh edge erosion is discussed in a later
section of this report.

### 4.6 Organic Matter Accumulation

|                                   |                                       |
|-----------------------------------|---------------------------------------|
| Subroutine: **organic_accretion** | Source code: WM_organic_accretion.f90 |

Once the mineral component of marsh accretion is calculated in the
previous subroutine, the **organic_accretion** subroutine is used to
determine the organic portion of marsh accretion. Organic accretion is
calculated as a function of vegetation species coverage (represented by
the weighted FFIBS score – as described above) and look-up tables of
organic matter accumulation rates (OMAR) that vary based on FFIBS
category and location. The OMAR values, like mineral sediment
deposition, are in units of mass per unit area (g/cm<sup>2</sup>); these
rates are converted to a vertical accretion depth by dividing OMAR by
the bulk density of the organic portions of the marsh soil. This organic
bulk density is represented by the self-packing density of organic
sediments, which was determined with the same ideal mixing model as used
for mineral sediment accretion (again, see [Attachment D2: ICM-Wetlands,
Vegetation & Soils Model
Improvements](https://coastal.la.gov/wp-content/uploads/2023/08/D2_2023ICM-Wetlands-Veg-Soils-Model-Improvements_Jun2020_v2.pdf)).

The lookup table (Table 2) for OMAR by FFIBS category was derived from
CRMS soil data which was partitioned into Deltaic Plain and Chenier
Plain zones. In addition, separate OMAR values were used for fresh and
intermediate marshes that are located within an active delta. The
*active delta* locations were added to the OMAR tables to represent
sites that are located in active delta splays with current river
connectivity. The regions of the model that are treated as *active
delta* are set via the active deltaic zone input file described above.

The OMAR lookup tables are derived from CRMS observations, and the rates
were averaged categorically for each FFIBS classification category.
However, ICM-LAVegMod calculates a weighted FFIBS score to better
represent the mixture of species present within each ICM-LAVegMod grid
cell. This weighted FFIBS score is used to linearly interpolate the
categorical OMAR data so that the organic accretionary processes will
represent the mixture of vegetation species present in each grid cell.
For example, if one ICM-LAVegMod grid cell has an even 50/50 mixture of
brackish and saline marsh vegetation than it will have a weighted FFIBS
score of 11.5; subsequently an interpolated median OMAR of 0.735
g/cm<sup>2</sup> would be assigned to this grid cell.

<caption>

Table 2. Organic matter accumulation rates (OMAR) and weighted FFIBS
ranges for categorical FFIBS data used in ICM-Morph
</caption>

| Organic Matter Accumulation Rates (g/cm<sup>2</sup>) by FFIBS Category and Location |  |  |  |  |
|----|----|----|----|----|
|  |  | Deltaic Plain | Chenier Plain | Active Delta |
| Forested Wetlands 0 ≤ FFIBS \< 0.15 | Lower | 0.079 | 0.079 | \- |
|  | Median | 0.093 | 0.093 | \- |
|  | Upper | 0.11 | 0.11 | \- |
| Fresh Marsh 0.15 ≤ FFIBS \< 1.5 | Lower | 0.073 | 0.04 | 0.145 |
|  | Median | 0.089 | 0.058 | 0.145 |
|  | Upper | 0.107 | 0.085 | 0.145 |
| Intermediate Marsh 1.5 ≤ FFIBS \< 5.0 | Lower | 0.068 | 0.033 | 0.145 |
|  | Median | 0.076 | 0.04 | 0.145 |
|  | Upper | 0.085 | 0.048 | 0.145 |
| Brackish Marsh 5.0 ≤ FFIBS \< 18.0 | Lower | 0.058 | 0.035 | \- |
|  | Median | 0.065 | 0.048 | \- |
|  | Upper | 0.074 | 0.05 | \- |
| Saline Marsh 18.0 ≤ FFIBS | Lower | 0.07 | 0.038 | \- |
|  | Median | 0.082 | 0.038 | \- |
|  | Upper | 0.097 | 0.038 | \- |

During model calibration tests, the quartile ranges for the OMAR values
were treated as a calibration parameter. Once the ICM-Hydro sediment
parameters were calibrated and adjusted to best match observed suspended
solids data, the modeled total accretion (e.g., mineral plus organic)
was compared to observed total accretion data in the CRMS network. Three
calculations of total accretion were made using the lower, median, and
upper estimates of OMAR in Table 2; all using the ICM-Hydro determined
rates of mineral sediment deposition. As seen in Figure 1, the median
OMAR values provided the best performance with respect to total
accretion. Therefore, for all 2023 Coastal Master Plan simulations, the
median OMAR values for each FFIBS category were used.

<figure>
<img
src="https://github.com/stacycalhoun/MP29ReportTemplates/blob/main/ICM-WetlandMorph_Figs/ICM-WetlandMorph-Fig1.png"
alt="A bar plot of CRMS total accretion compared to ICM total accretion coastwide average. Active deltaic accretion was greatest, followed by average delta plain accretion, and lastly average chenier plain accretion. Average active deltaic accretion included both fresh and intermediate marsh values." />
<figcaption aria-hidden="true">Figure 1. Calibration of total accretion using three different quartiles
for organic matter accumulation rates.</figcaption>
</figure>


Future improvements to ICM-Morph could allow for the OMAR quartile
selection to be made by model ecoregion, allowing for greater spatial
heterogeneity with how ICM-Morph models organic accretionary processes.
For 2023, all ecoregions are assigned the median OMAR rates based on
Figure 1.

### 4.7 Flotant Marsh Loss

|                         |                             |
|-------------------------|-----------------------------|
| Subroutine: **flotant** | Source code: WM_flotant.f90 |

The **flotant** subroutine is the first ICM-Morph subroutine that will
start to track changes to the landscape as a result of hydrologic
conditions for the given simulation year. This is done by updating a
raster array of *land change flags*. ICM-LAVegMod reports out a
percentage of flotant marsh within each grid cell that died during the
current simulation year (see Flotant Mortality and Flotant Establishment
sections above). ICM-Morph receives the percentage of flotant marsh
within each grid cell and calculates the total number of flotant pixels
within each respective grid cell that need to be converted to open
water. This subroutine scans the landtype pixel on a grid cell by grid
cell basis – flagging the appropriate number of flotant marsh pixels
that are to be converted to open water. For any pixel flagged for
conversion, the land change flag raster is updated to have a value of -2
for each respective dead flotant pixel (Table 3). This flag is later
used by the **update_landtype** subroutine.

<caption>

Table 3. Land change flag values.
</caption>

| Land change flag value | Description |
|----|----|
| 0 | No change |
| -1 | Conversion from vegetated wetland to open water due to inundation |
| -2 | Conversion from flotant marsh mat to open water |
| -3 | Conversion from marsh edge to open water due to erosion |
| 1 | Conversion from open water to new subaerial land eligible for vegetation |

### 4.8 Implement Shoreline Protection Projects

|  |  |
|----|----|
| Subroutine: **build_shoreline_projects** | Source code: WM_build_shoreline_projects.f90 |

This subroutine, **build_shoreline_projects**, is the first step in
ICM-Morph where the landscape processes are updated to account for a
master plan project; specifically the two projects types (shoreline
protection and bankline stabilization) that impact marsh edge erosion.
This is accomplished by reading in a new raster data file (in ASCI XYZ
format) that contains a map of multiplication factors, with a default
value of 1.0. The area of edge that is protected by the proposed project
will have mapped multiplication factors of 0.0; indicating that all edge
erosion processes are zeroed out in this footprint. This subroutine
reads in this map of erosion multiplication factors and multiplies them
by the historic marsh edge erosion rates read in during
**preprocessing**. This is done prior to the next subroutine,
**edge_erosion**, so that the implemented shoreline protection projects
are incorporated before marsh edge erosion is calculated for the
simulation year. The overall control program **ICM.py** tracks which
projects are active, and for all years after a shoreline protection
project has been implemented, these multiplication factors will be read
into ICM-Morph; this results in the impacts on marsh edge erosion rates
are continuously applied in all model years after implementation.

### 4.9 Marsh Edge Erosion

|                              |                                        |
|------------------------------|----------------------------------------|
| Subroutine: **edge_erosion** | Source code: WM_build_edge_erosion.f90 |

After the locations of marsh edge are identified for any given model
year by **edge_delineation**, these edge pixels are overlaid with a
static raster of historic marsh edge erosion rates that were quantified
from high resolution aerial imagery (see [Attachment B1: Landscape Input
Data](https://coastal.la.gov/wp-content/uploads/2023/07/B1_LandscapeInputData_Jul2023_v5.pdf)).
The historic rates are used, in conjunction with the number of elapsed
model years to determine whether any given edge pixel would be eroded
during a given model year. This subroutine reads the assigned historic
marsh edge erosion rate (m/yr) for a given pixel and calculates the
number of years needed to erode one 30-m pixel. For example, a pixel
with an assigned historic edge erosion rate of 15 m/yr would only be
subjected to edge erosion processes once every 2 years; a pixel with an
edge erosion rate of 5 m/yr would be subjected to edge erosion processes
once every 6 years.

If an edge pixel is identified as being subject to erosion loss during a
model year, the subroutine updates the land change flag value for the
respective pixel to an integer value of -3, which indicates conversion
from marsh edge to open water due to erosion (Table 3).

This subroutine is only capable of removing one whole edge pixel in a
given year, neither partial pixels, nor erosion events more frequent
than once per year are able to be incorporated. Therefore, durations to
erode one pixel are always rounded up to the next year; this has the
possibility to lead to an under-representation of historic edge retreat.
For example, if a location has an erosion rate of 29 m/yr, this would
require 1.03 years to erode one 30-m pixel. However, this code rounds
that up to erode one edge pixel every 2 years. Over a 50-year
simulation, the historic erosion rate would have resulted in 1,450
meters of erosion, whereas this subroutine would only remove 25 pixels,
resulting in 750 meters of modeled erosion. If the historic rate was 16
m/yr, this would result in 800 m over 50 years, or again 750 meters of
modeled erosion. This could be changed by simply rounding instead of
using a ceiling function on the year calculation.

### 4.10 Map Spatial Distribution of Bareground

|                                |                                    |
|--------------------------------|------------------------------------|
| Subroutine: **map_bareground** | Source code: WM_map_bareground.f90 |

One change to the overall ICM-Morph and ICM-LAVegMod processes for 2023
was the spatial handling of bareground within each ICM-LAVegMod grid
cell. ICM-LAVegMod tracks the percentage of a grid cell that is
unvegetated land, but there is no spatial knowledge in ICM-LAVegMod of
exactly which ICM-Morph pixels should be unvegetated. A decision was
made during the model improvements process to assume that all bareground
within a grid cell should be allocated to the lowest elevated pixels
classified as land within each grid cell. This subroutine loops through
the land pixels within each grid cell, and compares each respective
pixel elevation to all other pixels within that grid cell.

ICM-LAVegMod classifies two types of bareground, new bareground that was
vegetated in the previous year and lost vegetative cover due to
hydrologic conditions, and old bareground that was previously
unvegetated and remained so due to a lack of appropriate vegetation
available for establishment (see Update Water Area and Update Vegetation
Coverages sections of ICM-LAVegMod above).

Therefore this bareground assignment to the lowest elevation pixels is
performed twice; first for the portion of each grid cell that is *old
bareground* (this will be the lowest elevation land), and second for the
portion that is *new bareground* (this will be higher in elevation than
*old bareground* but lower in elevation than all vegetated land pixels).
If the pixel is identified as *old bareground*, a bareground flag will
be set to value of 1 for the given pixel; otherwise the *new bareground*
will have a flag value of 2. This bareground flag is used later in the
**update_elevation** subroutine.

### 4.11 Chronic Inundation Stress Thresholds

|  |  |
|----|----|
| Subroutine: **inundation_thresholds** | Source code: WM_inundation_thresholds.f90 |

Per the recommendations provided in Baustian et al. (2020), inundation
depth thresholds are calculated for each landscape pixel to determine
whether chronic inundation exists at a given location, which would
result in the inability for vegetated wetland to persist. As stated in
Appendix A of Baustian et al. (2020):

> The upper limit for vegetation occurrence can be defined as occurring
> at the upper bound of the 95th or 99th percentile confidence interval
> of the observed inundation depths, which varies as a function of the
> observed salinity at the same location over the same time frame. More
> specifically, this is the 97.5th and 99.5th quantile, respectively of
> the annual mean depth-salinity data pairs in CRMS from 2010 through
> 2017 and was calculated from a quantile regression analysis (Bolker,
> 2008) using the R statistical software (R Core Team, 2019). The 99.5th
> quantile of the data set was recommended as the water depth limitation
> which varies by salinity for vegetation occurrence. The depth limit,
> DL, at a given salinity, S, is:

$$DL_S = \mu_{depth,S} + Z * \sigma_{depth, S}$$

<div align="right">

<caption>

Eq. 1
</caption>

</div>

$$\mu_{depth,S} = \beta_0 +\beta_1 * S$$

<div align="right">

<caption>

Eq. 2
</caption>

</div>

$$\sigma_{depth, S} = \beta_2 + \beta_3e^{\beta_4*S}$$

<div align="right">

<caption>

Eq. 3
</caption>

</div>

> where a quantile regression was used to determine the mean and
> standard deviation of inundation depths as a function of salinity,
> $\mu_{depth,S}$ and $\sigma_{depth, S}$ , respectively. The upper
> bounds on the 95th and 99th percentile confidence intervals are
> calculated with Z<sub>95</sub> = 1.96 and Z<sub>99.5</sub> = 2.576,
> respectively.

So as to have flexibility with defining the salinity-inundation depth
threshold curve, the parameter values for $\beta_0$, $\beta_1$,
$\beta_2$, $\beta_3$, $\beta_4$, and $Z_{pct}$ are all provided as input
variables (written by **ICM.py** into the input parameters control text
file) read into the program during **preprocessing**. The values
selected by Baustian et al. (2020), representing the 99.5th percentile
regression curve (Figure 2) are provided in Table 4.

<caption>

Table 4. Values used to define the salinity-inundation depth threshold
relationships for the 2023 Coastal Master Plan
</caption>

| Variable  | Value for 2023 Coastal Master Plan |
|-----------|------------------------------------|
| $\beta_0$ | 0.0058                             |
| $\beta_1$ | -0.00207                           |
| $\beta_2$ | 0.0809                             |
| $\beta_3$ | 0.0892                             |
| $\beta_4$ | -0.19                              |
| $Z_{pct}$ | 2.57                               |

<figure>
<img
src="https://github.com/stacycalhoun/MP29ReportTemplates/blob/main/ICM-WetlandMorph_Figs/ICM-WetlandMorph-Fig2.png"
alt="A scatterplot of the distribution annual mean salinity (ppt) and annual mean inundation depth (m) at CRMS stations. The 99.5th quartile, represented by a dashed blue line, indicates the salinity inundation threshold curve above which vegetation will not be able to persist in the ICM." />
<figcaption aria-hidden="true">Figure 2. Salinity and water depth distribution of vegetation occurrence
in coastal Louisiana with the blue line representing the
salinity-inundation depth threshold curve used in ICM-Morph. Each point
represents annual mean salinity and annual mean depth for a single CRMS
station (2010-2017), n = 1944. The blue line is the 99.5th quantile, and
the grey line is the 0.05th quantile from a quantile regression
analysis.</figcaption>
</figure>


If a given pixel is classified as non-flotant, non-forested wetland and
it is inundated during both the present and previous model year to a
mean annual depth that is greater than the inundation-salinity
relationship defined by the above equations, then the pixel is flagged
for inundation loss. The subroutine updates the land change flag value
for the respective pixel to an integer value of -1 (Table 3), which
indicates conversion from wetland to open water due to chronic
inundation stress.

This chronic inundation stress criterion is not applied to wetland
pixels that are either flotant marsh or fresh forested wetlands. The
weighted FFIBS score from a pixel’s overlying ICM-LAVegMod grid cell is
examined, and if that score is less than or equal to 0.15, the
vegetative cover is assigned as fresh forested wetlands and the
inundation criterion is not applied to the respective pixel. If a
wetland pixel was already flagged for loss in edge_erosion, then the
land change flag value is not overwritten by this subroutine and the
pixel will be considered lost due to erosion during the year, not
inundation.

### 4.12 Update Elevation

|                                  |                                      |
|----------------------------------|--------------------------------------|
| Subroutine: **update_elevation** | Source code: WM_update_elevation.f90 |

Information derived in previous subroutines (**mineral_deposition**,
**organic_accretion**, **flotant**, **edge_erosion**, and
**map_bareground**) are used in this subroutine to update the elevation
of each pixel within the model domain. First, if a pixel is flagged as
flotant loss, the elevation of the pixel is lowered to an elevation 1.0
m below mean water surface elevation. Second, if a pixel is flagged as
edge erosion, the elevation of the pixel is lowered to an elevation 0.25
m below mean water surface elevation. Third if a pixel is flagged as
being old bareground, the elevation of the pixel is lowered 0.05 m.
These three values are each set in the input parameters control file
passed into ICM-Morph.

After the three special treatments described above, each pixel’s
elevation is updated based on the following logic:

- If a pixel is vegetated land, the change in elevation is the sum of:
  - calculated mineral accretion
  - calculated organic accretion
  - deep subsidence, and
  - shallow subsidence
- If a pixel is open water, the change in elevation is the sum of:
  - calculated mineral accretion and
  - deep subsidence.
- If a pixel is unvegetated land/bareground, the change in elevation is
  the sum of:
  - calculated mineral accretion
  - deep subsidence, and
  - shallow subsidence.
- If a pixel is upland/developed land, or located within a pre-defined
  polder the change in elevation is equal to deep subsidence only.
- If a pixel is flotant marsh, its elevation is unchanged (unless it was
  previously flagged as dead flotant, as described above).

The next step in this subroutine is to reset all elevation values for
water pixels that are located within maintained (e.g., dredged)
navigational channels. These locations are set to always maintain the
initial elevation; this is handled via the *maintained navigation
channels and rivers* input file.

Finally, if a pixel was identified as being located within the barrier
island model domain, the elevations as calculated by ICM-BI, and
processed by **ICM.py**, are used. These barrier island elevations are
passed into ICM-Morph and handled during **preprocessing**.

### 4.13 Classify New Subaerial Land

|                                    |                                        |
|------------------------------------|----------------------------------------|
| Subroutine: **new_subaerial_land** | Source code: WM_new_subaerial_land.f90 |

Once the landscape elevation has been updated, this subroutine compares
the annual mean inundation depths (from **inundation_depths**) for each
open water pixel for both the current and previous model years. If both
inundation depths are less than 0.1 m, than the pixel is assumed to be
shallow enough to allow for vegetation establishment and the land change
flag value (Table 4) for the pixel will be set to 1, indicating new
subaerial land eligible for vegetation establishment during the next
model year. The depth criterion for vegetation establishment is set in
the input parameters control file and passed into ICM-Morph during
**set_io**.

### 4.14 Update Landtype Classification

|                                  |                                     |
|----------------------------------|-------------------------------------|
| Subroutine: **updated_landtype** | Source code: WM_update_landtype.f90 |

This subroutine updates the landtype raster to account for the updated
*land change flags* (Table 4) that had been updated in **edge_erosion**,
**flotant**, **inundation_thresholds**, and **new_subaerial_land
subroutines**.

### 4.15 Calculate Inundation for Use in HSI Equations

|  |  |
|----|----|
| Subroutine: **inundation_HSI_bins** | Source code: WM_inundation_HSI_bins.f90 |

This subroutine is used solely to summarize seasonal water depths that
are required by some of the water fowl habitat suitability indices
within ICM-HSI (see [Attachment C10: 2023 Habitat Suitability Index
(HSI)
Model](https://coastal.la.gov/wp-content/uploads/2023/08/C10_2023HSIModel_Feb2023_v3.pdf)).
The various HSI equations use different seasons that are important to
the individual waterfowl species, so the monthly inundation depths
calculated in **inundation_depths** are used here to develop the
appropriate mean depths important to each species. The calculated mean
depths are then tabulated, for each grid cell, for a variety of depth
bins that are used directly in the HSI equations.

### 4.16 Implement Marsh Creation and Landbridge Projects

|  |  |
|----|----|
| Subroutine: **build_marsh_projects** | Source code: WM_build_marsh_project.f90 |

The last step in ICM-Morph for a given model year is to build landscape
projects so that the newly built projects are on the landscape for the
start of the next model year. This is first done by building marsh
creation and landbridge projects (see [Appendix F: Project Concepts for
descriptions](https://coastal.la.gov/wp-content/uploads/2023/05/F_ProjectConcepts_Apr2023_v2.pdf)).
These project types must have a pre-built raster file (in ASCI XYZ
format) and during the year of implementation **ICM.py** will pass the
filepath to these rasters into ICM-Morph. The format for these project
files should be a XYZ map that has a value of NoData for all non-project
areas, and a value of for the design elevation within the project
footprint. The design elevation is defined here as the *height above
mean water surface elevation*. When run, this subroutine will take this
value, determine the background water surface elevation at the project
site, and determine what the vertical elevation of the built marsh
should be relative to the vertical datum used by ICM-Morph.

Once the design elevation has been determined, this subroutine will
examine all pixels within and compare the water depth to the marsh
creation depth threshold value assigned to each project. For most marsh
creation projects, this depth threshold is set to 0.76 m; any pixel that
is water and has a depth greater than this threshold will be ineligible
for marsh creation. Landbridge projects due not have this same threshold
applied – they are modeled to build marsh throughout the entire project
footprint, regardless of water depth.

For all areas that are within the project footprint, and not flagged by
the depth threshold, the pixel elevation will be set to be equal to the
design elevation. The landtype of that pixel will also be set to 3,
indicating the land is newly built, unvegetated, and available for
vegetation to establish in ICM-LAVegMod during the next model year.

### 4.17 Implement Ridge Restoration and Levee Projects

|  |  |
|----|----|
| Subroutine: **build_ridge_projects** | Source code: WM_build_ridge_project.f90 |

This subroutine is very similar to **build_marsh_projects**, but is
instead used to update the landscape for ridge restoration and
structural risk reduction project types (e.g., levees). There are two
minor differences though. First, there is no water depth threshold
applied to ridge and levee project footprints; all pixels under the
project footprint will be built. Second, the input files for these
project types use a simple design elevation (relative to the ICM-Morph
vertical datum); unlike the marsh creation rasters which had provided an
elevation relative to the mean water surface.

### 4.18 Post-Processing/Summary Outputs

|  |  |
|----|----|
| Subroutine: **summaries** | Source code: WM_summaries.f90 |
| Subroutine: **write_output_summaries** | Source code: WM_write_output_summaries.f90 |
| Subroutine: **write_output_QAQC_points** | Source code: WM_write_output_QAQC_points.f90 |
| Subroutine: **write_output_binary_rasters** | Source code: WM_write_output_binary_rasters.f90 |
| Subroutine: **write_output_asci_rasters** | Source code: WM_write_output_asci_rasters.f90 |

After the implementation of marsh creation, ridge restoration, and levee
features, the landscape is done being modified for the model year. All
of the subroutines within ICM-Morph are simply post-processing data and
writing output files for use in other ICM subroutines. This includes
several CSV files that contain spatially averaged data for each
ICM-Hydro compartment that is used to define the landscape in ICM-Hydro
the following model year. Additional spatial averages are prepared for
the ICM-LAVegMod grid cells that are used in ICM-LAVegMod as well as
ICM-HSI. These various files, as well as additional summary output files
are written to files in **write_output_summaries** subroutine.

For 2023, ICM-Morph was updated to include a new output file that was
used for quality assurance/quality control (QAQC) during model
simulations. These files, referred to as QAQC save points, contain a
snapshot of all pertinent variables that are used throughout ICM-Morph
to define how a given pixel changes in land elevation, and/or landtype,
for any given model year. All data is saved for the year and appended as
a new row to the end of the file. Therefore, at the end of the
multidecade simulation, for each QAQC save point, there will be one file
that tracks all elevation and land change data at that point for all
model years. The locations used for the QAQC save points were randomly
placed throughout the entire model domain. In addition to the random
locations, QAQC save points were also placed at select transects that
had been pre-identified at areas of interest (Wax Lake Outlet, near
proposed diversions, etc.). Finally, QAQC save points were placed at
every CRMS location within the model domain. In all, there are 2,941
locations where 16 different variables from across ICM-Hydro,
ICM-LAVegMod, and ICM-Morph are all reported out. These files are
written during the **write_output_QAQC_points** subroutine.

In addition to the summary CSV files saved to disk, ICM-Morph also
exports raster based datasets to disk so that they can be used, in
binary format, as the starting conditions for the next model year in
ICM-Morph. These binary rasters can also be converted, using various
post-processsing scripts (described above), to file formats that are
GIS-readable. When output rasters are saved as binary files, the
**write_output_binary_rasters** subroutine completes in just 4 seconds
on Bridges-2 for the 2023 Master Plan. If chosen in the input parameters
control file, ICM-Morph can also be used to export the rasters in ASCI
XYZ format. This is considerably slower, and would also impact run times
in the subsequent model year. Therefore, all model simulations currently
export to binary, and utilize external post-processing script if and
when GIS-readable versions of the data are needed.

## 5.0 Modeling Submerged Aquatic Vegetation (SAV)

In addition to the modeling processes described in the sections above,
ICM-Morph is also where the newly developed model for submerged aquatic
vegetation (SAV) is housed. To incorporate the SAV model, two additional
subroutines were added to ICM-Morph. These are the **distance_to_land**
and **sav** subroutines. The first of these two subroutines is the most
time consuming subroutine of the entire ICM-Morph model; it takes
approximately 1 hour to complete on Bridges-2. This subroutine
calculates, for every pixel classified as open water, how far that pixel
is from land. This calculation must be omnidirectional; however, the SAV
model has an upper distance of 500 m from land, so the search can be
abbreviated once the subroutine determines a pixel is at least 500 m
from land.

The second subroutine is the Fortran representation of the statistical
model for SAV presence, as discussed in depth in [Attachment D3:
ICM-Wetlands – Submerged Aquatic Vegetation (SAV)
Updates](https://coastal.la.gov/wp-content/uploads/2023/07/D3_ICM-WetlandsSAVUpdate_Jul2023_v2.pdf).

In addition to the spatial data that is already used by and available in
ICM-Morph, there are several statistical parameters that must be read
into ICM-Morph. These are passed into the model via the input parameters
control file.

## 6.0 References

Baustian, M. M., Reed, D., Visser, J., Duke-Sylvester, S., Snedden, G.,
Wang, H., DeMarco, K., Foster-Martinez, M., Sharp, L. A., McGinnis, T.,
& Jarrell, E. (2020). 2023 Coastal Master Plan: Attachment D2:
ICM-Wetlands, Vegetation, and Soil Model Improvements. Version 2.
(p. 90). Baton Rouge, Louisiana: Coastal Protection and Restoration
Authority.

Couvillion, B. R. (2023). 2023 Coastal Master Plan: Attachment B1:
Landscape Input Data. Version 5. (p. 43). Baton Rouge, Louisiana:
Coastal Protection and Restoration Authority.

DeMarco, K., Schoolmaster, D., & Couvillion, B. R. (2023). 2023 Coastal
Master Plan: Attachment D3: ICM-Wetlands – Submerged Aquatic Vegetation
(SAV) Updates. Version 02. (pp. 1-58). Baton Rouge, Louisiana: Coastal
Protection and Restoration Authority.

Fitzpatrick, C., Jankowski, K. L., & Reed, D. (2021). 2023 Coastal
Master Plan: Attachment B3: Determining Subsidence Rates for Use in
Predictive Modeling. Version 3. (p. 71). Baton Rouge, Louisiana: Coastal
Protection and Restoration Authority.

Lindquist, D. C. (2023). 2023 Coastal Master Plan: Attachment C10: 2023
Habitat Suitability Index (HSI) Model. Version 3. (pp. 1-35). Baton
Rouge, Louisiana: Coastal Protection and Restoration Authority.

White, E., Meselhe, E., McCorquodale, A., Couvillion, B., Dong, Z.,
Duke-Sylvester, S. M., & Wang, Y. (2017). 2017 Coastal Master Plan:
Attachment C3-22: Integrated Compartment Model (ICM) Development
(pp. 1–49) \[Version Final\]. Coastal Protection and Restoration
Authority.

[^1]: <https://github.com/CPRA-MP/ICM_MorphRasters>

[^2]: [input params.csv with descriptions of variables, via
    Github](https://github.com/CPRA-MP/ICM/blob/da25a31d1107a9445f5cbabe4b8a5b4785af3ad3/ICM.py#L2586)

[^3]: <https://github.com/CPRA-MP/ICM_Morph>

[^4]: Previous versions of the ICM mapped mineral sediment deposition
    annually. This resulted in a uniform sediment deposition thickness
    over all landscape pixels that had been inundated at least once
    during the simulated year. Improvements to the models for the 2023
    plan converted this to the monthly depositional mapping that is
    described here.
