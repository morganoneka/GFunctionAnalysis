# Spatial Analysis using G-Function Analysis

The G-function is a useful tool for analyzing infiltration of cell types from mIHC, histology, and other spatial biological data. This code provides a workflow for analyzing and visualizing these results.

## Required to run
- Snakemake
- The following R libraries:
  - ggplot2
  - spatstat
  - matlab

( possibly several others I'll double-check )



## How to Run
There are three files you need to pay attention to:
- **Snakefile**: The snakemake file that runs all the code.
- **config.yaml**: A configuration file that helps snakemake run on your specific data. You will need to edit this to identify where your data is and what the various column names are.
- **Interactions.json**: A json file indicating phenotypes of interest.

There are examples of both the config and interactions files included in this repository to illustrate what they look like.

### Input File
For this analysis, you will need files with all of these columns:
- X coordinate
- Y coordinate
- Phenotype(s). This pipeline supports both "wide" data (many phenotype columns with true/false values) or "narrow" data (one phenotype column with many different possible labels).

This file can be a .csv, .txt, or other file that has one delimiter. It is important that all of your files have the exact same column names. They also must be in one folder.

### Editing config.yaml
Next, you'll want to **edit the config.yaml file** so the Snakemake file knows where all the relevant files are. These are the values in the Snakemake file:
- **CODE_LOCATION**: The directory that contains all of the Gcross code. This will be the location where this repo was cloned to.
- **DATA_DIR**: The directory where all your output will be saved
- **INPUT_DIR**: The path to the files being analyzed.
- **XPOS** and **YPOS**: The columns identifying the X and Y coordinates of each cell
- **PATIENT_COL** and **SAMPLE_COL**: The columns with identifiers for each patient and sample.
- **RADIUS_MAX**: The maximum radius to compute Gcross for. For example, if this value is 500, the Gcross value will be computed for radius=1 to radius=500.
- **JSON_LOCATION**: The location of the interaction .json file.
- **RADII**: The radii to calculate the AUC for.

### Editing Interactions.json
This indicates what interactions you are interested in analyzing.

#### Reference vs. Non-Reference
Each interaction contains a "reference" and "non-reference" cell type. The Gcross analysis will calculate the infiltration of non-reference cells into reference cells. For example, suppose you were interested in quantifying the infiltration of CTLs into epithelial cells. Epithelial cells would be the reference cells, and CTLs would be the non-reference cells.

Each cell type is a list of dictionaries.

#### Multiple Cell Types
Reference and non-reference cells can have require specific values of multiple columns.

The example below shows the definition of a CTL. For this sample, we want a cell to be considered a CTL if its value in the CD8 column is "CD8+" **AND** if its value in the CD3 column is "CD3+".

This results in the non-reference value mapping to a list of one dictionary. Each key in the dictionary maps to the value it must have to be considered a CTL.

``` json
"Non-Reference" : [ {"CD8" : ["CD8+"], "CD3" : ["CD3+"]} ]
```

Cell types can also require only one of multiple conditions to be true.

The example below shows the definition of an epithelial cell. For this sample, a cell is considered to be epithelial if its value in the Ecadherin column is "ECadherin+" **OR** its value in the Pancytokeratin column is "Pancytokeratin+".

In this case, the reference value maps to a list of multiple dictionary. Within each dictionary, each key in the dictionary maps to the value it must have to fulfill one of the definitions of an epithelial cell.

``` json
"Reference" : [ {"ECadherin" : ["ECadherin+"]} , {"Pancytokeratin" : ["Pancytokeratin+"]}]
```

### Running the snakemake pipeline
Once the above files have been edited, open the command line, navigate to the location of the snakemake file, and enter the command `snakemake`.

### Additional notes
- The snakemake, config, and phenotype file do not have to be in the same folder as the Gcross code. I find it helpful to save it in the project directory, but in a separate folder from the data.
- This current workflow automatically generates a ton of plots, which can take up a lot of space. I'm working on a version that excludes plotting.
- There is a file `CalculateAUCKbins.R` that I haven't included in the workflow yet. It calculates the AUC for "bins" of radii (for example, 0 to 20 pixels, then 20 to 40, etc.). The code is there and generally clean though, so it should be somewhat self-explanatory to run.

## Current Workflow
This is specifically what the code does:

### XY Plotting
For each sample, and for each interaction, a plot will be generated showing the XY coordinates of the reference and non-reference cells in that sample.

### Gcross Analysis
For each interaction, a .csv file will be generated containing the Gcross values for each sample.

### Gcross Plotting
For each sample and interaction, a plot will be generated of that Gcross curve.

### AUC Calculation
One .csv file will be generated containing the AUC for each interaction at the desired pixel values.
