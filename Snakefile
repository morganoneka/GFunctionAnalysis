configfile: "config.yaml"

CODE_LOCATION = config["CODE_LOCATION"]

GCROSS_LOCATION = CODE_LOCATION + "GcrossMultiFile.R"
GCROSSTHEO_LOCATION = CODE_LOCATION + "GcrossMultiFileTheo.R"
PLOTXY_LOCATION = CODE_LOCATION + "PlotXYMultiFile.R"
PLOTGC_LOCATION = CODE_LOCATION + "PlotGcross.R"
CALAUC_LOCATION = CODE_LOCATION + "CalculateAUC.R"
CALAUCTHEO_LOCATION = CODE_LOCATION + "CalculateAUCTheo.R"

DATA_DIR = config["DATA_DIR"]

DATA_DIR_GCROSS = DATA_DIR + "Gcross/"
DATA_DIR_PLOTXY = DATA_DIR + "PlotXY/"
DATA_DIR_PLOTGC = DATA_DIR + "PlotGC/"

INPUT_DIR = config["INPUT_DIR"]

XPOS = config["XPOS"]
YPOS = config["YPOS"]
RADIUS_MAX = config["RADIUS_MAX"]
SEP = config["SEP"]
JSON_LOCATION = config["JSON_LOCATION"]
RADII = config["RADII"]

rule all:
    input:
        DATA_DIR + "GcrossSummary.csv",
        DATA_DIR + "GcrossSummaryTheo.csv",
        DATA_DIR + "XYsummary.csv",
        DATA_DIR + "GcrossPlotSummary.csv",
        DATA_DIR + "GcrossAUC.csv",
        DATA_DIR + "GcrossAUC_THEORETICAL.csv"

rule create_dirs:
    input:
        directory(INPUT_DIR)
    output:
        directory(DATA_DIR_GCROSS),
        directory(DATA_DIR_PLOTXY),
        directory(DATA_DIR_PLOTGC)
    shell:
        "mkdir -p \"{DATA_DIR_GCROSS}\" ; mkdir -p \"{DATA_DIR_PLOTXY}\" ; mkdir -p \"{DATA_DIR_PLOTGC}\""

rule run_gcross:
    input:
        directory(INPUT_DIR),
        JSON_LOCATION,
        directory(DATA_DIR_GCROSS)
    output:
        DATA_DIR + "GcrossSummary.csv"
    shell:
        "Rscript \"{GCROSS_LOCATION}\" -i \"{input[0]}\" -o \"{DATA_DIR_GCROSS}\" --radiusmax {RADIUS_MAX} "
        "--xcol {XPOS} --ycol {YPOS} "
        "--sep \'{SEP}\' --jsonphenotype \"{input[1]}\" --summary \"{output}\" -v FALSE"


rule run_gcross_theo:
    input:
        directory(INPUT_DIR),
        JSON_LOCATION,
        directory(DATA_DIR_GCROSS)
    output:
        DATA_DIR + "GcrossSummaryTheo.csv"
    shell:
        "Rscript \"{GCROSSTHEO_LOCATION}\" -i \"{input[0]}\" -o \"{DATA_DIR_GCROSS}\" --radiusmax {RADIUS_MAX} "
        "--xcol {XPOS} --ycol {YPOS} "
        "--sep \'{SEP}\' --jsonphenotype \"{input[1]}\" --summary \"{output}\" -v FALSE"

rule plot_xyloc:
    input:
        directory(INPUT_DIR),
        DATA_DIR + "Interactions.json",
        directory(DATA_DIR_PLOTXY)
    output:
        DATA_DIR + "XYsummary.csv"
    shell:
        "Rscript \"{PLOTXY_LOCATION}\" -i \"{input[0]}\" -o \"{DATA_DIR_PLOTXY}\" "
        "--xcol {XPOS} --ycol {YPOS} "
        "--sep \'{SEP}\' --jsonphenotype \"{input[1]}\" --summary \"{output}\" -v FALSE"

rule plot_gcross:
    input:
        directory(DATA_DIR_GCROSS)
    output:
        DATA_DIR + "GcrossPlotSummary.csv"
    shell:
        "Rscript \"{PLOTGC_LOCATION}\" -i \"{input[0]}\" -o \"{DATA_DIR_PLOTGC}\" "
        "--summary \"{output}\""

rule calculate_auc:
    input:
        directory(DATA_DIR_GCROSS)
    output:
        DATA_DIR + "GcrossAUC.csv"
    shell:
        "Rscript \"{CALAUC_LOCATION}\" -i \"{input[0]}\" -o \"{DATA_DIR}\" --radius=\"{RADII}\""

rule calculate_auc_theo:
    input:
        directory(DATA_DIR_GCROSS)
    output:
        DATA_DIR + "GcrossAUC_THEORETICAL.csv"
    shell:
        "Rscript \"{CALAUCTHEO_LOCATION}\" -i \"{input[0]}\" -o \"{DATA_DIR}\" --radius=\"{RADII}\""
