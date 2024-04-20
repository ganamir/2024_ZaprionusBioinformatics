# 2024_ZaprionusBioinformatics
## Important sites:
https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/

https://hcc.unl.edu/docs/applications/app_specific/bioinformatics_tools/alignment_tools/bowtie2/

https://www.ebi.ac.uk/sites/ebi.ac.uk/files/content.ebi.ac.uk/materials/2014/140217_AgriOmics/dan_bolser_snp_calling.pdf

## Methods:
### Find .fastq.gz files and copy them to my directory

``` find /weka/data/lab/rudman/Schmidt2018/Zap/ -type f -name "*gz*" -exec cp {} /scratch/user/amir.gabidulin/20240418_231031/Practice2/FASTQ_Files \; ``` 

### Find .fastq.gz files, subset samples to n number of lines, and save to SubsetData folder, maintaining the same names and formats. 

``` find /weka/scratch/user/amir.gabidulin/20240418_231031/Practice2/FASTQ_Files/ -type f -name "*gz" -exec sh -c 'zcat "$1" | head -n 2000 | gzip >  "/weka/scratch/user/amir.gabidulin/20240418_231031/Practice2/SubsetData/$(basename "$1" .gz)".subset.fastq.gz' _ {} \; ```

### cutadaptTrimmer.sh >>> Search for my .fastq.gz files, both treatment and founders, and insert them appropriately into the cutadapt function with 20 filtering quality and my adapter sequences:

```
#!/bin/bash

module load cutadapt

# Define adapter sequences
ADAPTER_FWD="AGATCGGAAGAGC"
ADAPTER_REV="AGATCGGAAGAGC"

# Define quality cutoff
QUAL=20

# Define input and output directories
INPUT_DIR="/scratch/user/amir.gabidulin/20240418_231031/Practice2/SubsetData"
OUTPUT_DIR="/scratch/user/amir.gabidulin/20240418_231031/Practice2/TrimmedData"

# Process regular paired-end files
for R1_file in "$INPUT_DIR"/*_R1_001.fastq.subset.fastq.gz; do
    # Extract corresponding R2 filename
    R2_file="${R1_file/_R1_001/_R2_001}"

    # Extract basename of the files
    base_name=$(basename "$R1_file" "_R1_001.fastq.subset.fastq.gz")

    # Define output filenames with full path
    OUTPUT_FWD="$OUTPUT_DIR/${base_name}.trimmed_R1.fastq.gz"
    OUTPUT_REV="$OUTPUT_DIR/${base_name}.trimmed_R2.fastq.gz"

    # Run Cutadapt
    cutadapt -a $ADAPTER_FWD -A $ADAPTER_REV \
             --quality-cutoff $QUAL \
             -o $OUTPUT_FWD -p $OUTPUT_REV \
             $R1_file $R2_file

    # Optionally, add additional processing or logging here
done

# Process additional paired-end files with different naming convention
for T_file in "$INPUT_DIR"/T*_FOUND_R1.fastq.subset.fastq.gz; do
    # Extract corresponding R2 filename
    R2_file="${T_file/_R1/_R2}"

    # Extract basename of the files
    base_name=$(basename "$T_file" "_R1.fastq.subset.fastq.gz")

    # Define output filenames with full path
    OUTPUT_FWD="$OUTPUT_DIR/${base_name}.trimmed_R1.fastq.gz"
    OUTPUT_REV="$OUTPUT_DIR/${base_name}.trimmed_R2.fastq.gz"

    # Run Cutadapt
    cutadapt -a $ADAPTER_FWD -A $ADAPTER_REV \
             --quality-cutoff $QUAL \
             -o $OUTPUT_FWD -p $OUTPUT_REV \
             $T_file $R2_file

    # Optionally, add additional processing or logging here
done
 ```

### Build a Reference Genome from your Genome Assembly (dm6 for DGRP)

``` bowtie2-build dm6.fa.gz reference_genome ```


### bowtie2AlignScript.sh >>> Align reads to your reference index

```
#!/bin/bash

module load bowtie2

# Define paths
BOWTIE2_INDEX="/scratch/user/amir.gabidulin/20240418_231031/Practice2/DGRPassembly/reference_genome"
OUTPUT_DIR="/scratch/user/amir.gabidulin/20240418_231031/Practice2/AlignedReads"
INPUT_DIR="/scratch/user/amir.gabidulin/20240418_231031/Practice2/TrimmedData"
LOG_FILE="$OUTPUT_DIR/bowtie2_log.txt"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Function to run Bowtie2 alignment
run_bowtie2() {
    local SAMPLE_NAME="$1"
    local TRIMMED_FWD="$2"
    local TRIMMED_REV="$3"
    local OUTPUT_SAM="$OUTPUT_DIR/${SAMPLE_NAME}.sam"

    echo "Aligning $TRIMMED_FWD and $TRIMMED_REV..."

    # Run Bowtie2 alignment
    bowtie2 -x "$BOWTIE2_INDEX" \
            -1 "$TRIMMED_FWD" \
            -2 "$TRIMMED_REV" \
            -S "$OUTPUT_SAM" \
            2>> "$LOG_FILE"

    # Optionally, add additional processing or logging here
}

# Align regular paired-end files
for R1_file in "$INPUT_DIR"/*.trimmed_R1.fastq.gz; do
    if [ -f "$R1_file" ]; then
        SAMPLE_NAME=$(basename "$R1_file" ".trimmed_R1.fastq.gz")
        TRIMMED_REV="$INPUT_DIR/${SAMPLE_NAME}.trimmed_R2.fastq.gz"
        run_bowtie2 "$SAMPLE_NAME" "$R1_file" "$TRIMMED_REV"
    else
        echo "Warning: $R1_file does not exist or is not a file."
    fi
done

```





