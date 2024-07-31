# 2024_ZaprionusBioinformatics: 
**Email at amir.gabidulin@wsu.edu if you have questions about the pipeline / need help.**


## Important sites:
https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/

https://hcc.unl.edu/docs/applications/app_specific/bioinformatics_tools/alignment_tools/bowtie2/

https://www.ebi.ac.uk/sites/ebi.ac.uk/files/content.ebi.ac.uk/materials/2014/140217_AgriOmics/dan_bolser_snp_calling.pdf

https://www.htslib.org/workflow/wgs-call.html

https://marineomics.github.io/POP_03_poolseq.html

# Methods:

## Scripts for copying files and taking their subsets to practice running the pipeline:

### Find .fastq.gz files and copy them to my directory

``` find /weka/data/lab/rudman/Schmidt2018/Zap/ -type f -name "*gz*" -exec cp {} /scratch/user/amir.gabidulin/20240418_231031/Practice2/FASTQ_Files \; ``` 

### Find .fastq.gz files, subset samples to n number of lines, and save to SubsetData folder, maintaining the same names and formats. 

``` find /weka/scratch/user/amir.gabidulin/20240418_231031/Practice2/FASTQ_Files/ -type f -name "*gz" -exec sh -c 'zcat "$1" | head -n 2000 | gzip >  "/weka/scratch/user/amir.gabidulin/20240418_231031/Practice2/SubsetData/$(basename "$1" .gz)".subset.fastq.gz' _ {} \; ```
## Working with raw files (FASTQ - FASTQ.gz)

### FASTQ raw files:
````fastqc 18057XD-04-06_S0_L001_R1_001.fastq.gz````


### FASTQ QC filtering: trim_galore.sh >>> 
**Explanation:** Specify Directory with DIR var where your FASTQs are at -> Insert directories of R1 and R2 corresponding files into R1 and R2 variables. 

**Purpose:** Trimming Adapter Seqs, Filtering anything with sub 40 BPs length, and allowing a max of 1 Ns in the sequences. Trimming and truncating at quality 20 BPs (phred+33 score). This is stringent filtering, lower coverage (< 50) or lower length of reads (< 100) will require changes. ZaprionusEXP had length of reads ~150, and a coverage of 100-250. 

```
#!/bin/bash
#SBATCH --partition=cas
#SBATCH --time=3-00:0:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --mail-type=ALL
#SBATCH --mail-user=amir.gabidulin@wsu.edu
#SBATCH --output=outputTRIM.out
#SBATCH --error=errorTRIM.err

module load trimgalore

# Specify the directory containing the fastq files
DIR="/weka/scratch/user/amir.gabidulin/20240503_022212/FinalDestination/20240510_132726/FASTQ_Files"

# Iterate over all R1 files in the directory
for R1 in "$DIR"/*_R1_001.fastq.gz; do
    # Derive the corresponding R2 file
    R2="${R1/_R1_/_R2_}"

    # Run trim_galore on the file pair
    trim_galore --paired --length 40 --max_n 1 "$R1" "$R2" -a AGATCGGAAGAGC -q 20 -o /weka/scratch/user/amir.gabidulin/20240503_022212/FinalDestination/20240510_132726/FinalVersion/TrimmedFASTQ
done
 ```

### Build a Reference Genome from your Genome Assembly (dm6 for DGRP)

``` bwa index dm6.fa.gz ```


### bwa_aligner.sh >>> Align reads to your reference index

**Explanation:** Look for .gz files in the input_dir, properly find the R1 and R2 pairs, and run a paired-end read alignment of reads to reference. 

**Purpose:** To call variants at different sites you need to understand where the sequence genome deviates from the reference. So BWA alignment aligns your .fastq reads against the reference genome and outputs that alignment information to .SAM (Sequence Alignment Map) file that are used for all the downstream NGS (Poolseq & WGS) analyses.

```
#!/bin/bash
#SBATCH --partition=cas
#SBATCH --time=7-00:0:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=50
#SBATCH --mail-type=ALL
#SBATCH --mail-user=amir.gabidulin@wsu.edu
#SBATCH --output=outputBWA.out
#SBATCH --error=errorBWA.err

module load bwa

# Directory containing .gz files
input_dir="/weka/scratch/user/amir.gabidulin/20240503_022212/FinalDestination/20240510_132726/FinalVersion/TrimmedFASTQ"
output_dir="$input_dir/bwaAligned"

mkdir -p "$output_dir"

# Reference genome
reference_genome="/weka/scratch/user/amir.gabidulin/20240503_022212/FinalDestination/20240510_132726/FinalVersion/DGRPassembly/dm6.fa"

# Number of threads for BWA
threads=50

# Iterate over all R1 files in the specified directory
for file1 in "$input_dir"/*_R1_*_val_1.fq.gz; do
    # Construct the expected R2 filename
    file2="${file1/_R1_/_R2_}"
    file2="${file2/_val_1/_val_2}"
    echo "$file1"
    echo "$file2"

    # Check if the R2 file exists
    if [ -f "$file2" ]; then
        header=$(zcat "$file1" | head -n 1)
        id=$(echo "$header" | head -n 1 | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/_/g')
        sm=$(echo "$header" | head -n 1 | grep -Eo "[ATGCN]+$")
        read_group=$(echo "@RG\tID:$id\tSM:${id}_${sm}\tLB:${id}_${sm}\tPL:ILLUMINA")

        echo "Read Group: $read_group"

        # Output filename
        output_file=$(basename "${file1%_R1_*.fq.gz}_aligned.sam")

        bwa mem "$reference_genome" -M -t "$threads" -v 3 -R "$read_group" "$file1" "$file2" > "$output_dir/$output_file"
    else
        echo "Corresponding R2 file for $file1 not found."
    fi
done
```

### samtobam.sh >>> converts .sam files to .bam (Binary *sequence* Alignment Map) files:

```
#!/bin/bash

module load samtools

# Define paths
SAM_DIR="/scratch/user/amir.gabidulin/20240418_231031/Practice2/AlignedReads"
BAM_DIR="/scratch/user/amir.gabidulin/20240418_231031/Practice2/bamFiles"

# Create BAM directory if it doesn't exist
mkdir -p "$BAM_DIR"

# Convert SAM files to BAM
for sam_file in "$SAM_DIR"/*.sam; do
    if [ -f "$sam_file" ]; then
        base_name=$(basename "$sam_file" .sam)
        bam_file="$BAM_DIR/$base_name.bam"
        samtools view -bS "$sam_file" | samtools sort -o "$bam_file"
        echo "Converted $sam_file to $bam_file"
    else
        echo "Warning: $sam_file does not exist or is not a file."
    fi
done
```
### Important to check file integrity, what aligned and didn't; Command below allows to check the 10 most common sequences that didn't align
```
samtools view -f 4 18057XD-04-06_S0_L001.bam | cut -f 10 | sort | uniq -c | sort -nr | head
```


### markAndRemove.sh >>> Mark duplicates and remove them from the bam files.

```
#!/bin/bash
#SBATCH --partition=cas
#SBATCH --time=7-00:0:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mail-type=ALL
#SBATCH --mail-user=amir.gabidulin@wsu.edu
#SBATCH --output=outputREMOVAL.out
#SBATCH --error=errorREMOVAL.err


module load samtools
module load picard

# Specify the directory containing the BAM files
BAM_DIR="/weka/scratch/user/amir.gabidulin/20240503_022212/FinalDestination/20240510_132726/FinalVersion/TrimmedFASTQ/bwaAligned/bamFiles"

# Loop through each BAM file in the specified directory
for BAM_FILE in "$BAM_DIR"/*.bam; do
  # Extract the base name of the BAM file (without path and extension)
  BASE_NAME=$(basename "$BAM_FILE" .bam)

  # Step 1: Mark duplicates using Picard
  picard MarkDuplicates \
      I="$BAM_FILE" \
      O="$BAM_DIR/${BASE_NAME}_marked.bam" \
      M="$BAM_DIR/${BASE_NAME}_metrics.txt" \
      REMOVE_DUPLICATES=true
  echo "Processed $BAM_FILE: duplicates marked and unpaired reads removed."
done

echo "All BAM files have been processed."
```
### coverage.sh >>> Generate coverage metrics -- Required for filtering and other calculated metrics.
```
#!/bin/bash
#SBATCH --partition=cas
#SBATCH --time=2-00:0:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=15
#SBATCH --mail-type=ALL
#SBATCH --mail-user=amir.gabidulin@wsu.edu
#SBATCH --output=outputCOVERAGE.out
#SBATCH --error=errorCOVERAGE.err

module load samtools
# Directory containing .bam files
BAM_DIR="/weka/scratch/user/amir.gabidulin/20240503_022212/FinalDestination/20240510_132726/FinalVersion/TrimmedFASTQ/bwaAligned/bamFiles/DuplicateRemoved"

# Output directory for coverage files
OUTPUT_DIR="/weka/scratch/user/amir.gabidulin/20240503_022212/FinalDestination/20240510_132726/FinalVersion/TrimmedFASTQ/bwaAligned/bamFiles/DuplicateRemoved/coverage"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop over all .bam files in the directory
for BAM_FILE in "$BAM_DIR"/*.bam
do
    # Extract the base name of the file (without the .bam extension)
    BASE_NAME=$(basename "$BAM_FILE" .bam)

    # Define the output coverage file name
    OUTPUT_FILE="$OUTPUT_DIR/$BASE_NAME.coverage"

    # Run samtools coverage and redirect the output to the coverage file
    samtools coverage "$BAM_FILE" -o "$OUTPUT_FILE"
done
```

### Combine .bam files into a single mpileup
```
#!/bin/bash
#SBATCH --partition=cas
#SBATCH --time=7-00:0:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mail-type=ALL
#SBATCH --mail-user=amir.gabidulin@wsu.edu
#SBATCH --output=outputCOMB.out
#SBATCH --error=errorCOMB.err

module load samtools
# Directory containing .bam files
BAM_DIR="/weka/scratch/user/amir.gabidulin/20240503_022212/FinalDestination/20240510_132726/FinalVersion/TrimmedFASTQ/bwaAligned/bamFiles/DuplicateRemoved"

samtools mpileup -f /weka/scratch/user/amir.gabidulin/20240503_022212/FinalDestination/20240510_132726/FinalVersion/DGRPassembly/dm6.fa -R "$BAM_DIR"/*.bam > combined.mpileup

```

### Combine everything into a .sync file using Popoolation2 tools:

```
#!/bin/bash
#SBATCH --partition=cas
#SBATCH --time=7-00:0:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mail-type=ALL
#SBATCH --mail-user=amir.gabidulin@wsu.edu
#SBATCH --output=outputSYNC.out
#SBATCH --error=errorSYNC.err
#SBATCH --mem-per-cpu=12G

module load java

java -ea -Xmx12g -jar /weka/scratch/user/amir.gabidulin/20240503_022212/FinalDestination/20240510_132726/FinalVersion/popoolation2_1201/mpileup2sync.jar --input combined.mpileup --output combined.sync --fastq-type sange$
```

### Use Popoolation2 tools to generate known indel sites into a .gtf format from the .mpileup (Run concurrently with the combination of .mpileup into .sync)

```
#!/bin/bash
#SBATCH --partition=cas
#SBATCH --time=7-00:0:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mail-type=ALL
#SBATCH --mail-user=amir.gabidulin@wsu.edu
#SBATCH --output=outputFINDIND.out
#SBATCH --error=errorFINDIND.err

module load perl
perl /weka/scratch/user/amir.gabidulin/20240503_022212/FinalDestination/20240510_132726/FinalVersion/popoolation2_1201/indel_filtering/identify-indel-regions.pl --input combined.mpileup --output knownIndels.gtf
```

### Use Popoolation2 tools to remove known indel sites from a .sync file from the generated .gtf file (This one removes suprisingly a lot of SNPs, up to a million in zap dataset)

```
#!/bin/bash
#SBATCH --partition=cas
#SBATCH --time=7-00:0:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mail-type=ALL
#SBATCH --mail-user=amir.gabidulin@wsu.edu
#SBATCH --output=outputFILTIND.out
#SBATCH --error=errorFILTIND.err

module load perl
perl /weka/scratch/user/amir.gabidulin/20240503_022212/FinalDestination/20240510_132726/FinalVersion/popoolation2_1201/indel_filtering/filter-sync-by-gtf.pl --input combined.sync --gtf knownIndels.gtf --output FinalClean.sync
```

### Generate Population Allele Frequencies (important _rc file)

```
#!/bin/bash
#SBATCH --partition=cas
#SBATCH --time=2-00:0:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mail-type=ALL
#SBATCH --mail-user=amir.gabidulin@wsu.edu
#SBATCH --output=outputAF.out
#SBATCH --error=errorAF.err

perl /weka/scratch/user/amir.gabidulin/20240503_022212/FinalDestination/20240510_132726/FinalVersion/popoolation2_1201/snp-frequency-diff.pl --input /weka/scratch/user/amir.gabidulin/20240503_022212/FinalDestination/20240510_132726/FinalVersion/TrimmedFASTQ/bwaAligned/bamFiles/DuplicateRemoved/DownStreamAnalysis/FinalClean.sync --output-prefix pop_diff --min-count 20 --min-coverage 50 --max-coverage 200
```

### Generate Population Fst

```
#!/bin/bash
#SBATCH --partition=cas
#SBATCH --time=2-00:0:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mail-type=ALL
#SBATCH --mail-user=amir.gabidulin@wsu.edu
#SBATCH --output=outputFst.out
#SBATCH --error=errorFst.err

perl /weka/scratch/user/amir.gabidulin/20240503_022212/FinalDestination/20240510_132726/FinalVersion/popoolation2_1201/fst-sliding.pl --input /weka/scratch/user/amir.gabidulin/20240503_022212/FinalDestination/20240510_132726/FinalVersion/TrimmedFASTQ/bwaAligned/bamFiles/DuplicateRemoved/DownStreamAnalysis/FinalClean.sync --output ZapVsCont.fst --min-count 20 --min-coverage 50 --max-coverage 200 --min-covered-fraction 0.2 --window-size 10000 --step-size 10000 --pool-size 200 --suppress-noninformative
```
