# 2024_ZaprionusBioinformatics

### Find .fastq.gz files and copy them to my directory

``` find /weka/data/lab/rudman/Schmidt2018/Zap/ -type f -name "*gz*" -exec cp {} /scratch/user/amir.gabidulin/20240418_231031/Practice2/FASTQ_Files \; ``` 

### Find .fastq.gz files, subset samples to n number of lines, and save to SubsetData folder, maintaining the same names and formats. 

``` find /weka/scratch/user/amir.gabidulin/20240418_231031/Practice2/FASTQ_Files/ -type f -name "*gz" -exec sh -c 'zcat "$1" | head -n 2000 | gzip >  "/weka/scratch/user/amir.gabidulin/20240418_231031/Practice2/SubsetData/$(basename "$1" .gz)".subset.fastq.gz' _ {} \; ```
