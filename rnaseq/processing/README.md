
### `align_and_count_features.py`

Performs the processing steps necessary to take raw RNAseq reads, align them to a reference genome, and count
the aligned genomic features in all samples. This script assumes a bash environment and requires external programs
to be installed and accessible on the system:
* `hisat2`
* `samtools`
* `featureCounts`

_Usage_
 

This utility operates as a standalone script and can be called directly. First generate a template configuration file:
```bash
python3 align_and_count_features.py --make-config
```
This will produce `config.json` with contents similar to:
```json
{
    "input_raw_sequence_reads_1": [
        "ctl_A_1.fq.gz", "ctl_B_1.fq.gz", "trt_C_1.fq.gz", "trt_D_1.fq.gz"
    ],
    "input_raw_sequence_reads_2": [
        "ctl_A_2.fq.gz", "ctl_B_2.fq.gz", "trt_C_2.fq.gz", "trt_D_2.fq.gz", 
    ],
    "cpu_threads": 16,
    "hisat2_index_filename_prefix": "index/genome",
    "rm_sam_after_sorting": true,
    "featureCounts_gtf_annotation_file": "annotation.gtf",
    "featureCounts_output": "ctl_trt_counts_raw.txt"
}
```
After editing the configuration file to suit the desired run conditions and input files, the complete analysis can be 
run using the following command:
```bash
python3 align_and_count_features --config config.json
```
