"""
align_and_count_features.py
Dylan H. Ross

description:
    Performs the processing steps necessary to take raw RNAseq reads, align them to a reference genome, and count
    the aligned genomic features in all samples. This script assumes a bash environment and requires external programs
    to be installed and accessible on the system:
        * hisat2
        * samtools
        * featureCounts
"""


import os
from json import load as jload, dump as jdump
import sys
from subprocess import run, DEVNULL
from re import compile as recomp
from datetime import datetime


def _help_msg():
    """
    _help_msg

    prints information about the arguments expected by this script and what they do
    """
    msg = "DHRutil.RNAseq.align_and_count_features\n\n" \
          "Performs the processing steps necessary to take raw RNAseq reads, align them to a reference genome, " \
          "and count the aligned genomic features in all samples. This script assumes a bash environment and " \
          "requires external programs to be installed and accessible on the system:\n* hisat2\n* samtools\n* " \
          "featureCounts" \
          "\nA configuration file is used to define the input files and set all of the parameters for the processing " \
          "steps. A template can be generated for the configuration file by calling this script with the " \
          "--make-config option" \
          "\n\nUsage:\n\tpython3 -m DHRutil.RNAseq.align_and_count_features [--config config.json /" \
          " --make-config / --help]\n\nOptions:\n\t--config config.json (specify the configuration file to use)" \
          "\n\t--make-config (generate a template for the configuration file)" \
          "\n\t--help (print this message and exit)\n"
    print(msg)


def _make_sample_config_file():
    """
    _make_sample_config_file

    writes an example configuration file to 'config.json'
    """
    config = {
        "input_raw_sequence_reads_1": [
            "ctl_A_1.fq.gz", "ctl_B_1.fq.gz", "trt_C_1.fq.gz", "trt_D_1.fq.gz"
        ],
        "input_raw_sequence_reads_2": [
            "ctl_A_2.fq.gz", "ctl_B_2.fq.gz", "trt_C_2.fq.gz", "trt_D_2.fq.gz", 
        ],
        "cpu_threads": 16,
        "hisat2_index_filename_prefix": "index/genome",
        "rm_sam_after_sorting": True,
        "featureCounts_gtf_annotation_file": "annotation.gtf",
        "featureCounts_output": "ctl_trt_counts_raw.txt"
    }
    with open('config.json', 'w') as f:
        jdump(config, f, indent=4)


def _load_config(config_file):
    """
    _load_config

    loads the configuration file

    Paramters
    ---------
    config_file : str
        path to the configuration file

    Returns
    -------
    config : dict
        dictionary with configuration parameters
    """
    # ensure the config file exists before trying to load it
    if not os.path.isfile(config_file):
        e = "_load_config: configuration file {} does not exist".format(config_file)
        raise ValueError(e)

    # load the configuration
    with open(config_file, 'r') as j:
        return jload(j)


def _hisat2_align_seq_reads(input_fq_1, input_fq_2, index_prefix, threads):
    """
    _hisat2_align_seq_reads

    aligns paired sequence reads (gzipped) to reference genome using hisat2
    stores output with alignment results in <base_name>.align.log

    Paramters
    ---------
    input_fq_1 : str
        first paired sequence read file, gzipped, must end in _1.fq.gz
    input_fq_2 : str
        second paired sequence read file, gzipped, must end in _2.fq.gz
    index_prefix : str
        prefix for the genome index filename (i.e. everything but .X.ht2)
    threads : int
        CPU threads to use when running hisat2

    Returns
    -------
    sam : str
        output .sam file
    """
    align_rate_pat = recomp(r'([0-9]+[.][0-9]+)% overall alignment rate')
    print("aligning paired reads {} and {} ... ".format(input_fq_1, input_fq_2), end="", flush=True)
    if input_fq_1[-8:] != '_1.fq.gz' or input_fq_2[-8:] != '_2.fq.gz':
        e = '_hisat_2_align_seq_reads: -1 and -2 paired read inputs must be named as ' \
            '<desc>_1.fq.gz and <desc>_2.fq.gz (input 1: "{}", input 2: "{}"")'.format(input_fq_1, input_fq_2)
        raise ValueError(e)
    base_name = input_fq_1[:-8]
    output_sam = base_name + '.sam'
    log_file = base_name + '.align.log'
    cmd = ["hisat2", "-q", "-p", str(threads), "--pen-noncansplice", "1000000", 
           "-x", index_prefix, "-1", input_fq_1, "-2", input_fq_2, "-S", output_sam]
    res = run(cmd, text=True, capture_output=True)
    if res.returncode:
        e = "_hisat_2_align_seq_reads: hisat2 returned with nonzero exit code: {} ({}, {})"
        raise RuntimeError(e.format(res.returncode, input_fq_1, input_fq_2))
    align_rate = float(align_rate_pat.search(res.stderr).group(1))
    if align_rate < 70:
        e = "_hisat_2_align_seq_reads: alignment rate of {:.2f}% is below 70% threshold ({}, {})"
        raise RuntimeError(e.format(align_rate, input_fq_1, input_fq_2))
    print("alignment rate: {:.2f}%".format(align_rate), flush=True)
    # write the results (from stdout) to a log file
    with open(log_file, 'w') as f:
        f.write(' '.join(res.args) + '\n')
        f.write(res.stderr)
    return output_sam


def _samtools_sort(sam, threads, mem="1G", rm_sam_after_sort=False):
    """
    _samtools_sort

    converts and sorts the .sam files
    if rm_sam_after_sorting is set to True, removes sam file after sorting

    Paramters
    ---------
    sam : str
        input .sam file 
    threads : int
        CPU threads to use for sorting
    mem : str

    rm_sam_after_sort : bool
        indicates whether the input .sam file should be removed after sorting to save space
    
    Returns
    -------
    sort_bams : list(str)
        list of output sorted .bam files
    """
    print("sorting {} ... ".format(sam), end="", flush=True)
    base_name = os.path.splitext(sam)[0]
    output_bam = base_name + '.sort.bam'
    cmd = ["samtools", "sort", "-m", mem, "-T", "tmp", "--threads", str(threads), sam, "-o", output_bam]
    res = run(cmd, stdout=DEVNULL, stderr=DEVNULL)
    if res.returncode:
        e = "_samtools_sort: samtools sort returned with nonzero exit code: {} ({})"
        raise RuntimeError(e.format(res.returncode, sam))
    print("ok", flush=True)
    # (optional) remove .sam files to save space
    if rm_sam_after_sort:
        print("removing {} ... ".format(sam), end="", flush=True)
        os.remove(sam)
        print("ok", flush=True)
    return output_bam


def _featurecounts_all_samples(sort_bams, output, gtf_annotation, threads):
    """
    _featurecounts_all_samples

    uses featureCounts to count the features from all samples
    creates <output>.txt (counted features) and <output>.log (log of featureCounts output)

    Parameters
    ----------
    sort_bams : list(str)
        list of sorted .bam files
    output : str
        name of output file with counted features
    gtf_annotation : str
        GTF annotation file name
    threads : int
        CPU threads to use for feature counting
    """
    log_file = os.path.splitext(output)[0] + '.log'
    print("running featureCounts on all samples ...", end="", flush=True)
    with open(log_file, 'wb') as f:
        cmd = ["featureCounts", "-p", "-t", "exon", "-a", gtf_annotation, "-g", "gene_name", "-T", 
               str(threads), "-o", output]
        cmd += sort_bams
        cmd_str = " ".join(cmd) + "\n"
        f.write(cmd_str.encode())
        res = run(cmd, stdout=f, stderr=f)
        if res.returncode:
            e = "_featurecounts_all_samples: featureCounts returned with nonzero exit code: {}"
            raise RuntimeError(e.format(res.returncode))
    print(" done", flush=True)


def _main():
    """
    _main

    main execution sequence
    """
    # keep track of elapsed processing time
    start_time = datetime.now()

    # parse the arguments to this script, if there are any issues print the help message, then raise an error
    if len(sys.argv) < 2:
        # no option specified
        _help_msg()
        e = "_main: no options provided"
        raise ValueError(e)
    elif sys.argv[1] == "--help":
        # print the help message and exit
        _help_msg()
        exit()
    elif sys.argv[1] == "--make-config":
        # make a new configuration file and exit
        _make_sample_config_file()
        exit()
    elif sys.argv[1] == "--config":
        # make sure a configuration file is provided
        if len(sys.argv) < 3:
            # no configuration file is provided
            _help_msg()
            e = "_main: no configuration file provided"
            raise ValueError(e)
        else:
            # try to load the configuration file
            cfg = _load_config(sys.argv[2])
    else:
        # unrecognized option
        _help_msg()
        e = "_main: unrecognized option: {}".format(sys.argv[1])
        raise ValueError(e)

    # using the loaded configuration file, perform the analysis

    # 1. use hisat2 to align sequence reads to reference genome
    sams = []
    for fq_1, fq_2 in zip(cfg['input_raw_sequence_reads_1'], cfg['input_raw_sequence_reads_2']):
        sams.append(_hisat2_align_seq_reads(fq_1, fq_2, cfg['hisat2_index_filename_prefix'], cfg['cpu_threads']))

    # 2. sort and convert .sam files with samtools
    sort_bams = []
    for sam in sams:
        sort_bams.append(_samtools_sort(sam, cfg['cpu_threads'], rm_sam_after_sort=cfg['rm_sam_after_sorting']))

    # 3. count features from all samples
    _featurecounts_all_samples(sort_bams, cfg['featureCounts_output'], cfg['featureCounts_gtf_annotation_file'], 
                               cfg['cpu_threads'])

    # report the elapsed processing time
    print('total processing time: {:.1f} minutes'.format((datetime.now() - start_time).seconds / 60.))


if __name__ == '__main__':
    _main()

