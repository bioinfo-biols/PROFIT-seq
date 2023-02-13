# PROFIT-seq: Targeted and accurate whole transcriptome sequencing in a programmable manner

Programmable full-length isoform transcriptome sequencing

## Author

Authors: Jinyang Zhang(zhangjinyang@biols.ac.cn), Fangqing Zhao(zhfq@biols.ac.cn)

Maintainer: Jinyang Zhang

## Release Notes

- version 1.0: First released version

## License

The code is released under the MIT License. See the LICENSE file for more detail

## Citing PROFIT-seq

Under submission

# Documentation

## 1. Installation

### 1.1 Install MinKNOW and corresponding gpu version of ont-guppy from the [Nanopore community](https://community.nanoporetech.com/downloads)

**NOTE: A valid ONT customer account is needed to download the MinKNOW and guppy software.** 

### 1.2 Replace embedded cpu version of guppy in MinKNOW to the gpu version  

Use gpu version of guppy in `/etc/systemd/system/guppyd.service`

``` 
service guppyd stop
# Modify guppyd service to use ont-guppy-gpu and dna_r9.4.1_450bps_fast.cfg 
systemctl daemon-reload
service guppyd start
```

Change the following settings in `/opt/ont/minknow/conf/package/sequencing/sequencing_MIN106_DNA.toml`

```
[analysis_configuration.read_detection]
break_reads_after_seconds = 0.4
```

### 1.3 Create password for the minknow user, and install python3 using system package manager of anaconda

**NOTE: you need to run PROFIT-seq as user `minknow` to access the guppy and MinKNOW server correctly.** 

### 1.4 Install PROFIT-seq

- Log in into user `minknow`

```bash
ssh minknow@localhost
```

- Clone this repository to your computer

```bash
cd /opt
git clone --recursive https://github.com/bioinfo-biols/PROFIT-seq.git
cd PROFIT-seq
```

- Create a virtual environment for running these scripts

- **NOTE: please use python 3.8.10**

```bash
virtualenv venv
source ./venv/bin/activate
```

- Install requirements & the corresponding version of `ont-pyguppy-client-lib`

```bash
pip install -r requirements.txt
pip install ont-pyguppy-client-lib==5.1.15
```

Note: if you only want to perform PROFIT-seq analysis (e.g. run PROFIT-seq using HPC), use step0_requirements.txt instead

```bash
make lib
```

## 2. Usage

### 2.1 Connect the MinION device and insert the configuration test / FLO-MIN106D flow cell correctly

Open the MinKNOW software to make sure flow cell have been successfully recognized.

### 2.2 Run PROFIT-seq using the `minknow` user

- Run PROFIT-seq from `minknow` user

```bash
ssh minknow@localhost
cd /opt/PROFIT-seq && source venv/bin/activate
python3 PROFIT-seq.py --mm_idx path_to_genome.mmi
```

```bash
Usage: PROFIT-seq.py [OPTIONS]

Options:
  --minknow_host TEXT    ip address for MinKNOW host. Defaults to to 127.0.0.1. 
  --minknow_port TEXT    port for MinKNOW service. Defaults to 8000.
  --guppy_address TEXT   address for guppy server. Defaults to ipc::///tmp/.guppy/5555.
  --guppy_config TEXT    guppy basecalling config. Defaults to dna_r9.4.1_450bps_fast.
  --dashboard_port TEXT  guppy basecalling config. Defaults to 55280.
  --mm_idx TEXT          Minimap2 index of reference sequences  [required]
  --version              Show the version and exit.
  --help                 Show this message and exit.
```

**NOTE: always use the fast basecalling config when running PROFIT-seq to avoid performance issues**

If everything works fine, the prompt of url for dashboard will appear on your screen.

### 2.3 Specify enrichment target or upload a toml formatted config on the web interface.

- Example of a valid TOML config

```toml
[[jobs]]
name = "Unblock_mt"
time = [0, 240]
ch = [129, 256]
bc = "all"
target = [
    {region = ["chrM", "0", "16569"], action="enrich"},
    {region = "multi", action="unblock"},
    {region = "miss", action="unblock"},
    {region = "unmapped", action="wait"},
]
```

- Usage:

```
- name: User-specified name for each target job.
- time: range of start and end time (minutes) for the job.
- ch: range of start and end channel for the job. (0-512 for MinION).
- bc: barcode for the job.
- target: specify what actions should be performed when aligning to specific region.
```

- Available barcode options:

```
# bc:
- barcode01,barcode02 (comma-seperated list of barcode names, only reads with these barcodes will be processed)
- classified (all reads with classified barcodes will be processed)
- unclassified (all reads with unclassified barcodes will be processed)
- all (all reads will be processed)
```

- Available target region options:

```
# target:
- chrom:start-end (reads that mapped to spefici region will be processed)
- multi (reads that are multi-mapped will be processed)
- mapped (all reads mapped to the reference index will be processed)
- miss (reads mapped to the reference index, but missed any target regions will be processed)
- unmapped (all reads that could not be aligned to the reference index will be processed)
- all (all reads with be processed)

Priority: all > unmapped > mapped > multi > region 
```

- Available action options:

```
# action
- stop_receiving (finish sequencing this read) 
- unblock (reject this read)
- wait (wait for decision in the next chunk) 
- balance (balance coverage for all target regions with action `balance`)

At lease one of the following combination of actions are required for a valid job
1. unmapped + mapped
2. unmapped + regions + miss
```

## 3.PROFIT-seq data analysis

Dependencies:
- Reference genome & annotation
- Guppy
- Porechop
- Minimap2
- samtools
- StringTie2
- gffread
- Salmon

### 3.1 Re-Basecall reads using hac/sup basecalling model using `minknow` user  

```bash
# e.g.
/opt/ont-guppy_5.1.15/bin/guppy_basecall_client \
    -r --input_path path_to_input --save_path path_to_output \
    -c dna_r9.4.1_450bps_hac.cfg -x auto \
    --port ///tmp/.guppy/5555 \
    --barcode_kits "EXP-NBD114" \
    --compress_fastq
```

### 3.2 Demultiplex reads acording to channel number

For instance, if you want to demultiplex reads from channel 1 to 256

```bash
python3 scripts/step1_demultiplex.py -i input.fastq.gz -o sample1.fastq.gz --start 1 --end 256
```

```
Usage: step1_demultiplex.py [OPTIONS]

Options:
  -i, --infile TEXT      input gzipped fastq file.  [required]
  -o, --outfile TEXT     output gzipped fastq file.  [required]
  -st, --start INTEGER   start channel number.  [required]
  -en, --end INTEGER     end channel number.  [required]
  -t, --threads INTEGER  number of threads.
```

### 3.3 Adapter trimming

Trim nanopore sequencing adapters using Porechop

```bash
porechop -i sample1.fastq.gz -o sample1.trimmed.fastq.gz --threads 32 --check_reads 1000
```

### 3.4 Consensus calling & orientation

Generate full-length consensus reads and non full-length reads 

```bash
python3 scripts/step2_consensus.py \
    -i sample1.trimmed.fastq.gz \
    -s sequencing_summary_FAQ85160_399ee876.txt \
    -o ./output_sample1 \
    -p sample1 \
    -t 16 \
    --trimA
```

```
Usage: step2_consensus.py [OPTIONS]

Options:
  -i, --input PATH       input trimmed fastq.  [required]
  -s, --summary PATH     input sequencing summary generate by MinKNOW.  [required]
  -o, --outdir PATH      output directory.  [required]
  -p, --prefix TEXT      output prefix name.  [required]
  -r, --adapter PATH     Adapter sequences file. Defaults to embedded splint adapter sequences.
  -t, --threads INTEGER  number of threads. Defaults to number of cpu cores.
  --trimA                trim 3' poly(A) sequences
``` 

- The output `sample1.fl.fa` is the full-length consensus reads with both 5' and 3' primers.
- The output `sample1.recovered.fa` is the partial fragents that only have one 5'/3' primer.

### 3.5 Isoform assembly & quantification

Use full-length reads for full-length isoform assembly and quantification

```bash
python step3_analysis.py \
    -i ./output_sample1 \
    -p sample1 \
    -r GRCh38.primary_assembly.genome.fa \
    -a gencode.v37.annotation.gtf \
    -t 16 \
    --assemble \
    --bed ../cancer_panel.bed
```

```
Usage: step3_analysis.py [OPTIONS]

Options:
  -i, --workspace PATH   directory of step2_consensus.py output  [required]
  -p, --prefix TEXT      sample prefix for step2_consensus.py  [required]
  -r, --genome PATH      reference genome fasta.  [required]
  -a, --gtf PATH         gene annotation gtf.  [required]
  -b, --bed PATH         bed file for target regions.  [required]
  -t, --threads INTEGER  number of threads. Defaults to number of cpu cores.
  --assemble             perform transcript isoform assemble.
  --help                 Show this message and exit.
``` 

- The output `sample1_isoforms.gtf` is the assembled and annotated full-length transcript isoforms in GTF format (requires `--assemble`).
- The output `sample1_isoforms.transcripts.sf` is transcript-level quantification results in the Salmon tsv format. The output contains five columns: transcript name, transcript length, effective length, CPM, number of reads.
- The output `sample1_isoforms.genes.sf` is gene-level quantification results. The output contains three columns: gene name, CPM, number of reads.
