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
git clone https://github.com/bioinfo-biols/PROFIT-seq.git
cd PROFIT-seq
```

- Create a virtual environment for running these scripts

```bash
virtualenv venv
source ./venv/bin/activate
```

- Install requirements & the corresponding version of `ont-pyguppy-client-lib`

```bash
pip install -r requirements.txt
pip install ont-pyguppy-client-lib==5.1.15
```

- (Optional) Install requirements for analysis PROFIT-seq data

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
- Minimap2
- samtools
- StringTie2
- gffcompare
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
python3 scripts/step1_demultiplex.py -i input.fastq.gz -o sample1.fastq.gz --start 0 --end 256
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
# e.g.
porechop -i sample1.fastq.gz -o sample1.trimmed.fastq.gz --threads 32 --check_reads 1000
```

### 3.4 Consensus calling & orientation

Generate full-length consensus reads and non full-length reads 

```bash
python3 scripts/step2_consensus.py \
    -i sample1.trimmed.fastq.gz \
    -c sample1.trimmed.fl_reads.fa \
    -n sample1.trimmed.nonfl_reads.fa \
    -s sample1.summary.txt \
    -t 32
```

```
Usage: step2_consensus.py [OPTIONS]

Options:
  -i, --input TEXT       input trimmed fastq.  [required]
  -c, --ccs TEXT         output full-length consensus reads.  [required]
  -n, --non TEXT         output non full-length reads.  [required]
  -s, --summary TEXT     consensus calling summary file.  [required]
  -r, --adapter TEXT     Adapter sequences file. Defaults to an embedded PROFIT-seq adapter fasta.
  -t, --threads INTEGER  number of threads.
```

### 3.5 (Optional) Trim poly(A) sequences

Trim poly(A) tails for full-length consensus reads to increase alignment accuracy

```bash
python scripts/step3_trim_polyA.py \
    -i sample1.trimmed.fl_reads.fa \
    -o sample1.trimmed.fl_reads.nopA.fa 
```

```
Usage: step3_trim_polyA.py [OPTIONS]

Options:
  -i, --input TEXT       input full-length consensus reads.  [required]
  -o, --output TEXT      output trimmed reads.  [required]
```

### 3.6 Isoform assembly

Use full-length reads for full-length isoform assembly with StringTie2

```bash
# Align to reference genome
minimap2 -t 32 -a -x splice GRCh38.primary_assembly.genome.fa sample1.trimmed.fl_reads.nopA.fa  \
    | samtools sort -@ 32 -o sample1.sorted.bam -
samtools index -@ 32 sample1.sorted.bam

# Full-length isoform assembly
stringtie sample1.sorted.bam \
    -p 32 -G gencode.v37.annotation.gtf \
    -t -L -c 1.5 -s 1 -g 0 -f 0.05 \
    -o sample1_out.gtf \
    -A sample1_genes.list
    
# Annotated StringTie2 output
gffcompare -o sample1_gencode -r gencode.v37.annotation.gtf sample1_out.gtf

# Generate annotated isoform sequences
gffread sample1_gencode.annotated.gtf -g GRCh38.primary_assembly.genome.fa -w sample1_gencode.annotated.fa
```

The output `sample1_gencode.annotated.fa` is the assembled and annotated full-length transcript isoform sequences.

Detailed annotation information is in `sample1_gencode.annotated.gtf`

### 3.7 Quantification

PROFIT-seq uses a hybrid-quantification strategy to combine full-length consensus and non-fl fragments.

**NOTE: If you do not want to perform transcript assembly, replace `sample1_gencode.annotated.fa` with the reference transcriptome fasta instead** 

```bash
# Full-length consensus reads
minimap2 -ax map-ont -t 32 -p 1.0 -N 100 sample1_gencode.annotated.fa sample1.trimmed.fl_reads.nopA.fa \
    | samtools sort -@ 32 -o sample1_transcripts.fl.sorted.bam -
samtools index sample1_transcripts.fl.sorted.bam
salmon quant --noErrorModel --noLengthCorrection -p 32 -l U \
    -t sample1_gencode.annotated.fa -a sample1_transcripts.fl.sorted.bam \
    -o output/fl_reads 
    
# Recover non-fl partial fragments 
minimap2 -ax map-ont -t 32 -p 1.0 -N 100 sample1_gencode.annotated.fa sample1.trimmed.nonfl_reads.fa  \
    | samtools sort -@ 32 -o sample1_transcripts.nonfl.sorted.bam -
samtools index sample1_transcripts.nonfl.sorted.bam
salmon quant --noErrorModel --noLengthCorrection -p 32 -l U \
    -t sample1_gencode.annotated.fa -a sample1_transcripts.nonfl.sorted.bam \
    -o output/nonfl_reads 
```

Integration of hybrid quantification results:

```
python scripts/step4_quantification.py \
    --fl path/to/fl_reads/quant.sf \
    --nonfl path/to/nonfl_reads/quant.sf \
    --output merged.sf
```

```
Usage: step4_quantification.py [OPTIONS]

Options:
  --fl TEXT              quant.sf for full-length reads.  [required]
  --nonfl TEXT           quant.sf for non-fl reads.  [required]
  --output TXT           output quantification result.  [required]
```

The final output contains four columns: transcript name, transcript length, CPM, number of reads