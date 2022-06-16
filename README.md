# full-spectrum-sequencing
Scripts for full spectrum transcriptome sequencing

# Installation

Step 1. Clone this repository to your computer

```bash
git clone https://github.com/bioinfo-biols/full-spectrum-sequencing.git
cd full-spectrum-sequencing
```

Step 2. Create a virtual environment for running these scripts
```bash
virtualenv venv
source ./venv/bin/activate
```

Step 3. Install the modified version of ssw
```bash
# Install pyssw
git clone https://github.com/Kevinzjy/pyssw.git
cd pyssw & make 

# Install Porechop
git clone https://github.com/artic-network/Porechop.git
cd Porechop && python setup.py install
```

# Usage

```bash
# Trim adapters
porechop \
  -i sample.raw.fastq.gz \
  -o sample.trimmed.fastq.gz \
  --threads 16

# Detect circular consensus reads
ccs \
  -i sample.trimmed.fastq.gz \
  -o sample.ccs.fa \
  -r sample.ccs.raw.fa \
  -t 16

# Generate consensus sequences
python /path/to/consensus.py \
  -i sample.ccs.fa \
  -o sample.clean.fa \
  -s sample.qc_summary.txt \
  -t 16
```