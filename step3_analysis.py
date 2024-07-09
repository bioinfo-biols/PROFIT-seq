#!/usr/bin/env python3
import os
import sys
import subprocess
import click
import pandas as pd
from collections import defaultdict
from pathlib import Path
from edlib import align
from PROFIT_seq.logger import get_logger
from PROFIT_seq.parser import yield_gff
from PROFIT_seq.seqIO import yield_fastx, revcomp
from PROFIT_seq.utils import load_toml, check_dir
LOGGER = get_logger("PROFIT-seq", debugging=False)


def check_dependencies(toolset):
    """
    Check required softwares
    """
    for tool in toolset:
        status, output = subprocess.getstatusoutput(f"which {tool}")
        if status != 0:
            sys.exit(f"{tool} command not found")
    return 0


def load_step2_summary(summary_file):
    """
    Load consensus output summary
    """
    summary = pd.read_csv(summary_file, index_col=0, sep="\t")
    blacklist = summary.index[summary['Blacklist'] is False]
    return blacklist


def dedup_reads(input_fa, out_dir, prefix, log_file, threads):
    out_paf = f"{out_dir}/{prefix}.all_vs_all.paf"
    out_fa = f"{out_dir}/{prefix}.dedup.fa"
    LOGGER.info(f"Perform minimap2 alignment")
    mm2_cmd = f"minimap2 -t {threads} -x ava-ont {input_fa} {input_fa} > {out_paf}"
    run_cmd(mm2_cmd, log_file)

    # Get reads cluster
    LOGGER.info(f"Loading all vs all paf")
    reads_cluster = defaultdict(dict)
    error_rate = 0.25
    with open(out_paf, "r") as f:
        for line in f:
            content = line.rstrip().split('\t')
            q_name, q_len, q_st, q_en, strand, t_name, t_len, t_st, t_en = content[:9]
            if (int(q_en) - int(q_st)) / int(q_len) < 0.8:
                continue
            if (int(t_en) - int(t_st)) / int(t_len) < 0.8:
                continue
            umi1 = q_name.split("|")[1]
            umi2 = t_name.split("|")[1]
            base1 = len(umi1) - sum([i == ';' for i in umi1])
            base2 = len(umi2) - sum([i == ';' for i in umi2])
            mm = min(base1, base2) * error_rate

            if strand == "+" and align(umi1, umi2)['editDistance'] > mm:
                continue
            if strand == "-" and align(umi1, revcomp(umi2))['editDistance'] > mm:
                continue
            reads_cluster[q_name][t_name] = 1

    # Get duplicate reads
    LOGGER.info(f"Extracting RCA duplicates")
    umi_to_read, read_to_umi = defaultdict(list), {}
    for umi_idx, q_name in enumerate(reads_cluster):
        umi_reads = [q_name, ] + list(reads_cluster[q_name])

        assigned = []
        for read_id in umi_reads:
            if read_id in read_to_umi:
                assigned.append(read_to_umi[read_id])

        if len(assigned) == 0:
            for read_id in umi_reads:
                read_to_umi[read_id] = umi_idx
            umi_to_read[umi_idx] = umi_reads

        else:
            for read_id in umi_reads:
                read_to_umi[read_id] = assigned[0]
            umi_to_read[assigned[0]] += umi_reads

            if len(assigned) > 1:
                for _idx in assigned[1:]:
                    for read_id in umi_to_read[_idx]:
                        read_to_umi[read_id] = assigned[0]
                    umi_to_read[assigned[0]] += umi_to_read[_idx]

    _clusters = set([x for _, x in read_to_umi.items()])
    dup_reads = set(sum([umi_to_read[_cluster][1:] for _cluster in _clusters], []))
    LOGGER.info(f"RCA duplicate reads: {len(dup_reads)}")

    # Filter fastq
    with open(input_fa, 'r') as f, open(out_fa, 'w') as out:
        for line in f:
            seq = f.readline()
            read_id = line.rstrip().split(' ')[0].lstrip(">")
            if read_id in dup_reads:
                continue
            out.write(line + seq)

    return out_fa


def extract_hifl_fasta(ccs_fa, blacklist, out_fa):
    """
    Filter high-confidence full-length consensus reads
    """
    with open(out_fa, 'w') as out:
        for seq_id, seq, sep, qual in yield_fastx(ccs_fa):
            if seq_id in blacklist:
                continue
            out.write(f">{seq_id}\n{seq}\n")
    return 0


def run_cmd(cmd, log_file):
    """
    Run command and redirect output to log file
    """
    with open(log_file, 'a') as log:
        status, output = subprocess.getstatusoutput(cmd)
        log.write(output)
        if status != 0:
            sys.exit(f"Failed to run command: {cmd}")
    return status


def run_stringtie(input_fa, genome_fa, gene_gtf, out_dir, prefix, isoform_gtf, log_file, threads):
    """
    Run stringtie2 pipeline for transcript isoform assembly
    """
    sorted_bam = f"{out_dir}/{prefix}_genome.sorted.bam"
    assembled_gtf = f"{out_dir}/{prefix}_out.gtf"

    mm2_cmd = f"minimap2 -t {threads} -ax splice {genome_fa} {input_fa} | samtools sort -@ {threads} -o {sorted_bam}"
    index_cmd = f"samtools index -@ {threads} {sorted_bam}"
    stringtie_cmd = f"stringtie {sorted_bam} -p {threads} -G {gene_gtf} -t -L -c 1.5 -s 1 -g 0 -f 0.05 -o {assembled_gtf} -A {out_dir}/{prefix}_genes.list"
    merge_cmd = f"stringtie --merge -G {gene_gtf} {assembled_gtf} > {isoform_gtf}"

    run_cmd(mm2_cmd, log_file)
    run_cmd(index_cmd, log_file)
    run_cmd(stringtie_cmd, log_file)
    run_cmd(merge_cmd, log_file)

    return 0


def extract_isoform_sequence(genome_fa, gtf_file, fa_file, log_file):
    """
    Extract isoform sequences
    """
    gffread_cmd = f"gffread {gtf_file} -g {genome_fa} -w {fa_file}"
    run_cmd(gffread_cmd, log_file)
    return fa_file


def merge_recovered_sequence(workspace, prefix, log_file):
    """
    Merge all partial recovered reads
    """
    merged_fa = f"{workspace}/tmp/{prefix}_merged.fa"
    merge_cmd = f"cat {workspace}/{prefix}.fl.fa {workspace}/{prefix}.recovered.fa > {merged_fa}"
    run_cmd(merge_cmd, log_file)
    return merged_fa


def run_salmon(input_fa, reference_fa, out_dir, prefix, log_file, threads):
    """
    Run salmon for primitive quantification
    """
    sorted_bam = f"{out_dir}/{prefix}.sorted.bam"
    mm2_cmd = f"minimap2 -ax map-ont -t {threads} -p 1.0 -N 100 {reference_fa} {input_fa} | samtools sort -@ {threads} -o {sorted_bam} - "
    index_cmd = f"samtools index {sorted_bam}"
    salmon_cmd = f"salmon quant --noErrorModel --noLengthCorrection -p {threads} -l U -t {reference_fa} -a {sorted_bam} -o {out_dir}/{prefix}"

    run_cmd(mm2_cmd, log_file)
    run_cmd(index_cmd, log_file)
    run_cmd(salmon_cmd, log_file)

    return f"{out_dir}/{prefix}/quant.sf"


def get_targets(assembled_gtf, bed_file):
    """
    Get assembled transcripts within target regions
    """
    # Load bed file
    target_index = defaultdict(dict)
    n = 0
    with open(bed_file, 'r') as f:
        for line in f:
            content = line.strip().split('\t')
            chrom, start, end = content[0], int(content[1]), int(content[2])
            div_start, div_end = int(start/500), int(end/500)
            for i in range(div_start, div_end+1):
                target_index[chrom].setdefault(i, []).append((start, end))
            n += 1
    LOGGER.info(f"Loaded {n} target regions")

    # Parse gtf
    target_transcripts, tscp_to_gene = {}, {}
    for parser in yield_gff(assembled_gtf, is_gtf=True):
        if parser.type != 'transcript':
            continue

        transcript_id = parser.transcript_id
        gene_id = parser.gene_id
        tscp_to_gene[transcript_id] = gene_id

        if parser.contig not in target_index:
            continue
        div_start = int(parser.start / 500)
        div_end = int(parser.end / 500)
        for i in range(div_start, div_end+1):
            if i not in target_index[parser.contig]:
                continue
            for (st, en) in target_index[parser.contig][i]:
                if en < parser.start or parser.end < st:
                    continue
                target_transcripts[transcript_id] = 1
                break

            if transcript_id in target_transcripts:
                break
    LOGGER.info(f"Loaded {len(target_transcripts)} target transcripts")

    return target_transcripts, tscp_to_gene


def correct_expression_level(fl_sf, recovered_sf, target_transcripts, tscp_to_gene, transcript_sf, gene_sf):
    """
    A simplified version of quantification algorithm
    """
    # Load salmon results
    col = "NumReads"
    fl_df = pd.read_csv(fl_sf, sep="\t", index_col=0)
    nonfl_df = pd.read_csv(recovered_sf, sep="\t", index_col=0)

    # Merge dataframe
    merged_reads = nonfl_df[col]
    merged_sf = pd.DataFrame({
        "Length": fl_df['Length'],
        "EffectiveLength": fl_df['EffectiveLength'],
        "CPM": 1000 * 1000 * merged_reads / merged_reads.sum(),
        "NumReads": merged_reads}
    )

    # Correction of target reads
    targets = merged_sf.index[merged_sf.index.isin(target_transcripts)]
    s_i = fl_df.loc[targets, col]
    p_i = s_i / s_i.sum()

    b_i = nonfl_df.loc[targets, "NumReads"]
    merged_sf.loc[targets, "NumReads"] = b_i.sum() * p_i
    merged_sf.loc[targets, "CPM"] = 1000 * 1000 * merged_sf['NumReads'] / merged_sf['NumReads']

    merged_sf.to_csv(transcript_sf, sep="\t", index=True, index_label="Name")

    # Gene-level expression
    gene_cpm = merged_sf.groupby(merged_sf.index.map(tscp_to_gene)).apply(lambda df: df['CPM'].sum())
    gene_reads = merged_sf.groupby(merged_sf.index.map(tscp_to_gene)).apply(lambda df: df['NumReads'].sum())
    gene_df = pd.DataFrame({
        "CPM": gene_cpm,
        "NumReads": gene_reads
    })
    gene_df.to_csv(gene_sf, sep="\t", index=True, index_label="Name")
    return 0



@click.command()
@click.option('--workspace', '-i', type=Path, required=True,
              help="directory of step2_consensus.py output")
@click.option('--prefix', '-p', type=str, required=True,
              help="sample prefix for step2_consensus.py")
@click.option('--genome', '-r', type=Path, required=True,
              help="reference genome fasta.")
@click.option('--gtf', '-a', type=Path, required=True,
              help="gene annotation gtf.")
@click.option('--bed', '-b', type=Path, required=True,
              help="bed file for target regions.")
@click.option('--threads', '-t', type=int, default=os.cpu_count(),
              help='number of threads. Defaults to number of cpu cores.')
@click.option('--assemble', is_flag=True,
              help="perform transcript isoform assemble.")
def main(workspace, prefix, genome, gtf, bed, threads, assemble):
    # Check software dependencies
    LOGGER.info("Checking dependencies")
    tool_chains = ["minimap2", "samtools", "salmon"]
    check_dependencies(tool_chains)

    assemble_tools = ["stringtie", "gffread"]
    if assemble:
        check_dependencies(assemble_tools)

    # Load step2 result
    LOGGER.info('Loading output from step2_consensus.py')
    summary_file = workspace / f"{prefix}.consensus_summary.txt"
    blacklist = load_step2_summary(summary_file)

    # Prepare files
    tmp_dir = check_dir(workspace / "tmp")
    log_file = workspace / f"{prefix}_analysis.log"

    # Extract high-confidence reads
    ccs_fa = workspace / f"{prefix}.fl.fa"
    dedup_fa = dedup_reads(ccs_fa, tmp_dir, prefix, log_file, threads)
    hifl_fa = tmp_dir / f"{prefix}.hifl.fa"
    extract_hifl_fasta(dedup_fa, blacklist, hifl_fa)

    # Run stringtie sequence assembly
    if assemble:
        LOGGER.info('1/3 Running isoform assembly pipeline')
        isoform_gtf = workspace / f"{prefix}_isoforms.gtf"
        run_stringtie(hifl_fa, genome, gtf, tmp_dir, prefix, isoform_gtf, log_file, threads)
    else:
        LOGGER.info('1/3 Skipped isoform assembly')
        isoform_gtf = gtf

    # Get transcriptome sequences
    isoform_fa = tmp_dir / f"{prefix}_isoforms.fa"
    extract_isoform_sequence(genome, isoform_gtf, isoform_fa, log_file)

    # Merged fa
    LOGGER.info('2/3 Running salmon for primitive quantification')
    merged_fa = merge_recovered_sequence(workspace, prefix, log_file)
    hifl_sf = run_salmon(hifl_fa, isoform_fa, tmp_dir, f"{prefix}_hifl", log_file, threads)
    merged_sf = run_salmon(merged_fa, isoform_fa, tmp_dir, f"{prefix}_merged", log_file, threads)

    # Merge salmon results
    LOGGER.info('3/3 Calculate corrected expression levels')
    target_transcripts, tscp_to_gene = get_targets(isoform_gtf, bed)
    transcript_sf = f"{workspace}/{prefix}_isoforms.transcripts.sf"
    gene_sf = f"{workspace}/{prefix}_isoforms.genes.sf"
    correct_expression_level(hifl_sf, merged_sf, target_transcripts, tscp_to_gene, transcript_sf, gene_sf)

    LOGGER.info('Finished all analysis.')


if __name__ == '__main__':
    main()
