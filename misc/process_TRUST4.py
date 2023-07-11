from biofrost.env import *


def load_sample(sample_file, reads_file):
    # Load reads 
    total_reads = pd.read_csv(reads_file, index_col=0)
    total_reads.index = total_reads.index.map(lambda x: x.split('/')[-1].split(".")[0])

    sample_df = pd.read_csv(sample_file, sep='\s+', engine='python', header=None)
    sample_df.columns = ['Sample', 'Group']
    sample_df = sample_df.set_index('Sample')
    sample_df['Reads'] = total_reads['Reads']
    
    return sample_df

def load_trust4_result(trust4_dir, sample_df):
    bcr_reads, tcr_reads = {}, {}
    for sample, row in sample_df.iterrows():
        clinic = row['Group']

        # read the individual sample result
        cdr3 = pd.read_csv(trust4_dir / sample / f"{sample}_report.tsv", sep="\t")
        cdr3 = cdr3.rename({"#count": "count"}, axis="columns")

        # filter out count of cdr ==0 and add the phenotype info for downstream comparison
        cdr3 = cdr3.loc[cdr3['count'] > 0]
        cdr3['sample'] = sample
        cdr3['clinic'] = clinic

        # determine whether the cdr animo acid is complete or partial
        cdr3['is_complete'] = cdr3['CDR3aa'].map(lambda x: "N" if x == "partial" or x == "out_of_frame" or x[0] == "_" or x[0] == "\\" else "Y")

        # exact the TCR and BCR
        cdr3_bcr = cdr3.loc[cdr3.apply(lambda x: x["V"].startswith("IG") or x["J"].startswith("IG") or x["C"].startswith("IG"), axis=1)].copy()
        cdr3_tcr = cdr3.loc[cdr3.apply(lambda x: x["V"].startswith("TR") or x["J"].startswith("TR") or x["C"].startswith("TR"), axis=1)].copy()

        if cdr3_bcr.shape[0] != 0:
            # add lib size and clinic traits
            cdr3_bcr.loc[:, 'lib_size'] = cdr3_bcr['count'].sum()

            # split BCR into heavy chain and light chain
            cdr3_bcr_heavy = cdr3_bcr.loc[cdr3_bcr.apply(lambda x: x["V"].startswith("IGH") or x["J"].startswith("IGH") or x["C"].startswith("IGH"), axis=1)].copy()
            cdr3_bcr_light = cdr3_bcr.loc[cdr3_bcr.apply(lambda x: x["V"].startswith(("IGL", "IGK")) or x["J"].startswith(("IGL", "IGK")) or x["C"].startswith(("IGL", "IGK")), axis=1)].copy()

            #save BCR and TCR info for downsteam use
            cdr3_bcr_light.to_csv(trust4_dir / sample / f"{sample}_BCR_light.tsv", sep="\t", index=False)
            cdr3_bcr_heavy.to_csv(trust4_dir / sample / f"{sample}_BCR_heavy.tsv", sep="\t", index=False)

        if cdr3_tcr.shape[0] != 0:
            cdr3_tcr.loc[:, 'lib_size'] = cdr3_tcr['count'].sum()
            cdr3_tcr.to_csv(trust4_dir / sample / f"{sample}_TCR.tsv", sep="\t", index=False)

        bcr_reads[sample] = cdr3_bcr['count'].sum()
        tcr_reads[sample] = cdr3_tcr['count'].sum()

    cdr3_reads_fraction = pd.DataFrame({"BCR": bcr_reads, "TCR": tcr_reads})
    cdr3_reads_fraction.to_csv(trust4_dir / "cdr3_reads.tsv", sep="\t", index=True)

# Input files
prefix = "output_FS_merged"

reads_file = Path(f'./06.repertoire/{prefix}.csv')
sample_file = Path('./01.reads/sample.txt')
trust4_dir = Path(f"./06.repertoire/{prefix}")

# Run
sample_df = load_sample(sample_file, reads_file)

lut = dict(zip(['I', 'L', 'H'], ['#66bd63', '#fdae61', '#d6604d']))
weight = dict(zip(['I', 'L', 'H'], [0, 1, 2]))

load_trust4_result(trust4_dir, sample_df)