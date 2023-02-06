"""Functions for manipulating annotation files"""
import re

GFF_gene_id = re.compile(r'gene_id=(\S+);')
GFF_tscp_id = re.compile(r'transcript_id=(\S+);')
GFF_gene_type = re.compile(r'gene_type=(\S+);')
GFF_gene_name = re.compile(r'gene_name=(\S+);')
GFF_transcript_type = re.compile(r'transcript_type=(\S+);')
GFF_transcript_name = re.compile(r'transcript_name=(\S+);')


class GFFParser(object):
    """
    A class to parse GENCODE GTF record
    Parameters
    ----------
    content : list
        List of gtf records.
    Attributes
    ----------
    contig : str
        Chromosome name, chr[1-22,X,Y,M] or GRC accession
    source : str
        Annotation source, {ENSEMBL, HAVANA}
    type : str
        Feature type, {gene, transcript, exon, CDS, UTR, start_codon, stop_codon, Selenocysteine}
    start : int
        Genomic start location, integer-value (1-based)
    end : int
        Genomic end location, integer-value
    score : str
        Score (not used), {.}
    strand : strand
        Genomic strand, {+,-}
    phase : str
        Genomic phase (for CDS features), {0, 1, 2}
    attr_string : str
        Additional information as key-value pairs
    attr : dict
        Dict of additional information
    gene_id : str
        Gene id
    transcript_id : str
        Transcript id
    gene_type : str
        Gene type
    gene_name : str
        Gene name
    transcript_type : str
        Transcript type
    transcript_name : str
        Transcript name
    """

    def __init__(self, content):
        self.contig = content[0]
        self.source = content[1]
        self.type = content[2]
        self.start, self.end = int(content[3]), int(content[4])
        self.score = content[5]
        self.strand = content[6]
        self.phase = content[7]
        self.attr_string = content[8]

    @property
    def attr(self):
        """Get additional information"""
        field = {}
        for attr_values in [re.split(r'=', i.strip()) for i in self.attr_string.split(';')]:
            key, value = attr_values[0], attr_values[1:]
            field[key] = ' '.join(value).strip('"')
        return field

    @property
    def gene_id(self):
        """Get gene id"""
        match = re.findall(GFF_gene_id, self.attr_string)
        return match[0] if match else None

    @property
    def transcript_id(self):
        """Get transcript id"""
        match = re.findall(GFF_tscp_id, self.attr_string)
        return match[0] if match else None

    @property
    def gene_type(self):
        """Get gene_type"""
        match = re.findall(GFF_gene_type, self.attr_string)
        return match[0] if match else None

    @property
    def gene_name(self):
        """Get gene name"""
        match = re.findall(GFF_gene_name, self.attr_string)
        return match[0] if match else None

    @property
    def transcript_type(self):
        """Get transcript type"""
        match = re.findall(GFF_transcript_type, self.attr_string)
        return match[0] if match else None

    @property
    def transcript_name(self):
        """Get transcript name"""
        match = re.findall(GFF_transcript_name, self.attr_string)
        return match[0] if match else None

    def __repr__(self):
        return '{} {}:{}-{}:{}'.format(self.type, self.contig, self.start, self.end, self.strand)


GTF_gene_id = re.compile(r'gene_id "(\S+)";')
GTF_tscp_id = re.compile(r'transcript_id "(\S+)";')
GTF_gene_type = re.compile(r'gene_type "(\S+)";')
GTF_gene_name = re.compile(r'gene_name "(\S+)";')
GTF_transcript_type = re.compile(r'transcript_type "(\S+)";')
GTF_transcript_name = re.compile(r'transcript_name "(\S+)";')


class GTFParser(GFFParser):
    @property
    def attr(self):
        field = {}
        for attr_values in [re.split(r'\s+', i.strip()) for i in self.attr_string.split(';')[:-1]]:
            key, value = attr_values[0], attr_values[1:]
            field[key] = ' '.join(value).strip('"')
        return field

    @property
    def gene_id(self):
        match = re.findall(GTF_gene_id, self.attr_string)
        return match[0] if match else None

    @property
    def transcript_id(self):
        match = re.findall(GTF_tscp_id, self.attr_string)
        return match[0] if match else None

    @property
    def gene_type(self):
        match = re.findall(GTF_gene_type, self.attr_string)
        return match[0] if match else None

    @property
    def gene_name(self):
        match = re.findall(GTF_gene_name, self.attr_string)
        return match[0] if match else None

    @property
    def transcript_type(self):
        match = re.findall(GTF_transcript_type, self.attr_string)
        return match[0] if match else None

    @property
    def transcript_name(self):
        match = re.findall(GTF_transcript_name, self.attr_string)
        return match[0] if match else None


def yield_gff(gff_file, is_gtf=True):
    """Iter through GTF/GFF records
    Args:
        gff_file (str): genome annotation file in GTF/GFF format
        is_gtf (bool, optional): Whether input is in GTF format. Defaults to True.
    Yields:
        GTF/GFFParser: GTF/GFFParser object
    Examples
    --------
    >>> # Load GTF file
    >>> gtf_file = "/data/public/database/GENCODE/Human_Release_37/chrM.gtf"
    >>> print(len([record for record in yield_gff(gtf_file)]))
    >>> # Load GFF3 file
    >>> gff_file = "/data/public/database/GENCODE/Human_Release_37/chrM.gff3"
    >>> print(len([record for record in yield_gff(gff_file, is_gtf=False)]))
    """
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            content = line.rstrip().split('\t')
            if is_gtf:
                parser = GTFParser(content)
            else:
                parser = GFFParser(content)
            yield parser