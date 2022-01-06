# Python 2 & 3 Compatibility
from __future__ import print_function, unicode_literals, absolute_import
# Dependencies
import sys, argparse, time, multiprocessing, os, re, gzip
from collections import *
import pandas as pd
import numpy as np
from Bio import SeqIO

# Optional
try:
    from tqdm import tqdm
    tqdm_available = True
except ModuleNotFoundError:
    tqdm_available = False

"""
==================================================================================================================
NIF Finder v2022.01.05
______________________
Created by Josh L. Espinoza (08.28.2017)
==================================================================================================================
Discover features in putative nonsense-mediated inducing factors.
==================================================================================================================

kozak_type   [Default: 1]
             0: `ATG`
             1: `[A|G]..ATG or ATG[G]
             2: `[A|G]..ATG[G]
             3: `[A|G].[C]ATG[G]
"""

# Objects & Functions
# ==================================================================================================================
# Genetic Code
# ==================================================================================================================
D_codon_aa = {
     "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "TCT": "S",
     "TCC": "S", "TCA": "S", "TCG": "S", "TAT": "Y", "TAC": "Y",
     "TAA": "*", "TAG": "*", "TGA": "*", "TGT": "C", "TGC": "C",
     "TGG": "W", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
     "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P", "CAT": "H",
     "CAC": "H", "CAA": "Q", "CAG": "Q", "CGT": "R", "CGC": "R",
     "CGA": "R", "CGG": "R", "ATT": "I", "ATC": "I", "ATA": "I",
     "ATG": "M", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
     "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGT": "S",
     "AGC": "S", "AGA": "R", "AGG": "R", "GTT": "V", "GTC": "V",
     "GTA": "V", "GTG": "V", "GCT": "A", "GCC": "A", "GCA": "A",
     "GCG": "A", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
     "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"}


# ==================================================================================================================
# Functions
# ==================================================================================================================
# Translate DNA to Amino Acid
def translate(seq):
    seq = seq.upper()
    peptide= ""
    for codon in [seq[k:k+3] for k in range(0,len(seq),3)]:
        if len(codon) == 3:
            peptide += D_codon_aa[codon]
    return peptide

# GC-Content
def GC(seq):
    seq = seq.upper()
    try:
        return float(seq.count("G") + seq.count("C"))/len(seq)
    except ZeroDivisionError:
        return np.nan

# Open Reading Frames
def orf_scan(id, seq, kozak_type=1, tol_min_codons=10, max=False, file_errors=sys.stderr):
    """
    A function that scans all reading frames of a DNA sequence string and outputs either a list of open reading frames [ORF]
    or a single open reading frame (maximum length)
    Parameters:
        seq = DNA sequence string input
        kozak_type   [Default: 1]
                     0: `ATG`
                     1: `[A|G]..ATG or ATG[G]
                     2: `[A|G]..ATG[G]
                     3: `[A|G].[C]ATG[G]
        min_c = Minimum codon threshold for protein generation. Does not include STOP codon.
        (e.g. MAAAT* would be 5 codons, shorter than the default 10 codons and would not be in output.)
        max = If max is True, then only the maximum ORF is returned, if False then a list of all ORFs are output
    """
    seq = seq.upper()

    # Define ORFs
    stop_codons= ["TAA","TGA","TAG"]
    start_positions = findall_positions("ATG", seq)

    try:
        # Find ORFs
        ORFs = list()
        if len(start_positions) > 0:
            polypeptide = list()
            for start_site in start_positions:
                # Kozak sequence threshold
                if kozak_type == 0: kozak_condition = True
                if kozak_type == 1: kozak_condition = (seq[start_site - 3] in list("AG")) or (seq[start_site + 3] == "G")
                if kozak_type == 2: kozak_condition = (seq[start_site - 3] in list("AG")) and (seq[start_site + 3] == "G")
                if kozak_type == 3: kozak_condition = (seq[start_site - 3] in list("AG")) and (seq[start_site + 3] == "G") and (seq[start_site - 1] == "C")

                if kozak_condition == True:
                    # Reading Frames
                    triplets= [seq[i:i+3] for i in range(start_site,len(seq),3)]
                    for codon in triplets:
                        polypeptide.append(codon)
                        if codon in stop_codons:
                            break

                    # Check for stop codons
                    if polypeptide[-1] in stop_codons:
                        if (len(polypeptide) - 1) >= tol_min_codons:
                            ORFs.append("".join(polypeptide))
                # Reset
                polypeptide = list()
        # Output
        if max == True:
            if len(ORFs) > 0:
                return [max(ORFs, key=len)]
            else:
                return []
        if max == False:
            return ORFs
    except IndexError:
        print("%s contains non-`ATG` start codon"%id, file=file_errors)
        return []

#Shannon entropy
def entropy(seq):
    seq = seq.upper()
    entropy_units = []
    N = len(seq)
    for x, count in Counter(seq).items():
        P_x = float(count)/N
        entropy_units.append(P_x * np.log2(P_x))
    return -sum(entropy_units)

# Findall Positions
def findall_positions(pattern, sequence, relative=False):
    """
    re.findall w/ overlapping positions
    """
    seq_length = len(sequence)
    positions = list()
    while True:
        match = re.search(pattern, sequence)
        if match is not None:
            current = match.start()
            sequence = sequence[current + 1:]
            if len(positions) == 0:
                positions.append(current)
            else:
                positions.append(current  + positions[-1] + 1)
        if match is None:
            break
    Ar_positions = np.array(positions)
    if relative == False: return Ar_positions
    if relative == True: return Ar_positions/seq_length

# ==================================================================================================================
# Objects
# ==================================================================================================================
class mRNA(object):
    """
    mRNA object that contains all of the sequences, computes lengths/gc, identifies uORFs, and organizes the data into a useable format
    """
    def __init__(self, id, metadata=dict()):
        self.id = id
        self.sequences = dict()
        self.metadata = metadata
        self.errors = list()
        self.fatal_error = False

        # Default NIFs
        self.nif_uorf = np.nan
        self.nif_long3utr = np.nan
        self.nif_dej = np.nan
        self.nif_termination_gc = np.nan

        # uORFs
        self.uorf_fna = OrderedDict()
        self.uorf_faa = OrderedDict()

     # In-silico splicing
    def splice(self):
        """
        Identifies if 5'UTR, CDS, and 3'UTRs are available.
        Splits sequence by exon-junctions for splice-sites
        """
        # Exons
        exons = re.split("[a-z]+", self.sequences["pre-mRNA"])
        self.num_exons = len(exons)
        self.sequences["exons"] = exons
        self.sequences["mRNA"] = "".join(exons)
        # Introns
        introns = re.split("[A-Z]+", self.sequences["pre-mRNA"])[1:-1]
        self.sequences["introns"] = introns

        # 5'UTR, 3'UTR
        utr5_proceed = True
        utr3_proceed = True
        if self.sequences["mRNA"].startswith(self.sequences["cds"]):
            self.errors.append("Warning: Missing `5'UTR'")
            utr5_proceed = False
        if self.sequences["mRNA"].endswith(self.sequences["cds"]):
            self.errors.append("Warning: Missing `3'UTR'")
            utr3_proceed = False
        if all([utr5_proceed, utr3_proceed]):
            if self.sequences["cds"] in self.sequences["mRNA"]:
                utr5, utr3 = re.split(self.sequences["cds"], self.sequences["mRNA"])
                self.sequences.update({"5'utr":utr5, "3'utr":utr3})
            else:
                self.errors.append("Fatal: `CDS` not in `mRNA`")
                self.fatal_error = True
        if all([utr5_proceed == False, utr3_proceed == False]):
            self.fatal_error = True
        if utr5_proceed and not utr3_proceed:
            self.sequences["5'utr"] = self.sequences["mRNA"][:len(self.sequences["cds"])]
        if utr3_proceed and not utr5_proceed:
            self.sequences["3'utr"] = self.sequences["mRNA"][len(self.sequences["cds"]):]

        # Exon Junctions
        if self.num_exons > 1:
            self.splice_sites = [re.search(exon, self.sequences["mRNA"]).start() for exon in self.sequences["exons"][1:]]
        else:
            self.splice_sites = list()

    # Upstream Open Reading Frames
    def uORFs(self, kozak_type, tol_min_codons, file_errors):
        """
        Identifies the upstream orfs while considering the kozak sequence
        """
        # 5'UTR
        if "5'utr" in self.sequences:
            seq_5utr = self.sequences["5'utr"]
            # Open Reading Frames
            seqs_uorfs_list = orf_scan(self.id, seq_5utr, kozak_type=kozak_type, tol_min_codons=tol_min_codons, file_errors=file_errors)


            if len(seqs_uorfs_list) > 0:
                for i, seq_uorf in enumerate(seqs_uorfs_list):
                    aa_seq = translate(seq_uorf)
                    # Headers
                    fna_header = "%s-%d len=%d;gc=%s;H=%s"%(self.id,i,len(seq_uorf), str(GC(seq_uorf)), str(entropy(seq_uorf)))
                    faa_header = "%s-%d len=%d;gc=%s;H=%s"%(self.id,i,len(seq_uorf), str(GC(seq_uorf)), str(entropy(aa_seq)))
                    # Fasta
                    self.uorf_fna[fna_header] = seq_uorf
                    self.uorf_faa[faa_header] = aa_seq
                # Get uORF data
                max_uorf = max(seqs_uorfs_list, key=len)
                uorf_content = [len(seqs_uorfs_list),  [*map(GC, seqs_uorfs_list)], [*map(len, seqs_uorfs_list)], len(max_uorf)]
                self.nif_uorf = True
            else:
                # No uORFs
                uorf_content = [0] + 3*[np.nan]
                self.nif_uorf = False
        else:
            # No 5'UTR
            uorf_content = 4*[np.nan]
        return pd.Series(uorf_content, index=["uORFs|__|count", "uORFs|__|gc-content", "uORFs|__|lengths", "uORFs|__|max-length"], name=self.id)

    # Region upstream of termination codon
    def termination_gc(self, tol_window):
        """
        Computes GC content for window upstream of termination codon
        """
        seq_query = self.sequences["cds"][-tol_window:]
        self.nif_termination_gc = GC(seq_query)
        return self.nif_termination_gc

    # Downstream Exon Junction
    def dEJs(self, tol_dej):
        """
        Compute downstream exon-junctions for transcripts with 3'UTRs
        """
        # No 3'UTR data
        if "3'utr" not in self.sequences:
            return np.nan
        else:
            # If there is a 5'UTR or not
            if "5'utr" in self.sequences:
                utr5_cds_size = len(self.sequences["5'utr"]) + len(self.sequences["cds"])
            else:
                utr5_cds_size = len(self.sequences["cds"])
            # If there are Splice Sites
            if len(self.splice_sites) > 0:
                if self.splice_sites[-1] >= (utr5_cds_size + tol_dej):
                    # Update NIF
                    self.nif_dej = True
                    query_dej = int(self.splice_sites[-1] - utr5_cds_size)
                    return query_dej

                else:
                    self.nif_dej = False
                    return np.nan
            # No Splice Sites
            else:
                return np.nan
    # Compute lengths for sequences
    def sequence_stats(self, tol_long3utr):
        """
        Calculates general stats like GC-content and Length for all of the sequences
        """
        # Sequences
        # 5'UTR
        if "5'utr" in self.sequences:
            seq_5utr = self.sequences["5'utr"]
            len_5utr = len(seq_5utr)
        else:
            seq_5utr = ""
            len_5utr = np.nan

        # 3'UTR
        if "3'utr" in self.sequences:
            seq_3utr = self.sequences["3'utr"]
            len_3utr = len(seq_3utr)

            # Update NIF
            if len_3utr >= tol_long3utr:
                self.nif_long3utr = True
            else:
                self.nif_long3utr = False
        else:
            seq_3utr = ""
            len_3utr = np.nan

        # Compute
        self.gc_content = pd.Series( [GC(seq_5utr), GC(self.sequences["cds"]), GC(seq_3utr), [*map(GC, self.sequences["exons"])], [*map(GC, self.sequences["introns"])]], index = ["5'utr", "cds", "3'utr", "exons", "introns"], name="gc")
        self.length = pd.Series( [len_5utr, len(self.sequences["cds"]), len_3utr, [*map(len, self.sequences["exons"])], [*map(len, self.sequences["introns"])]], index = ["5'utr", "cds", "3'utr", "exons", "introns"], name="length")

        # Concatenate
        # GC
        Se_gc = self.gc_content.copy()
        Se_gc.index = Se_gc.index.map(lambda x:"GC|__|%s"%x)
        # Length
        Se_length = self.length.copy()
        Se_length.index = Se_length.index.map(lambda x:"LENGTH|__|%s"%x)

        return pd.Series(pd.concat([Se_gc, Se_length]), name=self.id)

    # Organize NIFs
    def nifs(self, tol_long3utr, tol_window):
        """
        Organize NIFs into a boolean format
        """
        return pd.Series([self.nif_dej, self.nif_uorf, self.nif_long3utr, self.nif_termination_gc], index=["NIF|__|dEJ", "NIF|__|uORF", "NIF|__|long-3'utr", "NIF|__|GC(cds[-%d:termination_codon])"%(tol_window)], name="nif")

    # Compile all data
    def compute(self, kozak_type, tol_min_codons, tol_window, tol_dej, tol_long3utr, file_errors):
        """
        Compute the entire pipeline
        """
        # Upstream Open Reading Frames
        Se_uorfs = self.uORFs(kozak_type=kozak_type, tol_min_codons=tol_min_codons, file_errors=file_errors)
        # Upstream GC of Termination Codon
        query_termination_gc = self.termination_gc(tol_window=tol_window)
        # Downstream Exon Junctions
        query_dej = self.dEJs(tol_dej=tol_dej)
        # Number of Exons
        num_exons = len(self.sequences["exons"])
        # Sequence Stats
        Se_stats = self.sequence_stats(tol_long3utr=tol_long3utr)
        # Misc
        Se_misc = pd.Series([num_exons, query_dej], index=["EXONS|__|count", "dEJs|__|penultimate-exon>=%d(nts)"%(tol_dej)])
        # NIFs
        Se_nifs = self.nifs(tol_long3utr=tol_long3utr, tol_window=tol_window)

        # Output
        return pd.Series(pd.concat([Se_nifs, Se_stats, Se_misc, Se_uorfs]), name=self.id)



def main():

    parser= argparse.ArgumentParser(description=__doc__)
    # I/O
    parser.add_argument("-p","--path_premrna",help="pre-mRNA .fasta filepath")
    parser.add_argument("-c","--path_cds",help="CDS .fasta filepath")
    parser.add_argument("-o", "--out_dir", type=str, default = os.getcwd(), help = "Path/to/existing-directory-for-output.  Warning: Do not use `~` to specify /home/[user]/. [Default: Current directory]" )
    parser.add_argument("-n", "--project_name",type=str, default = None,help="Project name for files (no spaces). [Default: Date]")
    # uORF Start Sites
    parser.add_argument("-k","--kozak_type",type=int,default=1,help="Stringency of kozak sequences ranging from 0 - 3. 0: `ATG`; 1: `[A|G]..ATG or ATG[G]; 2: `[A|G]..ATG[G]; 3: `[A|G].[C]ATG[G]. [Default: 1]")
    # NIFs
    parser.add_argument("--tol_min_codons",type=int,default=10,help="Mimumum codons considered in output (not including STOP) [Default: 10 codons]")
    parser.add_argument("--tol_dej",type=int,default=55,help="Minimum distance between dEJ and termination codon [Default: 55 nts]")
    parser.add_argument("--tol_window", type=int, default=10, help="Number of nucleotides at the end of CDS to check GC-Content")
    parser.add_argument("--tol_long3utr", type=int, default=1200, help="Minimum number of nucleotides to be considered a long 3'utr")

    # parser.add_argument("--metadata", default=True, help="[True | False] to contain metadata in output.  Requirements: (1) Connected to internet; (2) `bioservcies` module installed. [Default: True]")
    # parser.add_argument("--n_jobs", type=int, default=2, help="Number of threads to use for mapping metadata from `ensembl`. `-1` for maximum number of available cores. [Default: 2]")

    opts = parser.parse_args()

    # Start Time
    start_time = time.time()

    # Init vars
    if opts.project_name == None: opts.project_name = os.path.split(opts.path_premrna)[1].split(".")[0]
    if opts.out_dir.endswith("/"): opts.out_dir = opts.out_dir[::-1]

    # Check if directory exists, make one if not
    print("Checking if directory exists, make one if not...", sep="\t", file=sys.stderr)
    if not os.path.exists(opts.out_dir):
        os.makedirs(opts.out_dir)

    # NIF | Finder
    print("================================\n NIF Finder\n================================", file=sys.stdout)
    print(opts.project_name)
    print("pre-mRNA:", opts.path_premrna, sep="\t", file=sys.stdout)
    print("CDS:", opts.path_cds, sep="\t", file=sys.stdout)
    print("Path:", opts.out_dir, sep="\t", file=sys.stdout)
    print("================\n Parameters \n================", file=sys.stdout)
    print("kozak_type: ", opts.kozak_type, sep="\t", file=sys.stdout)
    print("tol_min_codons: ", opts.tol_min_codons, sep="\t", file=sys.stdout)
    print("tol_dej: ", opts.tol_dej, sep="\t", file=sys.stdout)
    print("tol_long3utr", opts.tol_long3utr, sep="\t", file=sys.stdout)
    print("tol_window", opts.tol_window, sep="\t", file=sys.stdout)
    print("================\n Summary \n================", file=sys.stdout)


    # Error Log
    file_filter = open(opts.out_dir + "/" + opts.project_name + "__filter.log", "w")
    print("transcript_id", "drop_status", "description", sep="\t", file=file_filter)
    # Console Log
    file_errors = open(opts.out_dir + "/" + opts.project_name + "__errors.log", "w")


    # Vessel for mRNA objects
    # =======================
    D_id_mRNA = OrderedDict()

    # Load pre-mRNA Sequences
    # =======================
    handle_premrna = {True:gzip.open(opts.path_premrna, "rt"), False:open(opts.path_premrna, "r")}[opts.path_premrna.endswith(".gz")]
    if tqdm_available:
        iter_premrna = tqdm(SeqIO.parse(handle_premrna, "fasta"), desc="Reading pre-mRNA. .. ... ..... .......")
    else:
        print("Reading pre-mRNA. .. ... ..... .......", file=sys.stderr)
        iter_premrna = SeqIO.parse(handle_premrna, "fasta")
    for record in iter_premrna:
        # Synthesize mRNA object
        mrna = mRNA(id=record.id)
        mrna.sequences["pre-mRNA"] = str(record.seq)
        D_id_mRNA[record.id] = mrna
    set_premrna = set(D_id_mRNA.keys())

    # Load CDS Sequences
    # =======================
    handle_cds = {True:gzip.open(opts.path_cds, "rt"), False:open(opts.path_cds, "r")}[opts.path_cds.endswith(".gz")]

    if tqdm_available:
        iter_cds = tqdm(SeqIO.parse(handle_cds, "fasta"), desc="Analyzing CDS. .. ... ..... .......")
    else:
        print("Analyzing CDS. .. ... ..... .......", file=sys.stderr)
        iter_cds = SeqIO.parse(handle_cds, "fasta")
    for record in iter_cds:
        # Update mRNA object
        if record.id in set_premrna:
            mrna = D_id_mRNA[record.id]
            mrna.sequences["cds"] = str(record.seq).upper()
            # In-silico splicing
            mrna.splice()
        # Synthesize mRNA object and flag
        else:
            D_id_mRNA[record.id] = mRNA(record.id)
            mrna = D_id_mRNA[record.id]
            mrna.errors.append("Fatal: Missing `pre-mRNA`")
            mrna.fatal_error = True

    # Error Checking
    # =======================
    if tqdm_available:
        iter_qc= tqdm(D_id_mRNA.items(), desc="Quality control. .. ... ..... ........")
    else:
        print("Quality control. .. ... ..... .......", file=sys.stderr)
        iter_qc = D_id_mRNA.items()

    for id, mrna  in iter_qc:
        # Check if CDS Exists
        if "cds" not in mrna.sequences:
            mrna.errors.append("Fatal: Missing `CDS`")
            mrna.fatal_error = True
        # Check for missing 5'UTR, 3'UTR, or if CDS is not in mRNA
        if len(mrna.errors) > 0:
            print(id, mrna.fatal_error, "; ".join(mrna.errors), sep="\t", file=file_filter)

    # Quality Control
    # =======================
    D_id_mRNA_qa = OrderedDict([*filter(lambda item:item[1].fatal_error == False, D_id_mRNA.items())])

    # Compute NIFs
    # =======================
    args = { "kozak_type":opts.kozak_type, "tol_min_codons":opts.tol_min_codons, "tol_window":opts.tol_window, "tol_dej":opts.tol_dej, "tol_long3utr":opts.tol_long3utr, "file_errors":file_errors}

    if tqdm_available:
        iter_nifs = tqdm(D_id_mRNA_qa.items(), desc="Computing NIFs. .. ... ..... ........")
    else:
        print("Computing NIFs. .. ... ..... ........", file=sys.stderr)
        iter_nifs = D_id_mRNA_qa.items()

    data = list()
    for id, mrna in iter_nifs:
        data.append(mrna.compute(**args))

    # Future: Parallel?
    # def compute_attributes(mrna, args):
    #     return mrna.compute(**args)
    # data = Parallel(n_jobs=n_jobs)(delayed(compute_attributes)(mrna=mrna, args=args) for id, mrna in D_id_mRNA_qa.items())
    print("Writing files. .. ... ..... ........", file=sys.stderr)
    df_attrs = pd.DataFrame(data)
    df_attrs["EXONS|__|count"] = df_attrs["EXONS|__|count"].astype(int)
    df_attrs.to_csv(opts.out_dir + "/" + opts.project_name + "__nifs.tsv.gz", sep="\t", compression="gzip")

    # Write uORF fasta files
    # =======================
    with open("%s/%s__uorfs.fna"%(opts.out_dir, opts.project_name), "w") as file_fna:
        for id, mrna in D_id_mRNA_qa.items():
            for id_uorf, seq_uorf in mrna.uorf_fna.items():
                print(">%s\n%s"%(id_uorf, seq_uorf), file=file_fna)
    with open("%s/%s__uorfs.faa"%(opts.out_dir, opts.project_name), "w") as file_faa:
        for id, mrna in D_id_mRNA_qa.items():
            for id_uorf, seq_uorf in mrna.uorf_faa.items():
                print(">%s\n%s"%(id_uorf, seq_uorf), file=file_faa)

    file_errors.close()
    file_filter.close()

    # Summary
    print("Transcripts Processed:", df_attrs.shape[0], sep="\t", file=sys.stdout)
    print("Transcripts Dropped:", len(set_premrna) - df_attrs.shape[0], sep="\t", file=sys.stdout)

    print("uORFs:", df_attrs["NIF|__|uORF"].sum(), sep="\t", file=sys.stdout)
    print("dEJs:",  df_attrs["NIF|__|dEJ"].sum(), sep="\t", file=sys.stdout)
    print("Long 3'UTRs:",  df_attrs["NIF|__|long-3'utr"].sum(), sep="\t", file=sys.stdout)

    print("Runtime: " + str(time.time() - start_time) + " seconds", file=sys.stderr)

if __name__ == "__main__":
    main()
