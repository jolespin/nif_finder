#!/usr/bin/env python
import sys, os, argparse, gzip
from collections import OrderedDict
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser
from tqdm import tqdm

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2022.01.05"

complement = str.maketrans('ACGTacgt','TGCAtgca')
def reverse_complement(sequence):
    """
    Inputs a sequence of DNA and returns the reverse complement se$
    Outputs a reversed string of the input where A>T, T>A, G>C, an$
    """

    reverse_complement = sequence.translate(complement)[::-1]
    return reverse_complement



def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -a <assembly.fa[.gz]> -g <gene_models.gtf[.gz]>  -o <output.fa>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-a","--assembly", type=str, help = "path/to/assembly.fa[.gz]")
    parser.add_argument("-g","--gtf", type=str, help = "path/to/gene_models.gtf[.gz]")
    parser.add_argument("-o","--output", type=str, default="stdout", help = "path/to/pre-mRNA.fa")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename
    

    if opts.output == "stdout":
        f_out = sys.stdout 
    else:
        f_out = open(opts.output, "w")



    # Assembly
    reference_to_sequence = OrderedDict()
    with {True:gzip.open, False:open}[opts.assembly.endswith(".gz")](opts.assembly, "rt") as f:
        for id_record, sequence in tqdm(SimpleFastaParser(f), desc="Reading: {}".format(opts.assembly), unit=" sequence"):
            id_record = id_record.split(" ")[0]
            reference_to_sequence[id_record] = sequence
    reference_to_sequence = pd.Series(reference_to_sequence)

    # GTF
    exon_info = list()
    with {True:gzip.open, False:open}[opts.gtf.endswith(".gz")](opts.gtf, "rt") as f:
        for line in tqdm(f.readlines(), desc="Reading: {}".format(opts.gtf), unit=" line"):
            line = line.strip()
            if "\texon\t" in line:
                fields = line.split("\t")
                id_reference = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                sense = fields[6]
                id_gene = fields[-1].split('gene_id "')[1].split('"')[0] # Need a regex instead
                id_transcript = fields[-1].split('transcript_id "')[1].split('"')[0] # Need a regex instead
                exon_number = fields[-1].split('exon_number "')[1].split('"')[0] # Need a regex instead
                transcript_version = fields[-1].split('transcript_version "')[1].split('"')[0] # Need a regex instead
                exon_info.append([id_reference, start, end, sense, exon_number, id_gene, "{}.{}".format(id_transcript,transcript_version)])
            
    df_exons = pd.DataFrame(exon_info, columns=["reference_id", "start", "end", "sense", "exon_number", "gene_id", "transcript_id"]).sort_values(["reference_id","start"]).reset_index(drop=True)



    # Get preMRNA
    transcript_to_premRNA = dict()

    error_transcripts = list()
    for id_transcript, df in tqdm(df_exons.groupby("transcript_id", axis=0), desc="Merging exon sequences", unit=" transcript"):
        values = df.values
        number_of_exons = values.shape[0]
        id_reference = str(values[0][0])
        sense = values[0][3]
        
        exon_sequences = list()
        if id_reference in reference_to_sequence:
            sequence_reference = reference_to_sequence[id_reference]
            
            # More than 1 exon
            if number_of_exons > 1:
                for i in range(0, number_of_exons):
                    row = values[i]
                    exon_start = row[1] - 1 
                    exon_end = row[2]
                    sequence_exon = sequence_reference[exon_start:exon_end].upper()
                    exon_sequences.append(sequence_exon)
                    
                    if i < (number_of_exons - 1):
                        intron_start = exon_end
                        intron_end = values[i+1][1] - 1
                        sequence_intron = sequence_reference[intron_start:intron_end].lower()
                        exon_sequences.append(sequence_intron)
                        
            # Only 1 exon
            else:
                exon_start = values[0][1] - 1
                exon_end = values[0][2]
                sequence_exon = sequence_reference[exon_start:exon_end]
                exon_sequences.append(sequence_exon)
                
            sequence_premRNA = "".join(exon_sequences)
            
            if sense == "-":
                sequence_premRNA = reverse_complement(sequence_premRNA)
                
            transcript_to_premRNA[id_transcript] = sequence_premRNA
        else:
            error_transcripts.append(values)
    transcript_to_premRNA = pd.Series(transcript_to_premRNA)

    # Write pre-mRNA
    for id, seq in tqdm(transcript_to_premRNA.items(), desc="Writing pre-mRNA sequences", unit=" sequence", total=len(transcript_to_premRNA)):
        print(">{}\n{}".format(id, seq), file=f_out, end="\n", sep="")
    if f_out != sys.stdout: 
        f_out.close()




    
if __name__ == "__main__":
    main()
    
                