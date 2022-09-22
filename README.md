### NMD Inducible Factor Finder

### Install: 
Create a conda environment with all dependencies: 
`conda create -n nif-finder_env -c conda-forge pandas numpy biopython joblib tqdm bioservices -y`

### Usage: 
```


$ python nif_finder.py --help
usage: nif_finder.py [-h] [-p PATH_PREMRNA] [-c PATH_CDS] [-o OUT_DIR]
                     [-n PROJECT_NAME] [-k KOZAK_TYPE]
                     [--tol_min_codons TOL_MIN_CODONS] [--tol_dej TOL_DEJ]
                     [--tol_window TOL_WINDOW] [--tol_long3utr TOL_LONG3UTR]

optional arguments:
  -h, --help            show this help message and exit
  -p PATH_PREMRNA, --path_premrna PATH_PREMRNA
                        pre-mRNA .fasta filepath where exons are upper case and introns are lower case
  -c PATH_CDS, --path_cds PATH_CDS
                        CDS .fasta filepath of upper case CDS sequences
  -o OUT_DIR, --out_dir OUT_DIR
                        Path/to/existing-directory-for-output. Warning: Do not
                        use `~` to specify /home/[user]/. [Default: Current
                        directory]
  -n PROJECT_NAME, --project_name PROJECT_NAME
                        Project name for files (no spaces). [Default: Date]
  -k KOZAK_TYPE, --kozak_type KOZAK_TYPE
                        Stringency of kozak sequences ranging from 0 - 3. 0:
                        `ATG`; 1: `[A|G]..ATG or ATG[G]; 2: `[A|G]..ATG[G]; 3:
                        `[A|G].[C]ATG[G]. [Default: 1]
  --tol_min_codons TOL_MIN_CODONS
                        Mimumum codons considered in output (not including
                        STOP) [Default: 10 codons]
  --tol_dej TOL_DEJ     Minimum distance between dEJ and termination codon
                        [Default: 55 nts]
  --tol_window TOL_WINDOW
                        Number of nucleotides at the end of CDS to check GC-
                        Content
  --tol_long3utr TOL_LONG3UTR
                        Minimum number of nucleotides to be considered a long
                        3'utr
```

### Input: 
* pre-mRNA sequences: Use `get_premrna_sequences_from_gtf.py` to generate a pre-mRNA fasta file with exons as upper case and introns as lower case.

```
python get_premrna_sequences_from_gtf.py -h
usage: get_premrna_sequences_from_gtf.py -a <assembly.fa[.gz]> -g <gene_models.gtf[.gz]>  -o <output.fa>

    Running: get_premrna_sequences_from_gtf.py v2022.01.05 via Python v3.8.5 

optional arguments:
  -h, --help            show this help message and exit
  -a ASSEMBLY, --assembly ASSEMBLY
                        path/to/assembly.fa[.gz]
  -g GTF, --gtf GTF     path/to/gene_models.gtf[.gz]
  -o OUTPUT, --output OUTPUT
                        path/to/pre-mRNA.fa
                        
# Usage: 
python get_premrna_sequences_from_gtf.py -a Mus_musculus.GRCm39.dna.primary_assembly.fa.gz -g Mus_musculus.GRCm39.105.gtf.gz  > Mus_musculus.GRCm39.105.gtf.premrna.fa
  
Reading: /Users/jespinoz/Documents/NMD/Mus_musculus.GRCm39/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz: 61 sequence [00:33,  1.84 sequence/s]
Reading: /Users/jespinoz/Documents/NMD/Mus_musculus.GRCm39/Mus_musculus.GRCm39.105.gtf.gz: 100%|███████████████████████████████████████████████████████████████████████████████████████████████████████| 1866443/1866443 [00:08<00:00, 230844.77 line/s]
Merging exon sequences: 100%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 142435/142435 [00:22<00:00, 6230.19 transcript/s]
Writing pre-mRNA sequences: 100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 142435/142435 [00:06<00:00, 21961.93 sequence/s]                        
```
* CDS sequences: Download this from Ensemble and make sure it's the same build as the GTF used for pre-mRNA sequences.

### License: 
BSD-3

### Cite:
* Espinoza, Josh L. (2022): NIF Finder. figshare. Software. https://doi.org/10.6084/m9.figshare.17775977.v1 

* Shum EY, Espinoza JL, Ramaiah M, Wilkinson MF. Identification of novel post-transcriptional features in olfactory receptor family mRNAs. Nucleic Acids Res. 2015 Oct 30;43(19):9314-26. doi: 10.1093/nar/gkv324. PubMed PMID: 25908788; PubMed Central PMCID: PMC4627058.

* Lou CH, Shao A, Shum EY, Espinoza JL, Huang L, Karam R, Wilkinson MF. Posttranscriptional control of the stem cell and neurogenic programs by the nonsense-mediated RNA decay pathway. Cell Rep. 2014 Feb 27;6(4):748-64. doi: 10.1016/j.celrep.2014.01.028. PubMed PMID: 24529710; PubMed Central PMCID: PMC3962089.

* Domingo D, Nawaz U, Corbett M, Espinoza JL, Tatton-Brown K, Coman D, Wilkinson MF, Gecz J, Jolly LA. UPF3B mutations including a novel synonymous variant associated with absent speech implicate nonsense mediated mRNA decay as a regulator of neurodevelopmental disorder gene networks. Hum Mol Genet. 2020 Jul 15:ddaa151. doi: 10.1093/hmg/ddaa151. PMID: 32667670

* Huang L, Shum EY, Jones SH, Lou CH, Dumdie J, Kim H, Roberts AJ, Jolly LA, Espinoza JL, Skarbrevik DM, Phan MH, Cook-Andersen H, Swerdlow NR, Gecz J, Wilkinson MF. A Upf3b-mutant mouse model with behavioral and neurogenesis defects. Mol Psychiatry. 2018 Aug;23(8):1773-1786. doi: 10.1038/mp.2017.173. PubMed PMID: 28948974; PubMed Central PMCID: PMC5869067.

* Shum EY, Jones SH, Shao A, Dumdie J, Krause MD, Chan WK, Lou CH, Espinoza JL, Song HW, Phan MH, Ramaiah M, Huang L, McCarrey JR, Peterson KJ, De Rooij DG, Cook-Andersen H, Wilkinson MF. The Antagonistic Gene Paralogs Upf3a and Upf3b Govern Nonsense-Mediated RNA Decay. Cell. 2016 Apr 7;165(2):382-95. doi: 10.1016/j.cell.2016.02.046. PubMed PMID: 27040500; PubMed Central PMCID: PMC4826573.




