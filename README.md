```
==================================================================================================
NMD | Attribute Finder
==================================================================================================

# 1: Install Python
	(1.0) Download an install Python 3 from https://www.continuum.io/downloads

# 2: Install Dependencies [`numpy`, `pandas`, `joblib`, `bioservices`, and `biopython`]
	(2.0) Open terminal
	(2.1) Navigate to directory of `install_modules.sh` [e.g. `cd ~/Downloads/nmd_attribute-finder/``]
	(2.2) Change permissions on `install_modules.sh` [e.g. `chmod 775 install_modules.sh`]
	(2.3) Run shell script to install modules (alternatively run commands individually) [e.g. `./install_modules.sh`]

# 3: Run program
	(3.0) `python nmd_attribute-finder.py [commands]`
	(3.1) Example command: `python nif_finder.py -p /path/to/pre-mrna.fasta -c /path/to/cds.fasta -o /path/to/output-directory`

$ python nif_finder.py --help
usage: nif_finder.py [-h] [-p PATH_PREMRNA] [-c PATH_CDS] [-o OUT_DIR]
                     [-n PROJECT_NAME] [-k KOZAK_TYPE]
                     [--tol_min_codons TOL_MIN_CODONS] [--tol_dej TOL_DEJ]
                     [--tol_window TOL_WINDOW] [--tol_long3utr TOL_LONG3UTR]

optional arguments:
  -h, --help            show this help message and exit
  -p PATH_PREMRNA, --path_premrna PATH_PREMRNA
                        pre-mRNA .fasta filepath
  -c PATH_CDS, --path_cds PATH_CDS
                        CDS .fasta filepath
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