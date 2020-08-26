
Introduction

omicR creates fasta files, downloads genomes from NCBI using the refseq number, creates databases to run BLAST+, runs BLAST+ and filters these results to obtain the best match per sequence. 

These scripts can be used to run BLAST alignment of short-read (DArTseq data) and long-read sequences (Illumina, PacBio… etc). You can use reference genomes from NCBI, genomes from your private collection, contigs, scaffolds or any other genetic sequence that you would like to use as reference. 

This version is suitable for MAC, Windows and Linux users. 

Watch the tutorial video in YouTube omicR with Graphical User Interface (~20 min) https://youtu.be/pdMio2vj-FM 

You need to install the following:

a.	Python V3 or latest: https://www.python.org/downloads/

b.	Biopython: https://biopython.org/

c.	omicR: https://github.com/BTalamantesBecerra/omicR Download the following files from GitHub: “omicR.py” and “Currito.ico” and save them in your python working directory. 

d. BLAST+ latest version: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ 


-Open the script “omicR.py” in Python and run it.

-This will open a window where you can run the scripts. 

As general practice avoid installing your software in directories such as “C://Program files/” as the space between words will cause problems. 
Remember to add all the software you install into your System Variables Path. 
