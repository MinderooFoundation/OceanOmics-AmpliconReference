# Minderoo Foundation OceanOmics amplicon databases

This repository contains the Minderoo Foundatin OceanOmics amplicon databases and the code to generate them.

Currently the repository is structured into three folders: one for 12S genes, one for 16S genes, and one for publicly available mitogenomes.

The final database formatted for BLAST with taxonomy IDs is in databases/12S.16S.Mitogenomes.fasta

# Scripts

- getAssemblies.py - for each 12S and 16S. Writes an entrez script that downloads 16S and 12S sequences for Australian sharks, bony fish, marine mammals and reptiles. The equivalent in the mitogenomes folder downloads entire genomes.
- doAllQC.py - self-blasts the genes and writes out potential mislabels and contaminants to build a final BLAST-formatted database with taxonomy IDs.

- mergeDatabases.py - builds a final BLAST database of the three folders. Removes duplicates like voucher sequences that contain a 12S and a 16S gene.


All larger files have been gzipped.
