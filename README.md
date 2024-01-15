# Minderoo Foundation OceanOmics amplicon databases

This repository contains the Minderoo Foundatin OceanOmics amplicon databases and the code to generate them.

Currently the repository is structured into three folders: one for 12S genes, one for 16S genes, one for COI genes, and one for publicly available mitogenomes.

The final database formatted for BLAST with taxonomy IDs is in databases/12S.v0.10.16S.v0.4.Mitogenomes.v0.1.fasta.gz

# Scripts

- getAssemblies.py - for each 12S, 16S, and COI. Writes an entrez script that downloads 12S, 16S, and COI sequences for Australian sharks, bony fish, marine mammals and reptiles. The equivalent in the mitogenomes folder downloads entire genomes.
- doAllQC.py - self-blasts the genes and writes out potential mislabels and contaminants to build a final BLAST-formatted database with taxonomy IDs.

- mergeDatabases.py - builds a final BLAST database of the three folders. Removes duplicates like voucher sequences that contain a 12S, a 16S, or a COI gene.

All larger files have been gzipped.
