import os
import sys
import argparse
from statistics import mean
'''
INPUT: tab-delimited blastn output. Assuming that taxonomy ID is in this format:
-outfmt "6 qseqid sseqid staxids sscinames scomnames sskingdoms pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp"

This script also assumes that input has been filtered by 90% identity.
$ awk '{if ($7 > 90) print}' all_results.tsv > all_results.90perc.tsv

And very vague hits have been removed, too:
$ grep -v -e uncultured -e Uncultured -e chloroplast -e Unidentified -e unidentified all_results.90perc.tsv > all_results.90perc.noUnculturedUnidentifiedChloroplast.tsv

This script also assumes that taxonkit is installed and that the data has been downloaded to the default directory.
mamba create -n taxonkit taxonkit -y
conda activate taxonkit
mkdir ~/.taxonkit
cd ~/.taxonkit
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar xzvf taxdump.tar.gz

Usage: python computeLCA.py filtered_blast.tsv > Table_of_TaxonIDs.tsv

Currently, the output looks like this:
ASV_0   6833    cellular organisms;Eukaryota;Opisthokonta;Metazoa;Eumetazoa;Bilateria;Protostomia;Ecdysozoa;Panarthropoda;Arthropoda;Mandibulata;Pancrustacea;Crustacea;Multicrustacea;Hexanauplia;Copepoda;Neocopepoda;Gymnoplea;Calanoida Calanoida
ASV_1   6833    cellular organisms;Eukaryota;Opisthokonta;Metazoa;Eumetazoa;Bilateria;Protostomia;Ecdysozoa;Panarthropoda;Arthropoda;Mandibulata;Pancrustacea;Crustacea;Multicrustacea;Hexanauplia;Copepoda;Neocopepoda;Gymnoplea;Calanoida Calanoida
ASV_2   6833    cellular organisms;Eukaryota;Opisthokonta;Metazoa;Eumetazoa;Bilateria;Protostomia;Ecdysozoa;Panarthropoda;Arthropoda;Mandibulata;Pancrustacea;Crustacea;Multicrustacea;Hexanauplia;Copepoda;Neocopepoda;Gymnoplea;Calanoida Calanoida
ASV_3   6833    cellular organisms;Eukaryota;Opisthokonta;Metazoa;Eumetazoa;Bilateria;Protostomia;Ecdysozoa;Panarthropoda;Arthropoda;Mandibulata;Pancrustacea;Crustacea;Multicrustacea;Hexanauplia;Copepoda;Neocopepoda;Gymnoplea;Calanoida Calanoida
ASV_4   349674  cellular organisms;Eukaryota;Opisthokonta;Metazoa;Eumetazoa;Bilateria;Protostomia;Ecdysozoa;Panarthropoda;Arthropoda;Mandibulata;Pancrustacea;Crustacea;Multicrustacea;Hexanauplia;Copepoda;Neocopepoda;Podoplea;Poecilostomatoida;Oncaeidae;Oncaea  Oncaea
ASV_5   349674  cellular organisms;Eukaryota;Opisthokonta;Metazoa;Eumetazoa;Bilateria;Protostomia;Ecdysozoa;Panarthropoda;Arthropoda;Mandibulata;Pancrustacea;Crustacea;Multicrustacea;Hexanauplia;Copepoda;Neocopepoda;Podoplea;Poecilostomatoida;Oncaeidae;Oncaea  Oncaea
....
'''

def main():
    DATA_DIR = '~/.taxonkit/' # Where taxonkit data is stored.

    parser = argparse.ArgumentParser(description='Parses blastn output, computes LCA for all queries using taxonkit.')
    parser.add_argument('blastn', metavar='FILE', type=str,
                        help='''Filename of blastn outfmt 6 output. Ideally
                        "6 qseqid sseqid staxids sscinames scomnames sskingdoms pident length
                        qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp"''')
    parser.add_argument('--data-dir', metavar='DIR', type=str, default=DATA_DIR,
                        help=f'Where taxonkit data is stored. Defaults to "{DATA_DIR}".')
    args = parser.parse_args()

    data_dir = args.data_dir

    current_query = None
    taxons = None

    list_of_queries = []
    list_of_percentages = []
    list_of_lengths = []
    with open(args.blastn) as fh, open('TEMP_TAXONS.txt', 'w') as out:

        for line in fh:
            ll = line.rstrip().split('\t')
            query = ll[0]
            taxid = ll[2]
            perc = float(ll[6])
            length = float(ll[7])
            if ';' in taxid:
                # in some cases, an NCBI-NT entry has several taxonomy IDs like 1;2;321;5123413.
                # It seems that these few entries have the entire taxonomy tree. Just use the last one.
                taxid = taxid.split(';')[-1]

            #if ll[3] == 'N/A':
                # this Taxonomy ID comes from Entrez but it's not in taxonkit. not sure how this happens
                #print(ll)
                #continue

            if current_query != query:
                if taxons:
                    # write out the taxons to temporary files, removing duplicates
                    out.write(current_query + '\t' + ' '.join(set(taxons)) + '\n')
                    list_of_percentages.append(mean(percentages))
                    list_of_lengths.append(mean(lengths))
                # printed all for the last query, now start collecting for the new query
                taxons = [taxid]
                percentages = [perc]
                lengths = [length]
                current_query = query
                list_of_queries.append(current_query)
            else:
                # current query is same as before, add to list
                taxons.append(taxid)
                percentages.append(perc)
                lengths.append(length)

        # print for the very last query, too
        if taxons:
            #out.write(' '.join(set(taxons)) + '\n')
            out.write(current_query + '\t' + ' '.join(set(taxons)) + '\n')
            list_of_percentages.append(mean(percentages))
            list_of_lengths.append(mean(lengths))

    all_taxa = os.popen('taxonkit lca -i 2  TEMP_TAXONS.txt | cut -f 3 | taxonkit lineage -n').readlines()
    #all_taxa = os.popen('taxonkit lca -D TEMP_TAXONS.txt | cut -f 2 | taxonkit lineage -n | taxonkit reformat --format "{k};{p};{c};{o};{f};{g};{s}"').read().split('\n')

    #all_taxa = list(filter(None, all_taxa)) # the last line is empty for some reason. Extra linebreak?
    assert len(all_taxa) == len(list_of_queries), f'ERROR: there are {len(list_of_queries)} queries but {len(all_taxa)} taxa in the TEMP_TAXONS.txt file'

    reformatted_taxa = os.popen('taxonkit lca -i 2 TEMP_TAXONS.txt | cut -f 3 | taxonkit lineage -n | taxonkit reformat --format "{o};{f};{g};{s}"').read().split('\n')

    for q, t, l, p, x  in zip(list_of_queries, all_taxa, reformatted_taxa, list_of_percentages, list_of_lengths):
        # 119195  cellular organisms;Eukaryota;Opisthokonta;Metazoa;Eumetazoa;Bilateria;Deuterostomia;Chordata;Craniata;Vertebrata;Gnathostomata;Chondrichthyes;Elasmobranchii;Selachii;Galeomorphii;Galeoidea                                        Galeoidea        Eukaryota;Chordata;Chondrichthyes;;;;
        newnames = t.split('\t')
        # ['55137', 'cellular organisms;Eukaryota;Opisthokonta;Metazoa;Eumetazoa;Bilateria;Deuterostomia;Chordata;Craniata;Vertebrata;Gnathostomata;Chondrichthyes;Elasmobranchii;Batoidea;Myliobatiformes;Myliobatidae;Myliobatinae;Myliobatis', 'Myliobatis']
        lineage = newnames[1].split(';')
        if 'Chondrichthyes' in lineage:
            lineage = ';'.join(lineage[:16])
        else:
            lineage = ';'.join(lineage[:17])
        #186625  cellular organisms;Eukaryota;Opisthokonta;Metazoa;Eumetazoa;Bilateria;Deuterostomia;Chordata;Craniata;Vertebrata;Gnathostomata;Teleostomi;Euteleostomi;Actinopterygii;Actinopteri;Neopterygii;Teleostei;Osteoglossocephalai;Clupeocephala       Clupeocephala   ;;
        family_genus_s = l.rstrip().split('\t')[-1]
        newlineage = lineage+';'+family_genus_s
        last_one = ''
        for i in newlineage.split(';'):
            if i:
                last_one = i

        print(q + '\t' + lineage+';'+family_genus_s + '\t' + last_one + '\t%.2f\t%.0f'%(p, x) )
        #print(f'{q}\t{newlist}\t{lastone}')

    os.remove("TEMP_TAXONS.txt")

main()
