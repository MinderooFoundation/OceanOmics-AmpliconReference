import os
from Bio import SeqIO
from collections import defaultdict, Counter
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('-i', '--infile', help='Name of input fasta downloaded',
        required = True)
args = parser.parse_args()

if __name__ == '__main__':
    all_seqs = []
    removed_classes = defaultdict(int)
    removed = {} # key: sequence id. value: reason why removed.
    
    if not os.path.exists('0-taxoncheck'):
        os.mkdir('0-taxoncheck')
    with open('0-taxoncheck/TAXONKIT_names.txt', 'w') as out:
        for s in SeqIO.parse(args.infile, 'fasta'):
            # remove vague species
            if 'cf.' in s.description or 'sp.' in s.description or ' x ' in s.description or ' X ' in s.description \
                    or 'aff.' in s.description:
                removed[s.id] =  f'vague species identification ({s.description})'
                continue
            # get the taxonomy ID
            # 16S_NC_015120.1:1086-2750 Merluccius merluccius mitochondrion, complete genome
            if 'PREDICTED' in s.description:
                species  = s.description.split(' ')[2:4]
            else:
                species  = s.description.split(' ')[1:3]
            if len(species) != 2:
                # 16S_gi|NC_057650.1|ref|NC_057650.1|tax|531329|Plectorhinchus chaetodonoides
                species = s.description.split('|')[-1]
            else:
                species = ' '.join(species)
            # some species like 'thynnus thynnus', need to uppercas
            species = species.capitalize()
            out.write(f'{species}\n')
            all_seqs.append(s)

    # get the newest taxonomy IDs.
    os.popen(f'taxonkit name2taxid {out.name} | taxonkit lineage -i 2  > 0-taxoncheck/TAXONKIT_names_taxids.txt').read()

    taxid_dict = {}
    with open('0-taxoncheck/TAXONKIT_names_taxids.txt') as fh:
        for line in fh:
            ll = line.rstrip().split('\t')
            try:
                name, taxid = ll[0], ll[1]
            except:
                continue
            taxid_dict[name] = taxid

    # we will have names without taxonomy IDs - remove those sequences too.
    # now self-blast - we need to make the blast database
    if not os.path.exists('1-selfblast'):
        os.mkdir('1-selfblast')
    with open('1-selfblast/selfblastdb.fasta', 'w') as fastaout\
            , open('1-selfblast/selfblastdb_taxIDs.txt','w') as taxout:
        for s in all_seqs:
            if s.id in removed:
                # skip already removed seqs
                continue

            # naming scheme:
            # ID, taxID, speciesname, OLD description
            species  = s.description.split(' ')[1:3]
            if 'PREDICTED' in s.description:
                species = s.description.split(' ')[2:4]

            if len(species) != 2:
                # 16S_gi|NC_057650.1|ref|NC_057650.1|tax|531329|Plectorhinchus chaetodonoides
                species = s.description.split('|')[-1]
            else:
                species = ' '.join(species)
            if species not in taxid_dict:
                removed[s.id] = f'No taxonomy ID ({species})'
                continue
            new_id = s.id
            if '|' in new_id:
                # 16S_gi|NC_057650.1|ref|NC_057650.1|tax|531329| 
                # becomes
                # NC_057650.1
                new_id = new_id.split('|')[3]
            if 'whole genome' in s.description.lower() or 'linkage group' in s.description.lower():
                removed[s.id] = f'Looks like nuclear genome ({len(s)} {s.description})'
                continue
            new_name = f'{new_id} {taxid_dict[species]} {species} {s.description}'
            fastaout.write(f'>{new_name}\n{str(s.seq)}\n')
            taxout.write(f'{new_id} {taxid_dict[species]}\n')

    # make the blast database
    os.popen('makeblastdb -dbtype nucl -in 1-selfblast/selfblastdb.fasta -parse_seqids -taxid_map 1-selfblast/selfblastdb_taxIDs.txt').read()
    if not os.path.exists('1-selfblast/taxdb.tar.gz'):
        os.popen('wget wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz -O 1-selfblast/taxdb.tar.gz').read()
        os.popen('tar xzvf 1-selfblast/taxdb.tar.gz --directory 1-selfblast/').read()

    # run blast
    if not os.path.exists('1-selfblast/selfblastdb.results.tsv'):
        os.popen('blastn -db 1-selfblast/selfblastdb.fasta -query 1-selfblast/selfblastdb.fasta -out 1-selfblast/selfblastdb.results.tsv -num_threads 100 -outfmt "6 qseqid sseqid staxids sscinames scomnames sskingdoms pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp" -perc_identity 100 -qcov_hsp_perc 98').read()


    # now calculate LCA with my script
    os.popen('python computeLCA.py 1-selfblast/selfblastdb.results.tsv > 1-selfblast/selfblastdb.LCAs.tsv').read()
    if not os.path.exists('2-LCAs/'):
        os.mkdir('2-LCAs/')

    # now find the weirdos

    weird_queries = set()

    with open('1-selfblast/selfblastdb.LCAs.tsv') as fh:
        for line in fh:
            ll = line.rstrip().split('\t')
            lineage = ll[1]
            query = ll[0]
            # "{k};{p};{c};{o};{f};{g};{s}"
            # cellular organisms;Eukaryota;Opisthokonta;Metazoa;Eumetazoa;Bilateria;Deuterostomia;Chordata;Craniata;Vertebrata;Gnathostomata;Chondrichthyes;Elasmobranchii;Selachii;Galeomorphii;Galeoidea;Carcharhiniformes;Triakidae;Mustelus;Mustelus manazo 
            lineage_split = lineage.split(';')

            family, genus, species = lineage_split[-3:]
            if not family and not genus:
                weird_queries.add(query)

    queries_to_weird_subjects = defaultdict(set) # key: weird query, subject: set of all hits for this query
    subjects_to_species = {}
    weird_subjects_to_queries = defaultdict(set)

    with open('1-selfblast/selfblastdb.results.tsv') as fh:
        for line in fh:
            ll = line.split('\t')
            q = ll[0]
            s = ll[1]
            taxid, spec = ll[2], ll[3]
            if q in weird_queries:
                queries_to_weird_subjects[q].add(s)
                weird_subjects_to_queries[s].add(q)
                subjects_to_species[s] = (taxid, spec)


    with open('2-LCAs/TEMP_TAX.txt', 'w') as out:
        for seq_id in subjects_to_species:
            taxid, spec = subjects_to_species[seq_id]
            out.write(f'{seq_id}\t{taxid}\t{spec}\n')

    seq_to_fam = {}

    os.popen('taxonkit lineage -i 2 2-LCAs/TEMP_TAX.txt | taxonkit reformat -i 4 > 2-LCAs/TEMP_TAX_LINEAGE.txt').read()
    with open('2-LCAs/TEMP_TAX_LINEAGE.txt') as fh:
        for line in fh:
            #gb|KC136576.1|  334875  Hyperoglyphe antarctica cellular organisms;Eukaryota;Opisthokonta;Metazoa;Eumetazoa;Bilateria;Deuterostomia;Chordata;Craniata;Vertebrata;Gnathostomata;Teleostomi;Euteleostomi;Actinopterygii;Actinopteri;Neopterygii;Teleostei;Osteoglossocephalai;Clupeocephala;Euteleosteomorpha;Neoteleostei;Eurypterygia;Ctenosquamata;Acanthomorphata;Euacanthomorphacea;Percomorphaceae;Pelagiaria;Scombriformes;Centrolophidae;Hyperoglyphe;Hyperoglyphe antarctica Eukaryota;Chordata;Actinopteri;Scombriformes;Centrolophidae;Hyperoglyphe;Hyperoglyphe antarctica
            ll = line.rstrip().split('\t')
            seq_id, tax_id, spec_id, lineage, reform_lineage = ll
            # "{k};{p};{c};{o};{f};{g};{s}"
            kingdom, phylum, thisclass, order, family, genus, species = reform_lineage.split(';')
            seq_to_fam[seq_id] = family
        
    with open('2-LCAs/TEMP_AGAIN.txt', 'w') as out:
        for line in open('1-selfblast/selfblastdb.fasta'):
            if line.startswith('>'):
                ll = line.lstrip('>').split(' ')
                thisid, thistaxid, thisgenus, thisspec = ll[0], ll[1], ll[2], ll[3]
                out.write(f'{thisid}\t{thisgenus} {thisspec}\n')

    os.popen('taxonkit name2taxid -i 2 2-LCAs/TEMP_AGAIN.txt | taxonkit lineage -i 3 > 2-LCAs/TEMP_AGAIN_LINEAGE.txt').read()
    query_to_fam = {}
    for line in open('2-LCAs/TEMP_AGAIN_LINEAGE.txt'):
        ll = line.rstrip().split('\t')
        if len(ll) == 2:
            continue
        seq_id, species_name,taxid,lineage = ll
        # "{k};{p};{c};{o};{f};{g};{s}"
        family = 'NA'
        for element in lineage.split(';'):
            if element.endswith('dae'):
                family = element
        query_to_fam[seq_id] = family

    family_count = Counter()
    seen_iffys = set()
    with open('2-LCAs/Iffy_sequences.tsv', 'w') as out:
        out.write('Sequence ID\tListed family\tShould-be family\n')
        for query in queries_to_weird_subjects:
            all_q = queries_to_weird_subjects[query]
            try:
                x = query_to_fam[query]
            except:
                continue
            family_count = Counter()
            fam_to_seq = defaultdict(list)
            for q in all_q:
                family_count[seq_to_fam[q]] += 1
                fam_to_seq[seq_to_fam[q]].append( q)
            most_common = family_count.most_common(100)
            iffy_ones = most_common[1:]
            majority = most_common[0]
            # [('Scombridae', 5), ('Balistidae', 1)]
            # becomes Balistidae alone
            # which sequence is this?
            iffy_sequences = []
            for i in iffy_ones:
                fam = i[0]
                iffy_sequences += fam_to_seq[fam]
            #print(iffy_sequences, iffy_ones, 'should be ', majority)
            # ['gb|HQ592314.1|', 'gb|HQ592315.1|', 'gb|HQ592313.1|'] [('Ophidiidae', 3)] should be  ('Centrolophidae', 5)
            for i in iffy_sequences:
                if i not in seen_iffys:
                    if '|' in i:
                        shorti = i.split('|')[1]
                    else:
                        shorti = i
                    removed[shorti] = f'Placed as {seq_to_fam[i]}, might be {majority[0]}'
                    seen_iffys.add(i)

    if not os.path.exists('3-Final'):
        os.mkdir('3-Final')

    with open('3-Final/Final_database.fasta', 'w') as outfasta, \
            open('3-Final/QC_Removal_stats.tsv', 'w') as out,\
            open('3-Final/Final_database_taxids.txt', 'w') as outtaxa:
        for i in removed:
            out.write(f'{i}\t{removed[i]}\n')
        for seq in SeqIO.parse('1-selfblast/selfblastdb.fasta', 'fasta'):
            if seq.id in removed:
                continue
            thisid, thistaxid, *bla = seq.description.split(' ')
            outtaxa.write(f'{thisid} {thistaxid}\n')
            outfasta.write(seq.format('fasta'))
    
    # now format the BLAST database
    os.popen('makeblastdb -dbtype nucl -in 3-Final/Final_database.fasta -parse_seqids -taxid_map 3-Final/Final_database_taxids.txt').read()

    # now make stats for removal
    for i in removed:
        removed_classes[removed[i].split('(')[0]] += 1

    print('Final output is in "3-Final/Final_database.fasta".')
    print('Final stats are in "3-Final/QC_Removal_stats.tsv".')
    print('Summary stats are in "3-Final/QC_Removal_summary_stats.tsv".')
    print(f'Overall stats: removed {sum(removed_classes.values())} sequences.')
    with open("3-Final/QC_Removal_summary_stats.tsv", 'w') as out:
        for c in removed_classes:

            out.write(f'{c}\t{removed_classes[c]}\n')
