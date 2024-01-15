from Bio import SeqIO
lineage_dict = {}
tax_to_species = {}
for line in open('species_present_with_taxid.txt'):
    ll = line.rstrip().split('\t')
    try:
        spec, taxid = ll
    except:
        # no taxid
        continue

    lineage_dict[spec] = taxid
    tax_to_species[taxid] = spec

# now get the synonyms too
for line in open('/home/pbayer/.taxonkit/names.dmp'):
    ll = [x.strip() for x in line.split('|')]
    if ll[-2] == 'synonym':
        taxid, species = ll[0], ll[1]
        lineage_dict[species] = taxid
#/home/pbayer/.taxonkit/names.dmp:2599772        |       Himantura gerrardii     |               |       synonym

with open('mitogenomes_fish_nuccore.renamedFiltered.fasta', 'w') as out1, open('mitogenomes_fish_nuccore.renamedFiltered.taxids.txt', 'w') as out2:
    for seq in SeqIO.parse('mitogenomes_fish_nuccore.fasta', 'fasta'):
        if ' sp. ' in seq.description or \
                ' aff. ' in seq.description or \
                ' x ' in seq.description:
            continue
        thisid = seq.description
        species = ' '.join(thisid.split(' ')[1:3])

        try:
            taxid = lineage_dict[species]
        except:
            print('No taxid? for %s'%(species))
            continue
        newspec = tax_to_species[taxid]
        if newspec != species:
            print(newspec, species)
        # now let's reformat the name of this.
        # old name: NCBI ID {SPACE} All info we have
        # new name: NCBI ID {SPACE} Taxonomy ID { SPACE } Species {SPACE} All info we have
        newid = f'{seq.id} {taxid} {newspec} {seq.name}'
        seq.id = newid
        out1.write(seq.format('fasta'))

        out2.write(f'{seq.id} {taxid}\n')
