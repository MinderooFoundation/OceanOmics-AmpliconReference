from Bio import SeqIO
import os
from collections import Counter

twelve = [x.id for x in SeqIO.parse('12S/3-Final/Final_database.fasta', 'fasta')]
sixteen = [x.id for x in SeqIO.parse('16S/3-Final/Final_database.fasta', 'fasta')]

taxids = {}
for line in open('12S/3-Final/Final_database_taxids.txt'):
    ll = line.split()
    name, taxid = ll
    taxids[name] = taxid

for line in open('16S/3-Final/Final_database_taxids.txt'):
    ll = line.split()
    name, taxid = ll
    taxids[name] = taxid

c = Counter(twelve + sixteen)
dups = set()
for i in c:
    if c[i] != 1:
        dups.add(i)

with open('12S.16S.fasta', 'w') as out, \
        open('12S.16S.taxids.txt', 'w') as taxout:

    for t in SeqIO.parse('12S/3-Final/Final_database.fasta', 'fasta'):
        if t.id in dups and t.id.startswith('NC'):
            # keep this but rename
            old_id =  t.id
            t.id = '12S_' + t.id
            t.description = '12S_' + t.description
            t.name = '12S_' + t.name
            out.write(t.format('fasta'))
            taxout.write(f'{t.id} {taxids[old_id]}\n')
        elif t.id in dups and not t.id.startswith('NC'):
            # keep this - it has 12S and 16S together
            out.write(t.format('fasta'))
            taxout.write(f'{t.id} {taxids[t.id]}\n')
        else:

            oldid = t.id
            if t.id.startswith('NC'):
                # avoid duplication with whole mitogenomes
                t.id = '12S_%s'%(t.id)
            # not a duplicate, keep
            taxout.write(f'{t.id} {taxids[oldid]}\n')
            out.write(t.format('fasta'))

    for t in SeqIO.parse('16S/3-Final/Final_database.fasta', 'fasta'):
        if t.id in dups and t.id.startswith('NC'):
            # keep this but rename
            old_id = t.id
            t.id = '16S_' + t.id
            t.description = '16S_' + t.description
            t.name = '16S_' + t.name
            out.write(t.format('fasta'))
            taxout.write(f'{t.id} {taxids[old_id]}\n')
        elif t.id in dups and not t.id.startswith('NC'):
            # skip this - we kept it with 12S
            continue
        else:
            # not a duplicate, keep
            oldid = t.id
            if t.id.startswith('NC'):
                t.id = '16S_%s'%(t.id)
            out.write(t.format('fasta'))
            taxout.write(f'{t.id} {taxids[oldid]}\n')


os.popen('cat 12S.16S.taxids.txt Mitogenomes/mitogenomes_fish_nuccore.renamedFiltered.taxids.txt > 12S.16S.Mitogenomes.taxids.txt')
os.popen('cat 12S.16S.fasta Mitogenomes/mitogenomes_fish_nuccore.renamedFiltered.fasta > 12S.16S.Mitogenomes.fasta')
