import microhapdb
import microhapulator
import pyfaidx

index = pyfaidx.Fasta('../microhapulator/hg38.fasta')
l = microhapulator.locus.default_panel()
loci = microhapdb.loci[microhapdb.loci.ID.isin(l)]

with open('default-panel.fasta', 'w') as fh:
    for i, row in loci.iterrows():
        c = microhapulator.LocusContext(row)
        coords = [c.global_to_local(x) for x in microhapdb.allele_positions(row.ID)]
        print('>', c.defline(), ' variants=', ':'.join(map(str, coords)), '\n', c.sequence(index), sep='', file=fh)
