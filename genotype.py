#!/usr/bin/env python3
import argparse
from collections import defaultdict
import pysam
import sys

cli = argparse.ArgumentParser()
cli.add_argument('refr')
cli.add_argument('bam')
args = cli.parse_args()


varlocs = dict()
with open(args.refr, 'r') as fh:
    for line in fh:
        if not line.startswith('>'):
            continue
        locusid, refrloc, varinfo = line[1:].strip().split()
        varoffsets = varinfo.split('=')[1]
        varloc = [int(x) for x in varoffsets.split(':')]
        varlocs[locusid] = varloc

discarded = 0
bam = pysam.AlignmentFile(args.bam, 'rb')
for locusid in sorted(varlocs):
    genotypes = defaultdict(int)
    gt = defaultdict(dict)
    varloc = set(varlocs[locusid])
    cov_pos = list()
    for column in bam.pileup(locusid):
        cov_pos.append(column.n)
        if column.pos not in varloc:
            continue
        for record in column.pileups:
            aligned_base = None
            if record.is_del or record.is_refskip:
                continue
            aligned_base = record.alignment.query_sequence[record.query_position]
            gt[record.alignment.query_name][column.pos] = aligned_base
    for readname, gtdict in gt.items():
        gtlist = [gtdict[pos] for pos in sorted(gtdict)]
        if len(gtlist) < len(varloc):
            discarded += 1
            continue
        gtstr = ','.join(gtlist)
        genotypes[gtstr] += 1

    cov_avg = sum(cov_pos) / len(cov_pos)
    cov_min = min(cov_pos)
    cov_max = max(cov_pos)
    print('Locus: {loc:s}'.format(loc=locusid), end='; ', file=sys.stderr)
    print('Coverage: avg={avg:.1f} min={mn:d} max={mx:d}'.format(avg=cov_avg, mn=cov_min, mx=cov_max), file=sys.stderr)

    avgcount = sum(genotypes.values()) / len(genotypes.values())
    for gtstr, gtcount in sorted(genotypes.items()):
        if gtcount < avgcount:
            continue
        print(locusid, '{gt:s}\t{n:d}'.format(gt=gtstr, n=gtcount))

print('Discarded', discarded, 'reads with gaps or missing data at positions of interest', file=sys.stderr)
