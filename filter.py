#!/usr/bin/env python3
import argparse
import pysam

cli = argparse.ArgumentParser()
cli.add_argument('refr')
cli.add_argument('bam')
cli.add_argument('newbam')
args = cli.parse_args()

extent_spans = dict()
with open(args.refr, 'r') as fh:
    for line in fh:
        if not line.startswith('>'):
            continue
        locusid, refrloc, varinfo = line[1:].strip().split()
        varoffsets = varinfo.split('=')[1]
        varloc = [int(x) for x in varoffsets.split(':')]
        varmin, varmax = min(varloc), max(varloc)
        extent_spans[locusid] = (varmin, varmax)

bam = pysam.AlignmentFile(args.bam, 'rb')
newbam = pysam.AlignmentFile(args.newbam, 'wb', template=bam)
kept = 0
discarded = 0
for record in bam:
    varmin, varmax = extent_spans[record.reference_name]
    if record.reference_start > varmin or record.reference_end < varmax:
        discarded += 1
        continue
    kept += 1
    _ = newbam.write(record)

print('Kept', kept, 'records, discarded', discarded)
