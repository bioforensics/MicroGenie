#!/usr/bin/env python3
import argparse
import microgenie
import sys

cli = argparse.ArgumentParser()
cli.add_argument('refr')
cli.add_argument('bam')
args = cli.parse_args()

genotyper = microgenie.genotype(args.bam, args.refr)
for locusid, cov_by_pos, genotypes in genotyper:
    cov_avg = sum(cov_by_pos) / len(cov_by_pos)
    cov_min = min(cov_by_pos)
    cov_max = max(cov_by_pos)
    print('Locus: {loc:s}'.format(loc=locusid), end='; ', file=sys.stderr)
    print('Coverage: avg={avg:.1f} min={mn:d} max={mx:d}'.format(avg=cov_avg, mn=cov_min, mx=cov_max), file=sys.stderr)
    avgcount = sum(genotypes.values()) / len(genotypes.values())
    for gtstr, gtcount in sorted(genotypes.items()):
        if gtcount < avgcount:
            continue
        print(locusid, '{gt:s}\t{n:d}'.format(gt=gtstr, n=gtcount))
