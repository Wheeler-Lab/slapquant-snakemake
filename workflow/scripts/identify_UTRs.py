import itertools
import geffa
import numpy as np
import pandas as pd
from collections import defaultdict
import sys

import logging

logger = logging.getLogger('UTR_identification')
logging.basicConfig(filename=snakemake.log[0])

gff = geffa.GffFile(snakemake.input.gff, postpone_validation=True, ignore_unknown_feature_types=True)

for seqreg in gff.sequence_regions.values():
    for gene in [feature for feature in seqreg.node_registry.values() if feature.type == 'gene']:
        mRNAs = [feature for feature in gene.children if feature.type == 'mRNA']
        if len(mRNAs) == 0:
            logger.warning(f"{gene.attributes['ID']} is not a protein coding gene, skipping UTR assignment.")
            continue
        elif len(mRNAs) > 1:
            logger.warning(f"{gene.attributes['ID']} has multiple mRNAs assigned, UTR assignment isn't implemented yet.")
            continue
        mRNA = mRNAs[0]
        CDSs = mRNA.CDS_children()

        slas_sites = [feature for feature in gene.children if feature.type == 'SLAS']
        pas_sites = [feature for feature in gene.children if feature.type == 'PAS']

        if slas_sites:
            slas = sorted(slas_sites, key=lambda x: -int(x.attributes['Usage']))[0]
            CDS = CDSs[0]
            if mRNA.strand == '+':
                start = slas.end
                end = CDS.start
            else:
                end = slas.start
                start = CDS.end
            UTR5 = geffa.geffa.FivePrimeUTRNode(-1, seqreg, 'RNASeq', 'five_prime_UTR', start, end, '.', mRNA.strand, '.', f'ID={mRNA.attributes["ID"]}_UTR5;Parent={mRNA.attributes["ID"]}')
        if pas_sites:
            pas = sorted(pas_sites, key=lambda x: -int(x.attributes['Usage']))[0]
            CDS = CDSs[-1]
            if mRNA.strand == '+':
                start = CDS.end
                end = pas.start
            else:
                end = CDS.start
                start = pas.end
            UTR3 = geffa.geffa.ThreePrimeUTRNode(-1, seqreg, 'RNASeq', 'three_prime_UTR', start, end, '.', mRNA.strand, '.', f'ID={mRNA.attributes["ID"]}_UTR3;Parent={mRNA.attributes["ID"]}')
    break

gff.save(snakemake.output[0])