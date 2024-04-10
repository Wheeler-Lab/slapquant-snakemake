import itertools
import geffa
import numpy as np
import pandas as pd
from collections import defaultdict
import sys

import logging

logger = logging.getLogger('site_assignment')
logging.basicConfig(filename=snakemake.log[0])

gff = geffa.GffFile(snakemake.input.gff, fasta_file=snakemake.input.fasta, postpone_validation=True, ignore_unknown_feature_types=True)

sites = pd.concat([
    pd.read_csv(snakemake.input.splice).assign(type='SLAS'),
    pd.read_csv(snakemake.input.polyA).assign(type='PAS'),
]).sort_values(['Chromosome', 'Strand', 'Position'])

counters = {
    'SLAS': itertools.count(1),
    'PAS': itertools.count(1),
}

for (contig, strand), g in sites.groupby(['Chromosome', 'Strand']):
    try:
        seqreg = gff.sequence_regions[contig]
    except KeyError:
        logger.warning(f"No sequence region with name {contig} found in the gene models given. Discarding the SLAS and PAS nodes.")
        continue

    for _, row in g.iterrows():
        if row.type == 'SLAS':
            node = geffa.geffa.SLASNode
        elif row.type == 'PAS':
            node = geffa.geffa.PASNode
        else:
            raise ValueError(f'Invalid node type {row.type}')
        usage = row.iloc[3:-1].sum()
        node(-1, seqreg, 'RNASeq', row.type, row.Position+1, row.Position+2, '.', strand, '.', f'ID={row.type}_{next(counters[row.type])};Usage={usage}')

for seqreg in gff.sequence_regions.values():
    for slas in (node for node in seqreg.node_registry.values() if node.type == 'SLAS'):
        nodes_to_search = sorted(
            (
                feature for feature in seqreg.node_registry.values()
                if (
                    feature.type in ['CDS', 'PAS'] and
                    (
                        (
                            (slas.strand == '+') and (feature.strand == '+') and (feature.start >= slas.end) or
                            (slas.strand == '-') and (feature.strand == '-') and (feature.end <= slas.start)
                        )
                    )
                )
            ),
            key=lambda x: x.start if slas.strand == '+' else -x.end
        )
        if len(nodes_to_search) == 0:
            logger.warning(f"No CDS found we could match {slas.attributes['ID']} to. Skipping.")
            continue

        closest_node = nodes_to_search[0]

        if closest_node.type == 'PAS':
            logger.warning(f"Closest node to {slas.attributes['ID']} is a PAS, no CDS could be assigned.")
        else:
            slas.add_parent(closest_node.parents[0].parents[0])

    for pas in (node for node in seqreg.node_registry.values() if node.type == 'PAS'):
        nodes_to_search = sorted(
            (
                feature for feature in seqreg.node_registry.values()
                if (
                    feature.type in ['CDS', 'SLAS'] and
                    (
                        (
                            (pas.strand == '+') and (feature.strand == '+') and (feature.end <= pas.start) or
                            (pas.strand == '-') and (feature.strand == '-') and (feature.start >= pas.end)
                        )
                    )
                )
            ),
            key=lambda x: -x.end if pas.strand == '+' else x.start
        )
        if len(nodes_to_search) == 0:
            logger.warning(f"No CDS found we could match {pas.attributes['ID']} to. Skipping.")
            continue

        closest_node = nodes_to_search[0]

        if closest_node.type == 'SLAS':
            logger.warning(f"Closest node to {pas.attributes['ID']} is a SLAS, no CDS could be assigned.")
        else:
            pas.add_parent(closest_node.parents[0].parents[0])

gff.save(snakemake.output[0])