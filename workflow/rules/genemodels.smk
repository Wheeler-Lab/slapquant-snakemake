rule splice_polyA_gff:
    input:
        splice_sites="filtered_splice_sites.csv",
        polyA_sites="filtered_polyA_sites.csv",
        genome="fasta/genome.fasta",
    output:
        "splice_polyA.gff",
    conda: "../envs/pandas.yaml"
    script: "../scripts/create_splice_polyA_gff.py"

rule assign_sites_to_gene_models:
    input:
        sites="filtered_{type}_sites.csv",
        gff="fasta/genome.gff",
    output:
        "{type}_sites.gff",
    run:
        import geffa
        import numpy as np
        gff = geffa.GffFile(input.gff, postpone_validation=True, ignore_unknown_feature_types=True)
        cdss = {}
        for seqreq in gff.sequence_regions.values():
            entry = {'+': [], '-': []}
            for cds in (feature for feature in seqreq.node_registry.values() if feature.type == 'CDS'):
                entry[cds.strand].append(cds)
            for strand in ['+', '-']:
                if entry[strand]:
                    cds_pos, cds_nodes = zip(*sorted([(cds.end if strand == '+' else cds.start, cds) for cds in entry[strand]],
                                                    key=lambda x: x[0]))
                    entry[strand] = (np.array(cds_pos)-1, cds_nodes)
            cdss[seqreq.name] = entry

        sites = pd.read_csv(input.sites)
        counts = defaultdict(lambda: 0)
        for i, (_, row) in enumerate(sites.iterrows()):
            pos, nodes = cdss[row.Chromosome][row.Strand]
            if ((row.Strand == '+') and row.Position > pos[-1]) or ((row.Strand == '-') and row.Position < pos[0]):
                cds = None
            else:
                closest_end_idx = np.searchsorted(pos, row.Position) - (1*(row.Strand == '-'))
                cds = nodes[closest_end_idx]
                gene = cds.parents[0].parents[0]
                seqreg = gene.sequence_region
                gene_id = gene.attributes['ID']
                counts[gene_id] += 1
                slas = geffa.geffa.SLASNode(-1, seqreg, 'RNASeq', 'SLAS', row.Position+1, row.Position+2, '.', row.Strand, '.', f'ID=slas_{gene_id}-{counts[gene_id]};Parent={gene_id}')
        gff.save(output[0])