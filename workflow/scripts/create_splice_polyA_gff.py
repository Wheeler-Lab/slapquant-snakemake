import pandas as pd

entries = []
for site, fn in [('SLAS', snakemake.input.splice_sites), ('PAS', snakemake.input.polyA_sites)]:
    df = pd.read_csv(fn)
    keys = df.columns[3:]
    for i, (_, row) in enumerate(df.iterrows()):
        attributes = dict(zip(keys, row.iloc[3:]))
        attributes['ID'] = f'{site}_{row.Chromosome}-{i}'
        attributes_str = ';'.join([f'{k}={v}' for k, v in attributes.items()])
        entry_str = f'{row.Chromosome}\tUD1\t{site}\t{row.Position+1}\t{row.Position+2}\t.\t{row.Strand}\t.\t{attributes_str}\n##\n'
        entries.append((row.Chromosome, row.Position, entry_str))
with open(snakemake.output[0], 'wt') as output:
    for _, _, entry in sorted(entries, key=lambda e: e[:-1]):
        output.write(entry)
    output.write('##FASTA\n')
    with open(snakemake.input.genome, 'rt') as genome:
        for l in genome:
            output.write(l)