import pandas as pd

threshold = len(snakemake.config["rnaseq-runs"])        
df = pd.read_csv(snakemake.input[0]).set_index(['Chromosome', 'Position', 'Strand'])
filtered = df[df.sum(axis=1) >= threshold]
filtered.reset_index().to_csv(snakemake.output[0], index=False)