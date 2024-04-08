import pathlib
import pandas as pd

sites = pd.concat([pd.read_table(fn).assign(sample=pathlib.Path(fn).stem) for fn in snakemake.input])
sites = sites.groupby(['Chromosome', 'Position', 'Strand', 'sample']).Count.sum().reset_index()
sites = sites.set_index(['Chromosome', 'Position', 'Strand', 'sample']).loc[:, ['Count']].unstack('sample', 0)
sites.droplevel(0, 1).reset_index().to_csv(snakemake.output[0], index=False)
