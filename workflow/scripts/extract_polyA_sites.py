import sys
import pysam
print(pysam.__version__)
from collections import defaultdict

# Helper function to reverse complement a sequence
complement_translation = str.maketrans("ACGT", "TGCA")
def reverse_complement(seq):
    return seq[::-1].translate(complement_translation)

def polyA_entry_factory():
    return defaultdict(lambda: {'count': 0})

def extract_polyA_sites(bam_file, genome_fasta, genome_index, output_file):
    polyA_sites = defaultdict(polyA_entry_factory)

    samfile = pysam.AlignmentFile(bam_file, 'rb')
    fastafile = pysam.FastaFile(genome_fasta, genome_index)
    for entry in samfile:
        if (entry.reference_name is None) or (entry.query_qualities is None) or (entry.reference_name == '*'):
            continue
        # if (entry.cigartuples is None) or (entry.query_length != dict(entry.cigartuples)[0]) or (entry.mapping_quality == 0):
        #     continue
        nr_A = int(entry.query_name.split('/')[0].split('_')[-1])
        if entry.is_forward:
            if entry.reference_end + nr_A > fastafile.get_reference_length(entry.reference_name):
                continue
            post_seq = fastafile.fetch(entry.reference_name, entry.reference_end, entry.reference_end + nr_A)
            position = entry.reference_end
        else:
            if entry.reference_start - nr_A < 0:
                continue
            post_seq = reverse_complement(fastafile.fetch(entry.reference_name, entry.reference_start - nr_A, entry.reference_start))
            position = entry.reference_start
        
        count_A = post_seq.count('A')
        if count_A < nr_A / 2:
            site = polyA_sites[entry.reference_name][(position, '+' if entry.is_forward else '-')]
            site['count'] += 1

    site_counts = {'+': 0, '-': 0}
    with open(output_file, 'wt') as outfile:
        outfile.write('Chromosome\tPosition\tCount\tStrand\n')
        for chromosome in sorted(polyA_sites):
            chromosome_polyA_sites = polyA_sites[chromosome]
            for position, strand in sorted(chromosome_polyA_sites, key=lambda x: x[0]):
                site = chromosome_polyA_sites[(position, strand)]
                site_counts[strand] += 1
                outfile.write(f'{chromosome}\t{position}\t{site["count"]}\t{strand}\n')
    print(f'{site_counts["+"]} forward sites, {site_counts["-"]} reverse sites')

extract_polyA_sites(snakemake.input.bam, snakemake.input.genome, snakemake.input.genome_idx, snakemake.output[0])
# import sys
# extract_polyA_sites(*sys.argv[1:])