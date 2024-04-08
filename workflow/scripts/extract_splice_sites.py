import pysam
from collections import defaultdict

# Helper function to reverse complement a sequence
complement_translation = str.maketrans("ACGT", "TGCA")
def reverse_complement(seq):
    return seq[::-1].translate(complement_translation)

def splice_entry_factory():
    return defaultdict(lambda: {'count': 0})

def extract_splice_sites(bam_file, genome_fasta, genome_index, output_file, splice_fragment):
    splice_sites = defaultdict(splice_entry_factory)

    samfile = pysam.AlignmentFile(bam_file, 'rb')
    fastafile = pysam.FastaFile(genome_fasta, genome_index)
    for entry in samfile:
        if (entry.query_qualities is None) or (entry.reference_name == '*'):
            continue
        if (entry.cigartuples is None) or (entry.query_length != dict(entry.cigartuples)[0]) or (entry.mapping_quality == 0):
            continue
        if entry.is_forward:
            if entry.reference_start < len(splice_fragment):
                continue
            pre_seq = fastafile.fetch(entry.reference_name, entry.reference_start-len(splice_fragment), entry.reference_start)
            position = entry.reference_start-2
        else:
            if entry.reference_end > fastafile.get_reference_length(entry.reference_name) - len(splice_fragment):
                continue
            pre_seq = fastafile.fetch(entry.reference_name, entry.reference_end, entry.reference_end+len(splice_fragment))
            pre_seq = reverse_complement(pre_seq)
            position = entry.reference_end
        
        splice_frag_similarity = sum([splice_fragment[i] == pre_seq[i] for i in range(len(splice_fragment))])
        if splice_frag_similarity < len(splice_fragment) - 2:
            dinucleotide = pre_seq[-2:]
            site = splice_sites[entry.reference_name][(position, '+' if entry.is_forward else '-')]
            site['count'] += 1
            site['dinucleotide'] = dinucleotide

    site_counts = {'+': 0, '-': 0}
    with open(output_file, 'wt') as outfile:
        outfile.write('Chromosome\tPosition\tCount\tStrand\tDinucleotide\n')
        for chromosome in sorted(splice_sites):
            chromosome_splice_sites = splice_sites[chromosome]
            for position, strand in sorted(chromosome_splice_sites, key=lambda x: x[0]):
                site = chromosome_splice_sites[(position, strand)]
                site_counts[strand] += 1
                outfile.write(f'{chromosome}\t{position}\t{site["count"]}\t{strand}\t{site["dinucleotide"]}\n')
    print(f'{site_counts["+"]} forward sites, {site_counts["-"]} reverse sites')

extract_splice_sites(snakemake.input.bam, snakemake.input.genome, snakemake.input.genome_idx, snakemake.output[0], snakemake.config['splice_leader_sequence'])
