# Helper function to reverse complement a sequence
complement_translation = str.maketrans("ACGT", "TGCA")
def reverse_complement(seq):
    return seq[::-1].translate(complement_translation)

# Helper function to read four lines of a FASTQ file at once
def chunk_4(iterator):
    while True:
        try:
            yield [next(iterator).strip() for _ in range(4)]
        except StopIteration:
            break

# Checks if a given read contains the splice leader sequence
def splice_check(seq, quality, minimum_fragment_length, splice_leader_sequence):
    for _ in range(2):
        idx = seq.find(splice_leader_sequence)
        if idx >= 0:
            fragment = seq[idx+len(splice_leader_sequence):]
            if (len(fragment) >= minimum_fragment_length):
                return fragment, quality[-len(fragment):]
            raise ValueError
        seq = reverse_complement(seq)
    raise ValueError

# Polyadenylation site identification regex
import re
polyA_re = re.compile(f'^(.+?)(A{{{config["minimum_polyA_length"]},}})$')

# Checks if a given read for polyadenylation
def polyA_check(seq, quality, minimum_fragment_length):
    for _ in range(2):
        m = polyA_re.match(seq)
        if m:
            fragment = m.group(1)
            if len(fragment) >= minimum_fragment_length:
                return fragment, quality[:len(fragment)], len(m.group(2))
            raise ValueError
        seq = reverse_complement(seq)
    raise ValueError

# Filter for reads containing the splice leader sequence and polyA sites
rule find_splice_polyA_reads:
    input:
        "trimmed/{sample_id}_1.fastq",
        "trimmed/{sample_id}_1.fastq",
        # multiext(
        #     "joined/{sample_id}",
        #     ".join",
        #     ".un1",
        #     ".un2",
        # ),
        # "trimmed/{sample_id}_1.unpaired.fastq",
        # "trimmed/{sample_id}_2.unpaired.fastq",
    output:
        "splice_reads/{sample_id}.fastq",
        "polyA_reads/{sample_id}.fastq",
        "other_reads/{sample_id}.fastq",
    params:
        minimum_fragment_length = config['minimum_fragment_length'],
        splice_leader_sequence = config['splice_leader_sequence']
    threads: 1
    run:
        import itertools
        counter = itertools.count(1)
        with (
            open(output[0], 'wt') as splice_outfile,
            open(output[1], 'wt') as polyA_outfile,
            open(output[1], 'wt') as other_outfile,
            open(input[0], 'rt') as infile,
        ):
            for _, seq, _, quality in chunk_4(infile):
                c = next(counter)
                has_site = False
                try:
                    fragment, fragment_quality = splice_check(seq, quality, params["minimum_fragment_length"], params["splice_leader_sequence"])
                    splice_outfile.write(f'@read_{c}/1\n{fragment}\n+\n{fragment_quality}\n')
                    has_site = True
                except ValueError:
                    pass
                try:
                    fragment, fragment_quality, nr_A = polyA_check(seq, quality, params["minimum_fragment_length"])
                    polyA_outfile.write(f'@read_{c}_{nr_A}/1\n{fragment}\n+\n{fragment_quality}\n')
                    has_site = True
                except ValueError:
                    pass
                if not has_site:
                    other_outfile.write(f'@read_{c}_{nr_A}/1\n{fragment}\n+\n{fragment_quality}\n')

# Read the mapped bam file and extract splice sites
rule read_splice_bam:
    input:
        bam="mapped_splice_reads/{sample_id}.bam",
        genome="fasta/genome.fasta",
        genome_idx="fasta/genome.fasta.fai",
    output:
        "splice_sites/{sample_id}.txt"
    threads: 1
    conda:
        "../envs/pysam.yaml"
    script: "../scripts/extract_splice_sites.py"

# Read the mapped bam file and extract polyadenylation sites
rule read_polyA_bam:
    input:
        bam="mapped_polyA_reads/{sample_id}.bam",
        genome="fasta/genome.fasta",
        genome_idx="fasta/genome.fasta.fai",
    output:
        "polyA_sites/{sample_id}.txt"
    threads: 1
    conda:
        "../envs/pysam.yaml"
    script: "../scripts/extract_polyA_sites.py"

# Combine the splice sites in each sample into one csv file
rule combine_sites:
    #input:
    #    expand('{{type}}_sites/{sample}.txt', sample=config["rnaseq-runs"])
    output:
        "{type,(splice|polyA)}_sites.csv"
    threads: 1
    conda: "../envs/pandas.yaml"
    script: "../scripts/combine_sites.py"

rule filter_sites:
    input:
        "{type,(splice|polyA)]}_sites.csv"
    output:
        "filtered_{type,(splice|polyA)}_sites.csv"
    conda: "../envs/pandas.yaml"
    script: "../scripts/filter_sites.py"
