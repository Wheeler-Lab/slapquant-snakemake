# Run faidx to create an index for the genome FASTA
rule index_genome:
    input:
        "fasta/genome.fasta",
    output:
        "fasta/genome.fasta.fai",
    log:
        "fasta/genome_index.log",
    threads: 1
    params:
        extra="",
    wrapper:
        "v1.17.1/bio/samtools/faidx"

# Let bowtie2 also create an index for the genome FASTA
rule build_mapping_database:
    input:
        ref="fasta/genome.fasta"
    output:
        multiext(
            "mapping_database/genome",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    log:
        "logs/bowtie2_build/build.log",
    params:
        extra="",
    threads: 8
    wrapper:
        "v1.17.0/bio/bowtie2/build"

# # Trim the reads for a given sample
# rule trimmomatic_pe:
#     input:
#         r1="reads/{sample_id}_1.fastq",
#         r2="reads/{sample_id}_2.fastq"
#     output:
#         r1="trimmed/{sample_id}_1.fastq",
#         r2="trimmed/{sample_id}_2.fastq",
#         r1_unpaired="trimmed/{sample_id}_1.unpaired.fastq",
#         r2_unpaired="trimmed/{sample_id}_2.unpaired.fastq"
#     log:
#         "logs/trimmomatic/{sample_id}.log"
#     params:
#         trimmer=[
#             "ILLUMINACLIP:../resources/adaptersall_adaptors.fasta:2:30:10",
#             "LEADING:10",
#             "TRAILING:10",
#             "SLIDINGWINDOW:5:15",
#             "MINLEN:50"
#         ],
#         extra="",
#         compression_level="-9"
#     threads:
#         4
#     resources:
#         mem_mb=24*1024
#     wrapper:
#         "v3.7.0/bio/trimmomatic/pe"

rule bbduk:
    input:
        ["reads/{sample_id}_1.fastq", "reads/{sample_id}_2.fastq"],
        adapters=workflow.source_path("../../resources/adaptersall_adaptors.fasta"),
    output:
        ["trimmed/{sample_id}_1.fastq", "trimmed/{sample_id}_2.fastq"],
        stats="stats/trimming_stats/{sample_id}.txt",
    log:
        "logs/bbduk/{sample_id}.log",
    threads: 4
    params:
        command="bbduk.sh",
        extra="tpe tbo ",
        ktrim="r ",
        ref=lambda w, input: [input.adapters, "adapters", "artifacts"],
        k=23,
        mink=11,
        hdist=1,
        trimpolygright=10,
        minlen=25,
        maxns=30,
        entropy=0.5,
        entropywindow=50,
        entropyk=5,
    resources:
        mem_mb=4000,
    wrapper:
        "v3.7.0/bio/bbtools"

# Join paired reads for a given sample
rule fastq_join:
    input:
        r1="trimmed/{sample_id}_1.fastq",
        r2="trimmed/{sample_id}_2.fastq",
    output:
        multiext(
            "joined/{sample_id}",
            ".join",
            ".un1",
            ".un2",
        ),
    conda:
        "../envs/ea-utils.yaml"
    resources:
        fastq_join_processes=1
    shell:
        "fastq-join {input.r1} {input.r2} -o joined/{wildcards.sample_id}."

# Map / align filtered reads to the genome using bowtie2
rule map_reads:
    input:
        sample=["{read_type}/{sample_id}.fastq"],
        idx=multiext(
            "mapping_database/genome",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    output:
        "mapped_{read_type}/{sample_id}.bam",
    log:
        "logs/bowtie2/{read_type}_{sample_id}.log",
    params:
        extra="",
    threads: 8  # Needs to be >= 2 because it runs bowtie2 and samtools
    wrapper:
        "v1.17.0/bio/bowtie2/align"