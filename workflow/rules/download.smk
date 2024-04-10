rule download_ERA_fastq_files:
    output:
        r"reads/ERR{acc_p1,\d{3}}{acc_p2,\d{3}}_{run,\d}.fastq",
    threads: 2
    resources:
        era_downloads=1
    shell:
        "curl -sS ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR{wildcards.acc_p1}/ERR{wildcards.acc_p1}{wildcards.acc_p2}/ERR{wildcards.acc_p1}{wildcards.acc_p2}_{wildcards.run}.fastq.gz | gunzip - > {output[0]}"

rule download_SRA_fastq_files:
    output:
        r"reads/{accession,SRR\d{8}}_1.fastq",
        r"reads/{accession,SRR\d{8}}_2.fastq",
    log:
        "logs/fastq/{accession}.log"
    params:
        extra="--skip-technical --split-files"
    threads: 6
    wrapper:
        "v1.32.0/bio/sra-tools/fasterq-dump"

rule download_tritrypdb_genome:
    output:
        r"fasta/tritrypdb-{organism,[^\-]+}-{version,[0-9\.]+}.fasta",
    threads: 1
    shell:
        "curl -o {output[0]} \"https://tritrypdb.org/common/downloads/release-{wildcards.version}/{wildcards.organism}/fasta/data/TriTrypDB-{wildcards.version}_{wildcards.organism}_Genome.fasta\""

rule download_ncbi_genome:
    output:
        r"fasta/{accession,GC[A-Z]_\d{9}\.\d}}.fasta",
    threads: 1
    shell:
        "curl -OJX GET \"https://api.ncbi.nlm.nih.gov/datasets/v1/genome/accession/{wildcards.accession}/download?include_annotation_type=GENOME_FASTA&filename=genome.zip\" -H \"Accept: application/zip\";"
        "unzip -j -p genome.zip '**/*.fna' > {output[0]}"

rule download_tritrypdb_gff:
    output:
        r"fasta/tritrypdb-{organism,[^\-]+}-{version,[0-9\.]+}.gff",
    threads: 1
    shell:
        "curl -o {output[0]} \"https://tritrypdb.org/common/downloads/release-{wildcards.version}/{wildcards.organism}/gff/data/TriTrypDB-{wildcards.version}_{wildcards.organism}.gff\""

rule link_download:
    input:
        f"fasta/{config['genome-accession']}.{{extension}}"
    output:
        "fasta/genome.{extension}"
    threads: 1
    shell:
        "ln {input[0]} {output[0]}"
