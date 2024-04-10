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
        splice="filtered_splice_sites.csv",
        polyA="filtered_polyA_sites.csv",
        gff="fasta/genome.gff",
        fasta="fasta/genome.fasta"
    output:
        "assigned_sites.gff",
    log: "logs/site_assignment.log"
    conda: "../envs/geffa.yaml"
    script: "../scripts/assign_sites_to_gene_models.py"

rule identify_UTRs:
    input:
        gff="assigned_sites.gff",
    output:
        "UTRs_identified.gff",
    log: "logs/UTR_identification.log"
    conda: "../envs/geffa.yaml"
    script: "../scripts/identify_UTRs.py"
