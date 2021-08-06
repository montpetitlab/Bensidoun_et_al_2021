SAMPLES=['AP_A_replicat1', 'AP_A_replicat2', 'AP_A_replicat3',
         'AP_B_replicat1', 'AP_B_replicat2', 'AP_B_replicat3',
         'AP_C_replicat1', 'AP_C_replicat2', 'AP_C_replicat3',
         'AP_D_replicat1', 'AP_D_replicat2', 'AP_D_replicat3',
         'AP_E_replicat1', 'AP_E_replicat2', 'AP_E_replicat3',
         'AP_F_replicat1', 'AP_F_replicat2', 'AP_F_replicat3',
         'AP_G_replicat1', 'AP_G_replicat2', 'AP_G_replicat3']

rule all:
    input: 
        expand("outputs_5mil/fastqc_trimmed/{sample}.trimmed_fastqc.html", sample = SAMPLES),
        "outputs_5mil/fastp_trimmed/multiqc_report.html",
        expand("outputs_5mil/kmer_counts/{sample}.txt", sample = SAMPLES),
        expand("outputs_5mil/abundtrim/{sample}.abundtrim.fq.gz", sample = SAMPLES),
        expand('outputs_5mil/star/{sample}Aligned.sortedByCoord.out.bam.bai', sample = SAMPLES),
        "outputs_5mil/counts/raw_counts.tsv"

rule subsample_to_5million_reads:
    input: "{sample}.trimmed.fq.gz"
    output: "outputs_5mil/subsampled/{sample}.trimmed.fq.gz"
    conda: "envs/bbmap.yml"
    shell:'''
    reformat.sh in={input} out={output} samplereadstarget=5000000
    #seqtk sample -s100 {input} 5000000 > {output}
    #zcat {input} | head -n 20000000 | gzip > {output} 
    '''

rule fastqc:
    input: "outputs_5mil/subsampled/{sample}.trimmed.fq.gz"
    output: 
        "outputs_5mil/fastqc_trimmed/{sample}.trimmed_fastqc.html",
        "outputs_5mil/fastqc_trimmed/{sample}.trimmed_fastqc.zip"
    conda: "envs/fastqc.yml"
    shell:'''
    fastqc -o outputs_5mil/fastqc_trimmed {input}
    '''

rule fastp_trimmed_reads:
    input: "outputs_5mil/subsampled/{sample}.trimmed.fq.gz"
    output: 
        json="outputs_5mil/fastp_trimmed/{sample}.trimmed.fastp.json",
        html="outputs_5mil/fastp_trimmed/{sample}.trimmed.fastp.html",
    conda: "envs/fastp.yml"
    shell:'''
    fastp -i {input} -j {output.json} -h {output.html}
    '''

rule multiqc_fastp_trimmed:
    input: expand("outputs_5mil/fastp_trimmed/{sample}.trimmed.fastp.json", sample=SAMPLES)
    output: "outputs_5mil/fastp_trimmed/multiqc_report.html"
    params: 
        indir = "outputs_5mil/fastp_trimmed",
        outdir = "outputs_5mil/fastp_trimmed"
    conda: "envs/multiqc.yml"
    shell:'''
    multiqc {params.indir} -o {params.outdir} 
    '''

rule kmer_trim_reads:
    input: "outputs_5mil/subsampled/{sample}.trimmed.fq.gz"
    output: "outputs_5mil/abundtrim/{sample}.abundtrim.fq.gz"
    conda: "envs/khmer.yml"
    shell:'''
    trim-low-abund.py --gzip -C 3 -Z 18 -M 16e9 -V {input} -o {output}
    '''

rule kmer_countgraph:
    input: "outputs_5mil/abundtrim/{sample}.abundtrim.fq.gz"
    output: "outputs_5mil/kmer_counts/{sample}_countgraph"
    conda: "envs/khmer.yml"
    shell:'''
    load-into-counting.py {output} {input}
    '''

rule counts_kmers:
    input: 
        graph="outputs_5mil/kmer_counts/{sample}_countgraph",
        reads="outputs_5mil/abundtrim/{sample}.abundtrim.fq.gz"
    output: "outputs_5mil/kmer_counts/{sample}.txt"
    conda: "envs/khmer.yml"
    shell:'''
    count-median.py {input.graph} {input.reads} {output}
    '''

#################################
## map and count
#################################

rule download_genome:
    output: 'inputs/genome/GCF_000146045.2_R64_genomic.fna.gz'
    shell:'''
    wget -O {output} https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Saccharomyces_cerevisiae/latest_assembly_versions/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz
    '''

rule decompress_genome:
    input: 'inputs/genome/GCF_000146045.2_R64_genomic.fna.gz'
    output: 'inputs/genome/GCF_000146045.2_R64_genomic.fna'
    shell: "gunzip {input}"

rule download_gtf:
    output: 'inputs/genome/GCF_000146045.2_R64_genomic.gtf.gz'
    shell:'''
    wget -O {output} https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Saccharomyces_cerevisiae/latest_assembly_versions/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gtf.gz
    '''

rule decompress_gtf:
    input: 'inputs/genome/GCF_000146045.2_R64_genomic.gtf.gz'
    output: 'inputs/genome/GCF_000146045.2_R64_genomic.gtf'
    shell: "gunzip {input}"

rule star_index_genome:
    input:
        genome = 'inputs/genome/GCF_000146045.2_R64_genomic.fna',
        gtf = 'inputs/genome/GCF_000146045.2_R64_genomic.gtf'
    params: input_dir = 'inputs/genome' 
    output: 'inputs/genome/SAindex'
    conda: 'envs/star.yml'
    shell:'''
    STAR --runThreadN 1 --runMode genomeGenerate --genomeDir {params.input_dir} \
         --genomeFastaFiles {input.genome} --sjdbGTFfile {input.gtf} --sjdbOverhang  99
    '''

rule gunzip_reads:
    input: "outputs_5mil/subsampled/{sample}.trimmed.fq.gz"
    output: 'outputs_5mil/gunzipped/{sample}.trimmed.fq'
    shell:'''
    gunzip -c {input} > {output}
    '''

rule star_align:
    #v2.5.2a
    input:
        reads = 'outputs_5mil/gunzipped/{sample}.trimmed.fq',
        genome_index = 'inputs/genome/SAindex'
    output: 'outputs_5mil/star/{sample}Aligned.sortedByCoord.out.bam' 
    params: 
        out_prefix = lambda wildcards: 'outputs_5mil/star/' + wildcards.sample,
        genome_dir = 'inputs/genome'
    conda: 'envs/star.yml'
    shell:'''
    STAR --runThreadN 2 --genomeDir {params.genome_dir}      \
        --readFilesIn {input.reads} --outFilterType BySJout  \
        --outFilterMultimapNmax 20 --alignSJoverhangMin 8    \
        --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverLmax 0.6 --alignIntronMin 20 \
        --alignIntronMax 1000000 --alignMatesGapMax 1000000  \
        --outSAMattributes NH HI NM MD --outSAMtype BAM      \
        SortedByCoordinate --outFileNamePrefix {params.out_prefix}
    '''

rule index_bam:
    input: 'outputs_5mil/star/{sample}Aligned.sortedByCoord.out.bam'
    output: 'outputs_5mil/star/{sample}Aligned.sortedByCoord.out.bam.bai'
    conda: 'envs/samtools.yml'
    shell:'''
    samtools index {input}
    '''
    
rule htseq_count:
    input:
        bam = 'outputs_5mil/star/{sample}Aligned.sortedByCoord.out.bam',
        gtf = 'inputs/genome/GCF_000146045.2_R64_genomic.gtf'
    output: "outputs_5mil/htseq/{sample}_readcounts.txt"
    conda: "envs/htseq.yml"
    shell:'''
    htseq-count -m intersection-nonempty -s yes -f bam -r pos {input.bam} {input.gtf} > {output}
    '''

rule make_counts:
    input: expand("outputs_5mil/htseq/{sample}_readcounts.txt", sample = SAMPLES)
    output: "outputs_5mil/counts/raw_counts.tsv"
    conda: "envs/tidyverse.yml"
    script: "scripts/make_raw_counts.R"
