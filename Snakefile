SAMPLES=['AP_A_replicat1', 'AP_A_replicat2', 'AP_A_replicat3',
         'AP_B_replicat1', 'AP_B_replicat2', 'AP_B_replicat3',
         'AP_C_replicat1', 'AP_C_replicat2', 'AP_C_replicat3',
         'AP_D_replicat1', 'AP_D_replicat2', 'AP_D_replicat3',
         'AP_E_replicat1', 'AP_E_replicat2', 'AP_E_replicat3',
         'AP_F_replicat1', 'AP_F_replicat2', 'AP_F_replicat3',
         'AP_G_replicat1', 'AP_G_replicat2', 'AP_G_replicat3']

rule all:
    input: 
        "outputs/fastp_trimmed/multiqc_report.html",
        expand('outputs/star/{sample}Aligned.sortedByCoord.out.bam.bai', sample = SAMPLES),
        "outputs/counts/raw_counts.tsv",
        "outputs/fastp_ribo/multiqc_report.html",
        "outputs/star/multiqc_report.html"

rule separate_R1:
    input: "{sample}.trimmed.fq.gz"
    output: "outputs/separated/{sample}_R1.trimmed.fq"
    shell:'''
    zcat {input} | grep -A 3 "1:N:0:" > {output}    
    '''

rule separate_R2:
    input: "{sample}.trimmed.fq.gz"
    output: "outputs/separated/{sample}_R2.trimmed.fq"
    shell:'''
    zcat {input} | grep -A 3 "2:N:0:" > {output}    
    '''

rule repair_trimmed_reads:
    input:
        R1="outputs/separated/{sample}_R1.trimmed.fq",
        R2="outputs/separated/{sample}_R2.trimmed.fq"
    output:
        R1="outputs/repair/{sample}_R1.trimmed.fq",
        R2="outputs/repair/{sample}_R2.trimmed.fq",
        singletons="outputs/repair/{sample}_singleton.fq"
    conda:"envs/bbmap.yml"
    shell:'''
    repair.sh in1={input.R1} in2={input.R2} out1={output.R1} out2={output.R2} outs={output.singletons} repair
    '''

rule fastp_trimmed_reads:
    input:
        R1= "outputs/repair/{sample}_R1.trimmed.fq",
        R2= "outputs/repair/{sample}_R2.trimmed.fq"
    output: 
        json="outputs/fastp_trimmed/{sample}.trimmed.fastp.json",
        html="outputs/fastp_trimmed/{sample}.trimmed.fastp.html",
    conda: "envs/fastp.yml"
    shell:'''
    fastp -i {input.R1} -I {input.R2} -j {output.json} -h {output.html}
    '''

rule multiqc_fastp_trimmed:
    input: expand("outputs/fastp_trimmed/{sample}.trimmed.fastp.json", sample=SAMPLES)
    output: "outputs/fastp_trimmed/multiqc_report.html"
    params: 
        indir = "outputs/fastp_trimmed",
        outdir = "outputs/fastp_trimmed"
    conda: "envs/multiqc.yml"
    shell:'''
    multiqc {params.indir} -o {params.outdir} 
    '''

#################################
## check ribo
#################################

rule bbduk_find_ribo:
    output:
        ribo_R1='outputs/ribo/{sample}_ribo_R1.trimmed.fq',
        ribo_R2='outputs/ribo/{sample}_ribo_R2.trimmed.fq',
    input: 
        R1= "outputs/repair/{sample}_R1.trimmed.fq",
        R2= "outputs/repair/{sample}_R2.trimmed.fq",
        ribo='inputs/ribokmers.fa.gz'
    conda: 'envs/bbmap.yml'
    shell:'''
    bbduk.sh -Xmx4g in={input.R1} in2={input.R2} \
        outm={output.ribo_R1} outm2={output.ribo_R2} \ 
        k=31 ref={input.ribo}
    '''

rule fastp_ribo_reads:
    input: 
        R1='outputs/ribo/{sample}_ribo_R1.trimmed.fq',
        R2='outputs/ribo/{sample}_ribo_R2.trimmed.fq',
    output: 
        json="outputs/fastp_ribo/{sample}_ribo.fastp.json",
        html="outputs/fastp_ribo/{sample}_ribo.fastp.html",
    conda: "envs/fastp.yml"
    shell:'''
    fastp -i {input.R1} -I {input.R2} -j {output.json} -h {output.html}
    '''

rule multiqc_fastp_ribo:
    input: expand("outputs/fastp_ribo/{sample}_ribo.fastp.json", sample=SAMPLES)
    output: "outputs/fastp_ribo/multiqc_report.html"
    params: 
        indir = "outputs/fastp_ribo",
        outdir = "outputs/fastp_ribo"
    conda: "envs/multiqc.yml"
    shell:'''
    multiqc {params.indir} -o {params.outdir} 
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

rule star_align:
    #v2.5.2a
    input:
        R1= "outputs/repair/{sample}_R1.trimmed.fq",
        R2= "outputs/repair/{sample}_R2.trimmed.fq",
        genome_index = 'inputs/genome/SAindex'
    output: 'outputs/star/{sample}Aligned.sortedByCoord.out.bam' 
    params: 
        out_prefix = lambda wildcards: 'outputs/star/' + wildcards.sample,
        genome_dir = 'inputs/genome'
    conda: 'envs/star.yml'
    threads: 2
    shell:'''
    STAR --runThreadN {threads} --genomeDir {params.genome_dir}      \
        --readFilesIn {input.R1} {input.R2} --outSAMtype BAM \
        SortedByCoordinate --outFileNamePrefix {params.out_prefix}
    '''

rule index_bam:
    input: 'outputs/star/{sample}Aligned.sortedByCoord.out.bam'
    output: 'outputs/star/{sample}Aligned.sortedByCoord.out.bam.bai'
    conda: 'envs/samtools.yml'
    shell:'''
    samtools index {input}
    '''

rule flagstat_bam:
    input: 'outputs/star/{sample}Aligned.sortedByCoord.out.bam'
    output: 'outputs/star/{sample}Aligned.sortedByCoord.out.bam.flagstat'
    conda: 'envs/samtools.yml'
    shell:'''
    samtools flagstat {input} > {output}
    '''

rule multiqc_flagstat:
    input: expand('outputs/star/{sample}Aligned.sortedByCoord.out.bam.flagstat', sample=SAMPLES)
    output: "outputs/star/multiqc_report.html"
    params: 
        indir = "outputs/star",
        outdir = "outputs/star"
    conda: "envs/multiqc.yml"
    shell:'''
    multiqc {params.indir} -o {params.outdir} 
    '''

rule htseq_count:
    input:
        bam = 'outputs/star/{sample}Aligned.sortedByCoord.out.bam',
        gtf = 'inputs/genome/GCF_000146045.2_R64_genomic.gtf'
    output: "outputs/htseq/{sample}_readcounts.txt"
    conda: "envs/htseq.yml"
    shell:'''
    htseq-count -m intersection-nonempty -s yes -f bam -r pos --stranded=reverse {input.bam} {input.gtf} > {output}
    '''

rule make_counts:
    input: expand("outputs/htseq/{sample}_readcounts.txt", sample = SAMPLES)
    output: "outputs/counts/raw_counts.tsv"
    conda: "envs/tidyverse.yml"
    script: "scripts/make_raw_counts.R"
