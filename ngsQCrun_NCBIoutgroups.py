import numpy as np

ref_path = '/space/s2/joana/refgenomes/V.lagopus/vulpes_lagopus_scilife_genome.fasta'

accessions = {
      #  'arctic_fox' : ['SRR7704840'],
       # 'red_fox_1': ['SRR5328115'],
        #'red_fox_2' : ['SRR5328114'],
        #'red_fox_3' : ['SRR5328113'],
#        'TAME_01' : ['SRR5280510'],
#        'TAME_02' : ['SRR5280509'],
#        'TAME_03' : ['SRR5280508'],
#        'TAME_04' : ['SRR5280507'],
#        'TAME_05' : ['SRR5280506'],
#        'TAME_06' : ['SRR5280505'],
#        'TAME_07' : ['SRR5280504'],
#        'TAME_08' : ['SRR5280503'],
#        'TAME_09' : ['SRR5280502'],
#        'TAME_10' : ['SRR5280501'],
#        'Agressive_1' : ['SRR5280495'],
#        'Agressive_2' : ['SRR5280494'],
#        'Agressive_3' : ['SRR5280493'],
#        'Agressive_4' : ['SRR5280492'],
#        'Agressive_5' : ['SRR5280491'],
#        'Agressive_6' : ['SRR5280490'],
#        'Agressive_7' : ['SRR5280489'],
#        'Agressive_8' : ['SRR5280488'],
#        'Agressive_9' : ['SRR5280487'],
#        'Agressive_10' : ['SRR5280486'],
#        'Conventional_1' : ['SRR5280485'],
#        'Conventional_2' : ['SRR5280484'],
#        'Conventional_3' : ['SRR5280483'],
#        'Conventional_4' : ['SRR5280482'],
#        'Conventional_5' : ['SRR5280481'],
#        'Conventional_6' : ['SRR5280480'],
#        'Conventional_7' : ['SRR5280479'],
#        'Conventional_8' : ['SRR5280478'],
#        'Conventional_9' : ['SRR5280477'],
#        'Conventional_10' : ['SRR5280476'],
#        'HIC_redfox' : ['SRR8616952', 'SRR8616953']
        #'dog_tibetan_mastiff': ['SRR1138360'],
        #'dog_german_shepherd': ['SRR1122359'],
        #'dog_basenji' : ['SRR2149861'],
        #'coyote_alska' : ['SRR8049186'],
        #'coyote_mexico' : ['SRR8049187'],
        #'grey_wolf_mexico' : ['SRR8049200'],
        #'grey_wolf_saudi_arabia' : ['SRR8049193'],
        #'grey_wolf_daneborg' : ['SRR8049195'],
        #'grey_wolf_ellesmere' : ['SRR8049197'],
        #'grey_wolf_israel' : ['SRR2149870'],
        #'grey_wolf_syria' : ['SRR8049194'],
        #'dingo' : ['SRR2149867'],
        #'african_golden_wolf' : ['SRR8049196'],
        #'golden_jackal_syria' : ['SRR8049192'],
        #'golden_jackal_israel' : ['SRR2149878'],
        #'ethiopian_wolf' : ['SRR8049190'],
        #'dhole' : ['SRR8049189'],
        #'african_wild_dog' : ['SRR7874817'],
        #'black_backed_jackal' : ['ERR3210523'],
        #'culpeo_fox' : ['SRR1066703'],
        #'gray_fox' : ['SRR5198019', 'SRR5198020', 'SRR5198021'],
        #'island_fox' : ['SRR5198018', 'SRR5198017', 'SRR5198016'],
        #'polar_bear' : ['SRR942310']
        #'tibetan_wolf_1': ['SRR2017901'],
        #'tibetan_wolf_2': ['SRR2017902'],
        #'tibetan_wolf_3' : ['SRR2017892'],
        #'himalayan_wolf': ['SRR11085387']
        }

samples = accessions.keys()


rule all:
    input: expand('bam_final/{sample}.sorted.bam', sample=samples)
    # expand('fastqs/{sample}_{number}.fastq.gz', sample=samples, number=[1,2])


rule prefetch_samples:
    output: temp('accessions/{accession}.sra')
    shell: """
    prefetch -X 200G -o {output} {wildcards.accession}
    """

rule fastq_dump:
    input: 'accessions/{accession}.sra'
    output:
        temp('accessions/{accession}_1.fastq.gz'),
        temp('accessions/{accession}_2.fastq.gz')
    shell: """
    cd accessions && fastq-dump --gzip --split-files {wildcards.accession}.sra
    """

def get_input_merge_fastq(wildcards):
    return expand('accessions/{{accession}}_{number}.fastq.gz'.format(number=wildcards.number),
            accession=accessions[wildcards.sample])

rule merge_fastq:
    input: get_input_merge_fastq
    output: temp('fastqs/{sample}_{number}.fastq.gz')
    shell: 'cat {input} > {output}'


rule cut_adapt:
    input:
        'fastqs/{sample}_1.fastq.gz',
        'fastqs/{sample}_2.fastq.gz'
    output:
        temp('fastqs/{sample}_1.trim.fastq.gz'),
        temp('fastqs/{sample}_2.trim.fastq.gz')
    shell: 'cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -g GCTCTTCCGATCT -G GCTCTTCCGATCT -n 8 -e 0.1 -O 1 -m 30 -q 20,20 --max-n 0.5 --pair-filter any -o {output[0]} -p {output[1]} {input[0]} {input[1]}'


rule trim_clean:
    input:
        'fastqs/{sample}_1.trim.fastq.gz',
        'fastqs/{sample}_2.trim.fastq.gz'
    output:
        temp('fastqs/{sample}_1.trim1P.fastq.gz'),
        temp('fastqs/{sample}_1.trim1U.fastq.gz'),
        temp('fastqs/{sample}_2.trim2P.fastq.gz'),
        temp('fastqs/{sample}_2.trim2U.fastq.g'),
    shell:'java -jar /usr/local/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -quiet -threads 8 -phred33 {input[0]} {input[1]} {output[0]} {output[1]} {output[2]} {output[3]} SLIDINGWINDOW:4:20 MINLEN:30'

rule map_filt:
    input:
        'fastqs/{sample}_1.trim1P.fastq.gz',
        'fastqs/{sample}_2.trim2P.fastq.gz'
    output: temp('bams/{sample}.flagfilt.bam')
    shell: 'bwa mem -t 2 {ref_path} {input[0]} {input[1]} | samtools view -q 15 -bT {ref_path} -F 780 -o {output}'

rule map_filt_sort:
    input: 'bams/{sample}.flagfilt.bam'
    output: temp('bams/{sample}.flagfilt.sorted.bam')
    shell: 'samtools sort -o {output} {input}'

rule picard_rmdup:
    input: 'bams/{sample}.flagfilt.sorted.bam'
    output: temp('bams/{sample}.flagfilt.rmdup.bam')
    shell: """
    mkdir bams/tmp_0_{wildcards.sample};
    java -Xmx31G -XX:ParallelGCThreads=2 -Djava.io.tmpdir=bams/tmp_0_{wildcards.sample} -jar /home/joana/software/ngsQC/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I={input}  O={output} M=bams/{wildcards.sample}.rmdup.metrics
    """

rule rmdup_sort:
    input: 'bams/{sample}.flagfilt.rmdup.bam'
    output: temp('bams/{sample}.flagfilt.rmdup.sorted.bam')
    shell: 'samtools sort -o {output} {input}'

rule add_readgroup:
    input: 'bams/{sample}.flagfilt.rmdup.sorted.bam'
    output: temp('bams/{sample}.flagfilt.rmdup.sorted.rg.bam')
    shell: 'samtools addreplacerg -r ID:S1 -r SM:{wildcards.sample} -r PU:L1 -r PL:ILLUMINA -o {output} {input}'



rule indel_remap:
    input: 'bams/{sample}.flagfilt.rmdup.sorted.rg.bam'
    output:
        temp('bam_final/{sample}.flagfilt.rmdup.sorted.rg.intervals'),
        temp('bam_final/{sample}.flagfilt.rmdup.sorted.realign.rg.bam')
    shell: """
    samtools index {input} &&
    mkdir bam_final/tmp_0_{wildcards.sample} &&
    /home/joana/software/jre1.8.0_211/bin/java -Xmx5G -XX:ParallelGCThreads=8 -Djava.io.tmpdir=bam_final/tmp_0_{wildcards.sample} -jar /home/joana/software/GenomeAnalysisTK-3.5-0-g36282e4/GenomeAnalysisTK.jar -T RealignerTargetCreator -R {ref_path} -I {input} -o {output[0]} &&
    /home/joana/software/jre1.8.0_211/bin/java -Xmx4g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=bam_final/tmp_0_{wildcards.sample} -jar /home/joana/software/GenomeAnalysisTK-3.5-0-g36282e4/GenomeAnalysisTK.jar -T IndelRealigner -I {input} -R {ref_path}  -targetIntervals {output[0]} -o {output[1]}
    """


#  Create final files for snpCleaner
rule bam_final:
    input: 'bam_final/{sample}.flagfilt.rmdup.sorted.realign.rg.bam'
    output: 'bam_final/{sample}.sorted.bam'
    shell: """
    samtools sort -o {output} {input} &&
    samtools index {output}
    """

