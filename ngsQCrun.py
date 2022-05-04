import numpy as np
# to get function np.unique (removes repetitions in a list)
# for example np.unique(['sample1', 'sample2', 'sample2', 'sample3']) -> ['sample1', 'sample2', 'sample3']

 # The following workflow applies to Paired-end Illumina NovaSeq data.
 # Please ensure you have the following working directories for script to work as it is. Otherwise just change it to incorporate your own folder name
 # ### fastq_raw -> folder containing all your fastq.gz folders
 # ### fastq_clean-> directory that will contain your trimmed folders
 # ### bam_fil -> where mapped, filtered and sorted files will be
 # ### bam_rmdup -> duplicate removed bam files
 # ### bam_indel_remap -> GATK output after running RealignTargetCreater and IndelRealigner
 # ### bam_map -> mapdamage outputs

READS = os.listdir('fastq_raw')
SAMPLE_RUNS = np.unique([('_').join(f.split('_')[:3]) for f in READS])
SAMPLES = np.unique([sr.split('_')[0] for sr in SAMPLE_RUNS])


ref_path = '/space/s1/joana/refgenomes/V.lagopus/vulpes_lagopus_scilife_genome.fasta'

#please make sure you ran the following comands on ref_path:

#bwa index /space/s1/joana/refgenomes/V.lagopus/vulpes_lagopus_scilife_genome.fasta
#samtools faidx /space/s1/joana/refgenomes/V.lagopus/vulpes_lagopus_scilife_genome.fasta
#awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' /space/s1/joana/refgenomes/V.lagopus/vulpes_lagopus_scilife_genome.fasta.fai > /space/s1/joana/refgenomes/V.lagopus/vulpes_lagopus_scilife_genome.bed
#java -jar /home/joana/software/ngsQC/picard.jar CreateSequenceDictionary R=/space/s1/joana/refgenomes/V.lagopus/vulpes_lagopus_scilife_genome.fastaa O=/space/s1/joana/refgenomes/V.lagopus/vulpes_lagopus_scilife_genome.dict


rule all:
    input:
        expand('bam_final/{sample}.sorted.bam', sample=SAMPLES)

# 1. TRIM READS

rule cut_adapt:
    input:
        'fastq_raw/{sample_run}_R1_001.fastq.gz',
        'fastq_raw/{sample_run}_R2_001.fastq.gz'
    output:
        temp('fastq_clean/{sample_run}_trim_R1_001.fastq.gz'),
        temp('fastq_clean/{sample_run}_trim_R2_001.fastq.gz')
    shell: 'cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -g GCTCTTCCGATCT -G GCTCTTCCGATCT -n 2 -e 0.1 -O 1 -m 30 -q 20,20 --max-n 0.5 --nextseq-trim=20 --pair-filter any -o {output[0]} -p {output[1]} {input[0]} {input[1]}'

# use --nextseq-trim only if you have NovaSeq data

rule trim_clean:
    input:
        'fastq_clean/{sample_run}_trim_R1_001.fastq.gz',
        'fastq_clean/{sample_run}_trim_R2_001.fastq.gz'
    output:
        'fastq_clean/{sample_run}.trim_1P.fastq.gz',
        'fastq_clean/{sample_run}.trim_1U.fastq.gz',
        'fastq_clean/{sample_run}.trim_2P.fastq.gz',
        'fastq_clean/{sample_run}.trim_2U.fastq.gz'
    shell:'java -jar /usr/local/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -quiet -threads 8 -phred33 {input[0]} {input[1]} {output[0]} {output[1]} {output[2]} {output[3]} SLIDINGWINDOW:4:20 MINLEN:30'

# 2. MAP READS WHILE FILTERING
rule map_filt:
    input:
        'fastq_clean/{sample_run}.trim_1P.fastq.gz',
        'fastq_clean/{sample_run}.trim_2P.fastq.gz'
    output: temp('bam_filt/{sample_run}.flagfilt.bam')
    shell: 'bwa mem -t 2 {ref_path} {input[0]} {input[1]} | samtools view -q 15 -bT {ref_path} -F 780 -o {output}'

rule map_filt_sort:
    input: 'bam_filt/{sample_run}.flagfilt.bam'
    output: 'bam_filt/{sample_run}.flagfilt.sorted.bam'
    shell: 'samtools sort -o {output} {input}'

# 3. MARK DUPLICATES
rule picard_rmdup:
    input: 'bam_filt/{sample_run}.flagfilt.sorted.bam'
    output: temp('bam_rmdup/{sample_run}.flagfilt.rmdup.bam')
    shell: """
    mkdir bam_rmdup/tmp_0_{wildcards.sample_run};
    java -Xmx31G -XX:ParallelGCThreads=2 -Djava.io.tmpdir=bam_rmdup/tmp_0_{wildcards.sample_run} -jar /home/joana/software/ngsQC/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I={input}  O={output} M=evaluation/{wildcards.sample_run}.rmdup.metrics
    """

rule rmdup_sort:
    input: 'bam_rmdup/{sample_run}.flagfilt.rmdup.bam'
    output: temp('bam_rmdup/{sample_run}.flagfilt.rmdup.sorted.bam')
    shell: 'samtools sort -o {output} {input}'


# 4. ADD READGROUPS AND MERGE LIBRARIES FROM DIFFERENT RUNS INTO SAMPLES
rule add_readgroup:
    input: 'bam_rmdup/{sample}_{read_group}_{lane}.flagfilt.rmdup.sorted.bam'
    output: 'bam_rmdup/{sample}_{read_group}_{lane}.flagfilt.rmdup.sorted.rg.bam'
    shell:
        'samtools addreplacerg -r ID:{wildcards.read_group}_{wildcards.lane} -r SM:{wildcards.sample} -r PU:{wildcards.lane} -r PL:ILLUMINA -o {output} {input}'

rule merge_to_samples:
    input: expand('bam_rmdup/{sample_run}.flagfilt.rmdup.sorted.rg.bam', sample_run=SAMPLE_RUNS)
    output: expand('bam_merge/{sample}.flagfilt.rmdup.sorted.rg.merged.bam', sample=SAMPLES)
    run:
        for sample in SAMPLES:
            input = 'bam_rmdup/{}_*.rg.bam'.format(sample)
            output = 'bam_merge/{}.flagfilt.rmdup.sorted.rg.merged.bam'.format(sample)
            cmd = 'samtools merge -f {} {}'.format(output, input)
            print(cmd)
            shell(cmd)


# 5. INDEL REAMMAPPING
# java -jar /home/joana/software/ngsQC/picard.jar CreateSequenceDictionary R=/space/s1/joana/refgenomes/V.lagopus/vulpes_lagopus_scilife_genome.fastaa O=/space/s1/joana/refgenomes/V.lagopus/vulpes_lagopus_scilife_genome.dict  is necessary for gatk to work

rule indel_remap:
    input: 'bam_merge/{sample}.flagfilt.rmdup.sorted.rg.merged.bam'
    output:
        'bam_indelremap/{sample}.flagfilt.rmdup.sorted.rg.merged.intervals',
        'bam_indelremap/{sample}.flagfilt.rmdup.sorted.rg.merged.realign.bam'
    shell:"""
    samtools index {input} &&
    mkdir bam_indelremap/tmp_0_{wildcards.sample} &&
    /home/joana/software/jre1.8.0_211/bin/java -Xmx5G -XX:ParallelGCThreads=8 -Djava.io.tmpdir=bam_indelremap/tmp_0_{wildcards.sample} -jar /home/joana/software/GenomeAnalysisTK-3.5-0-g36282e4/GenomeAnalysisTK.jar -T RealignerTargetCreator -R {ref_path} -nt 16 -I {input} -o {output[0]} &&
    /home/joana/software/jre1.8.0_211/bin/java -Xmx4g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=bam_indelremap/tmp_0_{wildcards.sample} -jar /home/joana/software/GenomeAnalysisTK-3.5-0-g36282e4/GenomeAnalysisTK.jar -T IndelRealigner -I {input} -R {ref_path} -targetIntervals {output[0]} -o {output[1]}
    """

#6. MAP DAMAGE

rule map_damage:
    input: 'bam_indelremap/{sample}.flagfilt.rmdup.sorted.rg.merged.realign.bam'
    output: 'bam_damage/results_{sample}'
    shell: """
    samtools index {input};
    /home/joana/.local/bin/mapDamage -d {output} plot  --merge-reference-sequences -n 0.5 -v -i {input} -r {ref_path}
    """

# 7. Create final files for snpCleaner
#the input either cames from bam_indelremap or bam_damage, depending on damage patterns
rule bam_final:
    input: 'bam_indelremap/{sample}.flagfilt.rmdup.sorted.rg.merged.realign.bam'
    output: 'bam_final/{sample}.sorted.bam'
    shell: """
    samtools sort -o {output} {input}' &&
    samtools index {output}
    """

