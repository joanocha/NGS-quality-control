import numpy as np
# to get function np.unique (removes repetitions in a list)
# for example np.unique(['sample1', 'sample2', 'sample2', 'sample3']) -> ['sample1', 'sample2', 'sample3']


READS = os.listdir('fastq_raw')
SAMPLE_RUNS = np.unique([('_').join(f.split('_')[:3]) for f in READS])
SAMPLES = np.unique([sr.split('_')[0] for sr in SAMPLE_RUNS])


ref_path = '/space/s1/joana/refgenomes/V.lagopus_10X_Stefan/vulpes_lagopus_scilife_genome.fasta'

rule all:
     input: expand('bam_indelremap/{sample}.flagfilt.rmdup.sorted.rg.merged.realign.bam', sample=SAMPLES)
 

# 1. TRIM READS
rule trim_reads:
    input: 'fastq_raw/{sample_run}_R1_001.fastq.gz'
    output:
        'fastq_clean/{sample_run}.trim_1P.fastq.gz',
        'fastq_clean/{sample_run}.trim_2P.fastq.gz'
    shell: """
    java -jar /usr/local/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 20 -basein {input} -baseout fastq_clean/{wildcards.sample_run}.trim.fastq.gz -trimlog fastq_clean/{wildcards.sample_run}.trim.log ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30
    """

rule run_fastqc_clean: # this was recently added and need to run it by hand after run finnishes (April 18th 2019)
    shell: 'fastqc -o evaluation/fastqc_clean -f fastq -t 24 fastq_clean/*.fastq.gz'

rule summarize_fastqc_clean:
    shell: 'cd evaluation/fastqc_clean && python /home/joana/scripts/summarize_fastqc.py' # this was recently added and need to run it by hand after run finnishes (April 18th 2019)



# 2. MAP READS WHILE FILTERING
rule map_filt:
    input:
        'fastq_clean/{sample_run}.trim_1P.fastq.gz',
        'fastq_clean/{sample_run}.trim_2P.fastq.gz'
    output:
        'bam_filt/{sample_run}.flagfilt.bam'  # temp means that file will be removed once all the rules that need it as input finish
    shell: """
    bwa mem -t 24 {ref_path} {input[0]} {input[1]} | samtools view -q 15 -bT {ref_path} -F 780 -o {output}
    """

rule map_filt_sort: #will have to remove unsorted files by hand after sorting as the temp file was added after the run started (April 18th 2019)
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
# java -jar /home/joana/software/ngsQC/picard.jar CreateSequenceDictionary R=vulpes_lagopus_scilife_genome.fasta O=vulpes_lagopus_scilife_genome.dict
# this command is necessary for gatk to work. this should be done with any ref genome of interest in the respective directory

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
