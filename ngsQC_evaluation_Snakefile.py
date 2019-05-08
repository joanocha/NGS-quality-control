import numpy as np
# to get function np.unique (removes repetitions in a list)
# for example np.unique(['sample1', 'sample2', 'sample2', 'sample3']) -> ['sample1', 'sample2', 'sample3']


READS_FILENAMES = os.listdir('fastq_raw')
READS = [('_').join(f.split('.')[0].split('_')[:5]) for f in READS_FILENAMES]
SAMPLE_RUNS = np.unique([('_').join(f.split('_')[:3]) for f in READS_FILENAMES])
SAMPLES = np.unique([sr.split('_')[0] for sr in SAMPLE_RUNS])


ref_path = '/space/s1/joana/refgenomes/V.lagopus_10X_Stefan/vulpes_lagopus_scilife_genome.fasta'

rule all:
    input:
        'evaluation/fastqc_raw/all_mod_scores.csv',
        'evaluation/fastqc_clean/all_mod_scores.csv',
        'evaluation/fastq_coverage.txt',
        'evaluation/before_merge_coverage.txt',
        'evaluation/after_merge_coverage.txt'


### SUMMARIZE STATS 
# fastq_raw
# fastq_clean # has trimmed reads
# bam_filt  # has mapped, filtered for mapping quality reads and sorted reads
# bam_rmdup # has duplicates removed and sorted reads
# bam_merge #has merged libraries into samples
# bam_indelremap #has indel realigned reads

# For fastq files raw and trimmed

#for filename in *.fastq.gz; do  mv "$filename" "${filename//-LP-2/}"; done; #to get rid of LP-2 and simplify names
#for filename in *.fastq.gz; do  mv "$filename" "${filename//-LP2/}"; done; #to get rid of LP2 and simplify names


rule run_fastqc_raw:
    input: expand('fastq_raw/{read}.fastq.gz', read=READS)
    output:
        expand('evaluation/fastqc_raw/{read}_fastqc.html', read=READS),
        expand('evaluation/fastqc_raw/{read}_fastqc.zip', read=READS),
    shell: 'fastqc -o evaluation/fastqc_raw -f fastq -t 24 fastq_raw/*.fastq.gz'

rule summarize_fastqc_raw:
    input: expand('evaluation/fastqc_raw/{read}_fastqc.zip', read=READS)
    output: 'evaluation/fastqc_raw/all_mod_scores.csv'
    shell: 'cd evaluation/fastqc_raw && python /home/joana/scripts/summarize_fastqc.py'

rule run_fastqc_clean:
    input:
        expand('fastq_clean/{sample_run}.trim_1P.fastq.gz', sample_run=SAMPLE_RUNS),
        expand('fastq_clean/{sample_run}.trim_2P.fastq.gz', sample_run=SAMPLE_RUNS)
    output:
        expand('evaluation/fastqc_clean/{sample_run}.trim_1P_fastqc.html', sample_run=SAMPLE_RUNS),
        expand('evaluation/fastqc_clean/{sample_run}.trim_2P_fastqc.html', sample_run=SAMPLE_RUNS),
        expand('evaluation/fastqc_clean/{sample_run}.trim_1P_fastqc.zip', sample_run=SAMPLE_RUNS),
        expand('evaluation/fastqc_clean/{sample_run}.trim_2P_fastqc.zip', sample_run=SAMPLE_RUNS),
    shell: 'fastqc -o evaluation/fastqc_clean -f fastq -t 24 fastq_clean/*.fastq.gz'

rule summarize_fastqc_clean:
    input:
        expand('evaluation/fastqc_clean/{sample_run}.trim_1P_fastqc.zip', sample_run=SAMPLE_RUNS),
        expand('evaluation/fastqc_clean/{sample_run}.trim_2P_fastqc.zip', sample_run=SAMPLE_RUNS),
    output: 'evaluation/fastqc_clean/all_mod_scores.csv'
    shell: 'cd evaluation/fastqc_clean && python /home/joana/scripts/summarize_fastqc.py'

rule nbases_fastq_raw:
    input: expand('fastq_raw/{read}.fastq.gz', read=READS)
    output: 'evaluation/nbases_raw.txt'
    shell: """
    for i in fastq_raw/*.fastq.gz
    do
    printf $i"\t" >> {output}
    zcat $i|paste - - - - | cut -f 2 | tr -d '\n' | wc -c >> {output}
    done
    """

rule nbases_fastq_clean:
    input:
        expand('fastq_clean/{sample_run}.trim_1P.fastq.gz', sample_run=SAMPLE_RUNS),
        expand('fastq_clean/{sample_run}.trim_2P.fastq.gz', sample_run=SAMPLE_RUNS)
    output: 'evaluation/nbases_clean.txt'
    shell: """
    for i in fastq_clean/*.fastq.gz
    do
    printf $i"\t" >> {output}
    zcat $i|paste - - - - | cut -f 2 | tr -d '\n' | wc -c >> {output}
    done
    """

rule cov_fastq:
    input:
        'evaluation/nbases_raw.txt',
        'evaluation/nbases_clean.txt'
    output:
        temp('evaluation/cov_raw.txt'),
        temp('evaluation/cov_clean.txt'),
        'evaluation/fastq_coverage.txt'
    shell:"""
    REF=/space/s1/joana/refgenomes/V.lagopus_10X_Stefan/vulpes_lagopus_scilife_genome.bed
    GENLEN=$(awk '{{genlen+=$3}} END {{print genlen}}' $REF)
    awk -v tot="$GENLEN" '{{print $0"\t"$2/tot}}' {input[0]} > {output[0]}
    awk -v tot="$GENLEN" '{{print $0"\t"$2/tot}}' {input[1]} > {output[1]}
    join -j1 {output[0]} {output[1]} | cut -f 1,3,5 -d ' ' | sed 's/ /\t/g' > {output[2]}
    """

# from bam onwards
rule nbases_cov_bam_filt:
    input:
        expand('bam_filt/{sample_run}.flagfilt.sorted.bam', sample_run=SAMPLE_RUNS)
    output:
        'evaluation/avgcov_bamfilt.txt',
        'evaluation/nbases_bamfilt.txt'
    shell: """
    REF=/space/s1/joana/refgenomes/V.lagopus_10X_Stefan/vulpes_lagopus_scilife_genome.bed
    ALLCOV={output[0]} #output file for coverage
    ALLNBASES={output[1]} #output file for nbases
    for i in bam_filt/*sorted.bam
    do
    samtools index $i
    samtools bedcov $REF $i >> $i.cov.txt
    printf "$i\t" >> $ALLCOV
    awk '{{total+=$3; nbases+=$4}} END {{print nbases/total}}' $i.cov.txt >> $ALLCOV
    printf "$i\t" >> $ALLNBASES
    awk '{{nbases+=$4}} END {{print nbases}}' $i.cov.txt >> $ALLNBASES
    done
    """

rule nbases_cov_bam_rmdup:
    input:
        expand('bam_rmdup/{sample_run}.flagfilt.rmdup.sorted.rg.bam', sample_run=SAMPLE_RUNS)
    output:
        'evaluation/avgcov_bamrmdup.txt',
        'evaluation/nbases_bamrmdup.txt'
    shell: """
    REF=/space/s1/joana/refgenomes/V.lagopus_10X_Stefan/vulpes_lagopus_scilife_genome.bed
    ALLCOV={output[0]} #output file for coverage
    ALLNBASES={output[1]} #output file for nbases
    for i in bam_rmdup/*rg.bam
    do
    samtools index $i
    samtools bedcov $REF $i >> $i.cov.txt
    printf "$i\t" >> $ALLCOV
    awk '{{total+=$3; nbases+=$4}} END {{print nbases/total}}' $i.cov.txt >> $ALLCOV #
    printf "$i\t" >> $ALLNBASES
    awk '{{nbases+=$4}} END {{print nbases}}' $i.cov.txt >> $ALLNBASES
    done
    """


rule summarize_steps_before_merge:
    input:
        'evaluation/avgcov_bamfilt.txt',
        'evaluation/avgcov_bamrmdup.txt'
    output:
        'evaluation/before_merge_coverage.txt'
    shell: """
    printf 'sample\t' > {output}
    ls evaluation/avgcov* | paste - - - - - - >> {output}
    paste evaluation/avgcov* | cut -f 1,2,4,6,8,10,12 | tail -n +2>> {output}
    """

rule nbases_cov_bam_merge:
    input:
        expand('bam_merge/{sample}.flagfilt.sorted.rmdup.merged.bam', sample=SAMPLES)
    output:
        'evaluation/after_merge_avgcov_bammerge.txt',
        'evaluation/after_merge_nbases_bammerge.txt'
    shell: """
    REF=/space/s1/joana/refgenomes/V.lagopus_10X_Stefan/vulpes_lagopus_scilife_genome.bed
    ALLCOV={output[0]} #output file for coverage
    ALLNBASES={output[1]} #output file for nbases
    for i in bam_merge/*merged.bam
    do
    samtools index $i
    samtools bedcov $REF $i >> $i.cov.txt
    printf "$i\t" >> $ALLCOV
    awk '{{total+=$3; nbases+=$4}} END {{print nbases/total}}' $i.cov.txt >> $ALLCOV
    printf "$i\t" >> $ALLNBASES
    awk '{{nbases+=$4}} END {{print nbases}}' $i.cov.txt >> $ALLNBASES
    done
    """

rule nbases_cov_bam_indelremap:
    input:
        expand('bam_indelremap/{sample}.flagfilt.rmdup.sorted.rg.merged.realign.bam', sample=SAMPLES)
    output:
        'evaluation/after_merge_avgcov_bamindelremap.txt',
        'evaluation/after_merge_nbases_bamindelremap.txt'
    shell: """
    REF=/space/s1/joana/refgenomes/V.lagopus_10X_Stefan/vulpes_lagopus_scilife_genome.bed
    ALLCOV={output[0]}
    ALLNBASES={output[1]}
    for i in bam_indelremap/*realign.bam
    do
    samtools index $i
    samtools bedcov $REF $i >> $i.cov.txt
    printf "$i\t" >> $ALLCOV
    awk '{{total+=$3; nbases+=$4}} END {{print nbases/total}}' $i.cov.txt >> $ALLCOV
    printf "$i\t" >> $ALLNBASES
    awk '{{nbases+=$4}} END {{print nbases}}' $i.cov.txt >> $ALLNBASES
    done
    """


rule summarize_steps_after_merge:
    input:
        'evaluation/after_merge_avgcov_bammerge.txt',
        'evaluation/after_merge_avgcov_bamindelremap.txt'
    output:
        'evaluation/after_merge_coverage.txt'
    shell: """
    printf 'sample\t' > {output}
    ls evaluation/after_merge_avgcov* | paste - - - - - - >> {output}
    paste evaluation/after_merge_avgcov* | cut -f 1,2,4,6,8,10,12 | tail -n +2>> {output}
    """
