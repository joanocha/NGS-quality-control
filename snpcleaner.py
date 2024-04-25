import numpy as np
scaffold_ids = np.arange(0, 4048)

### First task was to create a bam file list as this will be a required input in ANGSD
# ls bam_final/ | grep ".bam$" > tmp1.txt
# sed 's/^/\/space\/s1\/joana\/data\/bam_final\//' tmp1.txt > all_vulpes_bams.txt
# the same was performed for all other lists of bam
# bam list
# ls bam_final/ | grep ".bam$" > tmp1.txt
# sed 's/^/\/global\/scratch\/users\/joana_rocha\/PANPAN\/CHIMPS_lowcov\/bam_final\//' tmp1.txt > allchimps_bams.txt
 

 
 
rule all:
    input:
        expand('sites/all_vulpes_scaffold_{scaffold_id}_qc_pass.rf', scaffold_id=scaffold_ids),
        expand('sites/all_vulpes_scaffold_{scaffold_id}_qc_fail.bz', scaffold_id=scaffold_ids),
        expand('sites/all_vulpes_scaffold_{scaffold_id}_qc_pass.vcf', scaffold_id=scaffold_ids)

rule run_snp_cleaner:
    input:
        'all_vulpes_bams.txt'
    output:
        'sites/all_vulpes_scaffold_{scaffold_id}_qc_pass.rf',
        'sites/all_vulpes_scaffold_{scaffold_id}_qc_fail.bz',
        'sites/all_vulpes_scaffold_{scaffold_id}_qc_pass.vcf'
    shell: """
    bcftools mpileup --threads 20 -f {ref_fasta} -b {input} -r scaffold_{wildcards.scaffold_id} -a SP,DP | bcftools call --skip-variants indels -f GQ -c -| /home/joana/software/: -Q 30 -k 43 -u 1 -h 0 -H 1e-4 -b 1e-20 -S 1e-4 -f 1e-4 -e 1e-4 -v -B {output[0]} -p {output[1]} | bzip2 -c > {output[2]}
    """

#old option: 
#bcftools mpileup --threads 20 -f {ref_fasta} -b {input} -r scaffold_{wildcards.scaffold_id} -a SP,DP | bcftools call --skip-variants indels -f GQ -c -| /home/joana/software/ngsQC/snpCleaner/snpCleaner.pl -D 4250 -k 68 - u 1 -h 0 -H 1e-4 -b 1e-20 -S 1e-4 -f 1e-4 -e 1e-4 -v -B {output[0]} -p {output[1]} > {output[2]}
 
 #snpCleaner - After taliking with Tyler
#-D  50x times 85 individuals -> max raw site read depth 
#-k 1 site is considered covered if at least 80% individuals are covered 
#-u 1 # individual considered covered if has >=1 reads at a site 
#-H 10^-4 (10^-6 is what Deb uses); setting it to lower will make it more relaxed
#-h 0 min p-value for exact test of HWE [0.0001]
#-b 1e-20 min p-value for base quality bias] # -b 1e-10 Deb, while Tyler used default; Tyler advised me to use -b 1e-20
#-S 1e-4  min p-value for strand bias [0.0001]
#-f 1e-4  min p-value for map quality bias [0]  
#-e 1e-4  min p-value for end distance bias [0.0001]
#-v option to keep all sites

rule group_scaffolds:
    input: 
        'sites/{wildcards.scaffold_id}_qc_pass.rf',
        'sites/{wildcards.scaffold_id}_qc_fail.bz2',
        'sites/{wildcards.scaffold_id}_qc_pass.vcf.bz2'
    output
        temp('sites/unsorted_sites.rf'),
        'sites/all_vulpes_sites_qc_fail.bz2',
        'sites/all_vulpes_sites_qc_pass.vcf.bz2',
        'sites/all_vulpes_sites_qc_pass.rf'
    shell:  """
    cat {input[0]} > {output[0]} &&
    cat {input[1]} > {output[1]} &&
    cat {input[2]} > {output[2]} &&
    sort -k1,1 -V -s {output[0]} > {output[3]}
    """


#Before running angsd 2 things are needed:
#position file it will have the list of scaffolds and the length of the scaffolds in column 1 and 2, respectively. This will be given to angsd -sites option
#the following rules can be done in sites and snps directory, depending on either one choose to retain or not non variants

rule create_pos_file: 
    input: 
        'sites/all_vulpes_sites_qc_pass.rf'
    output:
        'sites/all_vulpes_sites_qc_pass.pos'
    shell: """
    /home/joana/software/ngsQC/vcfCleaner/expandrf pos {input} > {output} &&   
    angsd sites index {output}
    """
#bzcat allchimps_stacy_sites_qc_pass.vcf.bz2 | bcftools query -f '%CHROM %POS\n'  > allchimps_stacy_sites_qc_pass.pos




#list of contigs/scaffolds file -> list of scaffolds/chromosomes 
rule create_list_scaffolds:
    input: 
        'sites/all_vulpes_sites_qc_pass.pos'
    output:
        'sites/all_vulpes_sites_qc_pass_scaffolds.rf'
    shell: 'cut -f 1 {input} | uniq | sed 's/$/:/g' > {output}'
