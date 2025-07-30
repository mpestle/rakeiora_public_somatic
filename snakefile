# Somatic workflow for Cris Print/ Peter Tsai kinghorn1 dataset.
# Authors: Ben Curran/Peter Tsai/Matt Pestle
# Date: Jul 2025
#
# Try running this on a VM with N cpus and giving snakemake 2N cores
# So if cpuinfo shows a count of 4, try
#      snakemake -c8 --use-singularity
# Timings:
# 4cpu Rakeiora system with 8 snakemake cores  : Under 5hrs for a tumour/normal pair
# 16cpu Rakeiora system with 16 snakemake cores: Under 3hrs for a tumour/normal pair

configfile: "inputs.yaml"
#
rawDataDirectory      = config["rawDataDirectory"]
intermediateDirectory = config["intermediateDirectory"]
resultsDirectory      = config["resultsDirectory"]
referenceDirectory    = config["referenceDirectory"]
referenceFile         = config["referenceFile"]
singularityDirectory  = config["singularityDirectory"]
libDirectory          = config["libDirectory"]
metricsDirectory      = config["metricsDirectory"]
logsDirectory         = config["logsDirectory"]
workflowSingularity   = config["workflowSingularity"]
varscanLib            = config["varscanLib"]
vcfannoConfig         = "resources/conf.toml"
#
normal = "normal"
tumour = "tumour"
#
#
# Outputs:
varscanSnpVcf         = intermediateDirectory + "/{tumour}.{chr}.varscan.somatic.snp.vcf",
varscanSnpVcfgz       = intermediateDirectory + "/{tumour}.{chr}.varscan.somatic.snp.vcf.gz",
varscanSnpVcfgzCsi    = intermediateDirectory + "/{tumour}.{chr}.varscan.somatic.snp.vcf.gz.csi",
varscanIndelVcf       = intermediateDirectory + "/{tumour}.{chr}.varscan.somatic.indel.vcf",
varscanIndelVcfgz     = intermediateDirectory + "/{tumour}.{chr}.varscan.somatic.indel.vcf.gz",
varscanIndelVcfgzCsi  = intermediateDirectory + "/{tumour}.{chr}.varscan.somatic.indel.vcf.gz.csi",
mergedVCF             = resultsDirectory      + "/tumour_merged.vcf.gz"
decomposedVCF         = intermediateDirectory + "/tumour.decomposed.vcf"
normalisedVCF         = intermediateDirectory + "/tumour.normalised.vcf"
normalisedVCFgz       = intermediateDirectory + "/tumour.normalised.vcf.gz"
normalisedVCFgztbi    = intermediateDirectory + "/tumour.normalised.vcf.gz.tbi"
annotatedVCF          = intermediateDirectory + "/tumour.annotated.vcf"
filteredVCF           = intermediateDirectory + "/tumour.filtered.vcf"


#################################################################################
################ SINGULARITY_BIND environment variable ##########################
#################################################################################
# The Singularity needs access to some file systems, which need to be specified
# here in the SINGULARITY_BIND environment variable.
# It needs to see the /shared file system and the data: /rv/kinghorn1 (or wherever mounted)
# so that the symlinks to files on that data system can be followed.
# So in the calling environment:
# export SINGULARITY_BIND=/shared,/rv/kinghorn1
envvars:
    "SINGULARITY_BIND"
#################################################################################

# We will achieve parallelisation by partitioning into the 25 different chromosomes
# and letting snakemake handle processes that can be run concurrently.
CHROMOSOMES = [f"chr{i}" for i in range(1, 23)]
CHROMOSOMES = CHROMOSOMES + ["chrM","chrX","chrY"]

# Currently desired outputs - probably more than is required
all_results = [
	 expand(varscanSnpVcf, tumour=tumour, chr=CHROMOSOMES)
	,expand(varscanSnpVcfgz, tumour=tumour, chr=CHROMOSOMES)
	,expand(varscanSnpVcfgzCsi, tumour=tumour, chr=CHROMOSOMES)
	,expand(varscanIndelVcf, tumour=tumour, chr=CHROMOSOMES)
	,expand(varscanIndelVcfgz, tumour=tumour, chr=CHROMOSOMES)
	,expand(varscanIndelVcfgzCsi, tumour=tumour, chr=CHROMOSOMES)
	,mergedVCF
	,decomposedVCF 
	,normalisedVCF
	,normalisedVCFgz 
	,normalisedVCFgztbi
	,annotatedVCF
	,filteredVCF
	]

rule all:
    input:
        all_results

####################### pileup and varscan pipe ###################################
rule varscanSomatic:
    input:
        normalBam = rawDataDirectory   + "/normal.bam",
        tumourBam = rawDataDirectory   + "/tumour.bam",
        reference = referenceDirectory + "/" + referenceFile,
        varscanlib = libDirectory + "/" + varscanLib
    output:
        varscanSnpVcf   = intermediateDirectory + "/tumour.{chr}.varscan.somatic.snp.vcf",
        varscanIndelVcf = intermediateDirectory + "/tumour.{chr}.varscan.somatic.indel.vcf"
    log:
        logsDirectory + "/tumour.{chr}.varscanSomatic.log"
    threads:
        2
#    resources:
#        mem = "16GB",
#        time = "10:00:00"
    container: singularityDirectory + "/" + workflowSingularity
    shell:
        """
	samtools mpileup \
	-r {wildcards.chr} \
	-B \
	-d 9001 \
	-q 1 \
	-f {input.reference} \
	{input.normalBam} {input.tumourBam} \
	| \
	java \
	-jar {input.varscanlib} somatic \
	--min-var-freq 0.1 \
	--p-value 1.00 \
	--somatic-p-value 1.0 \
	--strand-filter 0 \
	--tumour-purity 0.5 \
	--output-vcf 1 \
	--min-coverage-normal 10 \
	--min-coverage-tumor 10 \
	--mpileup 1 \
	--output-snp {output.varscanSnpVcf} \
	--output-indel {output.varscanIndelVcf} \
	>{log} 2>&1
        """
###################################### Gzip snp varscans ######################################
rule GzipSnp:
    input:
        f = rules.varscanSomatic.output.varscanSnpVcf
    output:
        vcfgz =  rules.varscanSomatic.output.varscanSnpVcf + ".gz"
    container: singularityDirectory + "/" + workflowSingularity
    shell:
        """
            bgzip -c {input.f} > {input.f}.gz
        """
###################################### Gzip indel varscans ######################################
rule GzipIndel:
    input:
        f = rules.varscanSomatic.output.varscanIndelVcf
    output:
        vcfgz =  rules.varscanSomatic.output.varscanIndelVcf + ".gz"
    container: singularityDirectory + "/" + workflowSingularity
    shell:
        """
            bgzip -c {input.f} > {input.f}.gz
        """
###################################### Gzip snp index ######################################
rule GzipSnpIndex:
    input:
        f = rules.GzipSnp.output.vcfgz
    output:
        vcfcsi =  rules.GzipSnp.output.vcfgz + ".csi"
    container: singularityDirectory + "/" + workflowSingularity
    shell:
        """
            bcftools index {input.f}
        """
###################################### Gzip indel index ######################################
rule GzipIndelIndex:
    input:
        f = rules.GzipIndel.output.vcfgz
    output:
        vcfcsi =  rules.GzipIndel.output.vcfgz + ".csi"
    container: singularityDirectory + "/" + workflowSingularity
    shell:
        """
            bcftools index {input.f}
        """
###################################### Merge Indel ######################################
rule Merge:
    input:
        fs  = expand("int/tumour.{chr}.varscan.somatic.snp.vcf.gz",chr=CHROMOSOMES),
        fi  = expand("int/tumour.{chr}.varscan.somatic.indel.vcf.gz",chr=CHROMOSOMES),
        fndx = expand("int/tumour.{chr}.varscan.somatic.indel.vcf.gz.csi",chr=CHROMOSOMES),
        fsdx = expand("int/tumour.{chr}.varscan.somatic.snp.vcf.gz.csi",chr=CHROMOSOMES)
    output:
        vcf =  resultsDirectory + "/tumour_merged.vcf.gz"
    container: singularityDirectory + "/" + workflowSingularity
    shell:
        """
            bcftools concat -a {input.fs} {input.fi} | \
		bcftools sort --write-index -Oz -o {output.vcf}
        """
##################################### Decompose ######################################
rule decompose:
    input:
        vcf = rules.Merge.output.vcf
    output:
        vcf = decomposedVCF
    log:
        decompositionLog = "logs/tumour.decomposition.log"
    threads:
        2
    resources:
        mem = "8GB",
        time = "01:00:00"
    container: singularityDirectory + "/" + workflowSingularity
    shell:
        """
        vt decompose -s {input.vcf} -o {output.vcf} > {log.decompositionLog}
        """
##################################### Normalise ######################################
rule normalise:
    input:
        reference = referenceDirectory + "/" + referenceFile,
        vcf = rules.decompose.output.vcf
    output:
        vcf = normalisedVCF
    log:
        normalisationLog = "logs/tumour.normalised.log"
    threads:
        2
    resources:
        mem = "8GB",
        time = "01:00:00"
    container: singularityDirectory + "/" + workflowSingularity
    shell:
        """
        vt normalize {input.vcf} -r {input.reference} -o  {output.vcf} > {log.normalisationLog}
        """
##################################### Normalise_index vcf calls ######################################################
rule normalise_index:
    input:
        vcf = rules.normalise.output.vcf
    output:
        vcfgz = normalisedVCFgz,
        vcftbi = normalisedVCFgztbi
    log:
        normalisationLog = "logs/tumour.normalisedndx.log"
    threads:
        2
#    resources:
#        mem = "80GB",
#        time = "01:00:00"
    container: singularityDirectory + "/" + workflowSingularity
    shell:
        """
        bgzip -@ {threads} -c {input.vcf} > {output.vcfgz} && tabix -p vcf {output.vcfgz}
        """
##################################### Annotate ######################################
# function needs to go here to enable this stage to pick up inpupt from either germline or somatic threads.
# input within the rule will accept a function call/values returned by a function call as inputs.
# This should possibly go in the decompose stage?

def get_input(wildcards):
    input_list = []
    if config["DEG"]["exec"]:
          input_list.append("DEG/DEG.txt")
    if config["DTU"]["exec"]:
          input_list.append("DTU/DTU.txt")
    return input_list

rule annotate:
    input:
        vcf = rules.normalise_index.output.vcfgz
    output:
        annotatedVcf = annotatedVCF
    log:
        annotationLogs = "logs/tumour.annotation.log"

    # 8 threads is about twice as fast as 2: 27 min vs 15
    threads:
        8
    resources:
        mem = "8GB",
        time = "01:00:00"
    container: singularityDirectory + "/" + workflowSingularity
    shell:
        """
            vcfanno -p {threads} {vcfannoConfig} {input.vcf} > {output.annotatedVcf} 2>{log.annotationLogs}
        """
##################################### Filter ######################################
rule filter:
    input:
        vcf = rules.annotate.output.annotatedVcf
    output:
        filteredVcf = filteredVCF
    log:
        filteringLogs = "logs/tumour.filtering.log"
    threads:
        2
    resources:
        mem = "8GB",
        time = "01:00:00"
    container: singularityDirectory + "/" + workflowSingularity
    shell:
        """
            bcftools view -e'clinvarDiseaseName=""' -f'* [ %clinvarDiseaseName]\n | gnomADAlleleCount=""' -f'* [ %gnomadADAlleCount]\n ' {input.vcf} -o {output.filteredVcf} > {log.filteringLogs}
            #bcftools view {input.vcf} -o {output.filteredVcf} > {log.filteringLogs}
        """
