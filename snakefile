#
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
javaSingularity       = config["javaSingularity"]
samtoolsSingularity   = config["samtoolsSingularity"]
#variantUtilsSingularity   = config["variantutilsSingularity"]
varscanLib            = config["varscanLib"]
picardLib             = config["picardLib"]

#
sampleList   = ["normal","tumour"]
normal       = "normal"
tumourSamples = ["tumour"]

#
#
#
# Outputs:
markDuplicateBams   = intermediateDirectory + "/{sample}.markDuplicates.bam"
somaticSnps         = resultsDirectory + "/{tumourSample}.varscan.somatic.snp.vcf"
markedSortedIndexes = resultsDirectory + "/{sample}.markDuplicates.sorted.bam.bai"


somaticMergedVCF = intermediateDirectory + "/{sample}.somatic.merged.vcf"
decomposedVCF = intermediateDirectory + "/{sample}.decomposed.vcf"
normalisedVCF = intermediateDirectory + "/{sample}.normalised.vcf"
normalisedVCFgz = intermediateDirectory + "/{sample}.normalised.vcf.gz"
annotatedVCF = intermediateDirectory + "/{sample}.annotated.vcf"
filteredVCF = intermediateDirectory + "/{sample}.filtered.vcf"


vcfannoConfig = "resources/conf.toml"


#envvars:
#    "SINGULARITY_BIND"

report: "report/workflow.rst"

# Currently desired outputs
all_results = expand([filteredVCF,annotatedVCF,normalisedVCFgz,normalisedVCF,decomposedVCF,somaticMergedVCF,markDuplicateBams, markedSortedIndexes, somaticSnps], sample = sampleList, tumourSample = tumourSamples )


rule all:
    input:
        all_results

######################### Sort Bams ##########################################
rule sort:
	input:
		unsortedBam = rawDataDirectory + "/{sample}.bam"
	output:
		sortedBam = intermediateDirectory + "/{sample}.sorted.bam"
#		sortedBam = report(intermediateDirectory + "/{sample}.sorted.bam",category="BAM files")
	log:
		logsDirectory + "/{sample}.sort.log"
	threads:
		2
#	resources:
#		mem = "16GB",
#		time = "04:00:00"
	container:
		singularityDirectory + "/" + samtoolsSingularity
	shell:
		"samtools sort -@ {threads} -m 400M -O BAM -o {output.sortedBam} {input.unsortedBam} >{log} 2>&1"


######################### index ##########################################
rule index:
    input:
        sortedBam = intermediateDirectory  + "/{sample}.sorted.bam"
    output:
        bamIndex = intermediateDirectory  + "/{sample}.sorted.bam.bai"
    wildcard_constraints:
        sample = '|'.join([x for x in sampleList ]) # prevents sample being read as a regular expression
    log:
        logsDirectory + "/{sample}.index.log"
    threads:
       2
#    resources:
#        mem = "8GB",
#        time = "00:30:00"
    container:
        singularityDirectory + "/" + samtoolsSingularity
    shell:
        "samtools index -@ {threads} -m 400M {input.sortedBam} {output.bamIndex} >{log} 2>&1"


######################### Mark Duplicates ####################################
rule markDuplicates:
	input:
		sortedBam = intermediateDirectory + "/{sample}.sorted.bam",
		index = intermediateDirectory + "/{sample}.sorted.bam.bai"
	output:
		bam = intermediateDirectory + "/{sample}.markDuplicates.bam",
#		metrics = metricsDirectory + "/{sample}.metrics.txt"
		metrics = report(metricsDirectory + "/{sample}.metrics.txt",category="Metrics")
	log:
		logsDirectory + "/{sample}.markDuplicates.log"
	threads:
		2
#	resources:
#		mem = "32GB",
#		time = "01:00:00",
#		partition = "bigmem",
#		queue = "bigmem"
	container:
		singularityDirectory + "/" + javaSingularity
	shell:
		"java $JVM_ARGS -jar {libDirectory}/picard.jar MarkDuplicates I={input.sortedBam} O={output.bam} METRICS_FILE={output.metrics} QUIET=true SORTING_COLLECTION_SIZE_RATIO=0.125 MAX_RECORDS_IN_RAM=250000 >{log} 2>&1"

######################### Mark Duplicates Sort ####################################
rule markedDuplicatesSort:
	input:
		unsortedBam = rules.markDuplicates.output.bam
	output:
		sortedBam = resultsDirectory + "/{sample}.markDuplicates.sorted.bam",
	log:
		logsDirectory + "/{sample}.markDuplicatesSort.log"
	threads:
		2
#	resources:
#		mem = "8GB",
#		time = "01:00:00"
	container:
		singularityDirectory + "/" + samtoolsSingularity
	shell:
		"samtools sort -@ {threads} -m 400M -O BAM -o {output.sortedBam} {input.unsortedBam} >{log} 2>&1"


######################### Mark Duplicates index ####################################
rule markedDuplicatesIndex:
    input:
        sortedBam = rules.markedDuplicatesSort.output.sortedBam
    output:
        bamIndex = resultsDirectory  + "/{sample}.markDuplicates.sorted.bam.bai"
    log:
        logsDirectory + "/{sample}.markDuplicatesIndex.log"
    container:
        singularityDirectory + "/" + samtoolsSingularity
    threads:
       2
#    resources:
#        mem = "16GB",
#        time = "04:00:00"
    shell:
        "samtools index -@ {threads} -m 400M {input.sortedBam} {output.bamIndex} >{log} 2>&1"

########################################### Pile up #############################################
rule somaticPileUp:
    input:
        normalBam = resultsDirectory + "/" + normal + ".markDuplicates.sorted.bam",
        tumourBam = resultsDirectory + "/{tumourSample}.markDuplicates.sorted.bam",
        reference = referenceDirectory + "/" + referenceFile
    output:
        somaticPileUp = intermediateDirectory + "/{tumourSample}.somatic.pileup"
#        somaticPileUp = report(intermediateDirectory + "{tumourSample}.somatic.pileup",category="Somatic Pileup")
    wildcard_constraints:
        tumourSample = '|'.join([x for x in sampleList ]) # prevents sample being read as a regular expression
    log:
        logsDirectory + "/{tumourSample}.somaticPileup.log"
    threads:
        2
#    resources:
#        mem = "16GB",
#        time = "10:00:00"
    container:
        singularityDirectory + "/" + samtoolsSingularity
    shell:
        "samtools mpileup -P {threads} -B -d 9001 -q 1 -f {input.reference} -o {output.somaticPileUp} {input.normalBam} {input.tumourBam} >{log} 2>&1"

########################################## Call Variants #############################################

rule varscanSomatic:
    input:
        pileup = rules.somaticPileUp.output.somaticPileUp
    output:
        varscanSnpVcf = resultsDirectory + "/{tumourSample}.varscan.somatic.snp.vcf",
        varscanIndelVcf = resultsDirectory + "/{tumourSample}.varscan.somatic.indel.vcf"
#        varscanSnpVcf = report(intermediateDirectory + "/{tumourSample}.varscan.somatic.snp.vcf",caption="report/caption_snpvcf.rst",category="VCF files"),
#        varscanIndelVcf = report(intermediateDirectory + "/{tumourSample}.varscan.somatic.indel.vcf",caption="report/caption_indelvcf.rst",category="VCF files")
    log:
        logsDirectory + "/{tumourSample}.varscanSomatic.log"
    threads:
        2
#    resources:
#        mem = "16GB",
#        time = "10:00:00"
    container:
        singularityDirectory + "/" + javaSingularity
    shell:
         "java $JVM_ARGS -jar {libDirectory}/{varscanLib} somatic {input.pileup} --min-var-freq 0.1 --p-value 1.00 --somatic-p-value 1.0 --strand-filter 0 --tumour-purity 0.5 --output-vcf 1 --min-coverage-normal 10 --min-coverage-tumor 10 --mpileup --output-snp {output.varscanSnpVcf} --output-indel {output.varscanIndelVcf} >{log} 2>&1"


#
# New stuff
###################################### Decompose ######################################

#include: "../../tools/mergeVcf.snk
rule mergeVcf:
    input:
        snp = rules.varscanSomatic.output.varscanSnpVcf,
        indel = rules.varscanSomatic.output.varscanIndelVcf
    output:
        vcf = intermediateDirectory + "/{tumourSample}.somatic.merged.vcf"
    log:
        mergeVcfLog = "logs/{tumourSample}.mergeVcf.log"
    threads:
        2
    resources:
        mem = "8GB",
        time = "01:00:00"
    container:
        #"docker://mcfonsecalab/variantutils:0.6"
        singularityDirectory + "/variantutils_0.6.sif"
    shell:
        """
        bgzip -@ {threads} -c {input.snp} > {input.snp}.gz && \
        tabix -p vcf {input.snp}.gz && \
        bgzip -@ {threads} -c {input.indel} > {input.indel}.gz && \
        tabix -p vcf {input.indel}.gz && \
        bcftools merge --force-samples {input.snp}.gz {input.indel}.gz -o {output.vcf} > {log.mergeVcfLog}
        """

###################################### Decompose ######################################

#include: "../../tools/decompose.snk"
rule decompose:
    input:
        reference = referenceDirectory + "/" + referenceFile,
        vcf = rules.mergeVcf.output.vcf
    output:
        vcf = intermediateDirectory + "/{tumourSample}.decomposed.vcf"
    log:
        decompositionLog = "logs/{tumourSample}.decomposition.log"
    threads:
        2
    resources:
        mem = "8GB",
        time = "01:00:00"
    container:
        #"docker://mcfonsecalab/variantutils:0.6"
        singularityDirectory + "/variantutils_0.6.sif"
    shell:
        """
        vt decompose -s {input.vcf} -o {output.vcf} > {log.decompositionLog}
        """

###################################### Normalise ######################################

#include: "../../tools/normalise.snk"
rule normalise:
    input:
        reference = referenceDirectory + "/" + referenceFile,
        vcf = rules.decompose.output.vcf
    output:
        vcf = intermediateDirectory + "/{tumourSample}.normalised.vcf"
    log:
        normalisationLog = "logs/{tumourSample}.normalisation.log"
    threads:
        2
    resources:
        mem = "8GB",
        time = "01:00:00"
    container:
        #"docker://mcfonsecalab/variantutils:0.6"
        singularityDirectory + "/variantutils_0.6.sif"
    shell:
        """
        vt normalize {input.vcf} -r {input.reference} -o  {output.vcf} > {log.normalisationLog}
        """

###################################### Normalise_index vcf calls ######################################################
rule normalise_index:
    input:
        vcf = rules.normalise.output.vcf
    output:
        vcfgz = intermediateDirectory + "/{tumourSample}.normalised.vcf.gz",
        vcftbi = intermediateDirectory + "/{tumourSample}.normalised.vcf.gz.tbi",
    log:
        normalisationLog = "logs/{tumourSample}.normalisation.log"
    threads:
        2
#    resources:
#        mem = "80GB",
#        time = "01:00:00"
    container:
        singularityDirectory + "/vt.simg"
    shell:
        """
        bgzip -@ {threads} -c {input.vcf} > {output.vcfgz} && tabix -p vcf {output.vcfgz}
        """

###################################### Annotate ######################################

#include: "../../tools/annotate.snk"
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
        annotatedVcf = intermediateDirectory + "/{tumourSample}.annotated.vcf",
    log:
        annotationLogs = "logs/{tumourSample}.annotation.log"
    threads:
        2
    resources:
        mem = "8GB",
        time = "01:00:00"
    container:
        #"docker://mcfonsecalab/variantutils:0.6"
        singularityDirectory + "/variantutils_0.6.sif"
    shell:
        """
            vcfanno -p {threads} {vcfannoConfig} {input.vcf} > {output.annotatedVcf} 2>{log.annotationLogs}
        """


###################################### Filter ######################################

#include: "../../tools/filter.snk"
rule filter:
    input:
        vcf = rules.annotate.output.annotatedVcf
    output:
        filteredVcf = intermediateDirectory + "/{tumourSample}.filtered.vcf",
    log:
        filteringLogs = "logs/{tumourSample}.filtering.log"
    threads:
        2
    resources:
        mem = "8GB",
        time = "01:00:00"
    container:
        #"docker://mcfonsecalab/variantutils:0.6"
        singularityDirectory + "/variantutils_0.6.sif"
    shell:
            #bcftools view -e'clinvarDiseaseName=""' -f'* [ %clinvarDiseaseName]\n | gnomADAlleleCount=""' -f'* [ %gnomadADAlleCount]\n ' {input.vcf} -o {output.filteredVcf} > {log.filteringLogs}
        """
            bcftools view {input.vcf} -o {output.filteredVcf} > {log.filteringLogs}
        """

