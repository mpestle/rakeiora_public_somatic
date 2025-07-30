# Somatic Workflow

[Rakeiora](http://rakeiora.ac.nz)

---
Somatic workflow - tumour vs normal.

---

Unsure of the function of this workflow, but it does a bunch
of analysis of a tumour sample vs a normal sample. Perhaps
someone with more knowledge of these tools can update this
readme file somewhere along the line.

Utilises a singularity built for this workflow which uses
- SAMTOOLS_VERSION=1.19
- BCFTOOLS_VERSION=1.19
- VCFANNO_VERSION=v0.3.6
- Java 20

See singularity.def
To build this singularity:
```singularity build --fakeroot rakeiora-apps-20250730.sif singularity.def```

Two libraries are used from the /shared/lib area,
VarScan.v2.4.6.jar and picard.jar (and actually the picard isn't actually used),
along with a shared
reference file in /shared/reference/,
GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna

The two input files expected are named tumour.bam and normal.bam,
in the resources directory, so the definition of the inputs
in the Rakeiora sandbox should define these accordingly,
presumably with the tumour sample related to/dependent on
the previously selected normal one.

Note that a bunch of vcfs come back into the results directory,
and the contents of that directory will be copied to the jupyter hub
area when the researcher runs a workflow in Rakeiora. So before
this is run on real data, the kaitiaki of the data
and the workflow reviewers need to be OK with that.

NeSI does not claim ownership or authorship of this workflow.
Authorship credits will follow later.
