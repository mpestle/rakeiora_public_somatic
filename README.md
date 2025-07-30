# Somatic Workflow - kinghorn1 dataset version

[Rakeiora](http://rakeiora.ac.nz)

---
Somatic workflow - tumour vs normal.

---

This workflow is a combination of Ben Curran's original somatic workflow
(https://github.com/nesi/rakeiora-public-somatic)
and Peter Tsai's scripts to do a somatic analysis.

This workflow utilises a pipe between the
pileup and varscan functions and therefore eliminates a massive amount of
interim storage on disk that Ben's workflow required.
Parallelism is achieved by partitioning the mpileup and varscan steps
by chromosome (1-22, X, Y, M).

This is running in the Rakeiora environment on a 16cpu32ram VM
in about 3 hours, when snakemake is given 16 cores.

Before running, the calling environment needs to set (and export) the
SINGULARITY_BIND environment variable, which needs to include the
/shared resources and the attached volume with the dataset, so
after you get the volume attached:

```export SINGULARITY_BIND=/shared,/rv/kinghorn1```

Then

```snakemake -c16 --use-singularity```

Utilises a singularity specifically built for this workflow which contains:
- SAMTOOLS_VERSION=1.19
- BCFTOOLS_VERSION=1.19
- VCFANNO_VERSION=v0.3.6
- Java 20

See singularity.def

To build this singularity:

```singularity build --fakeroot rakeiora-apps-20250730.sif singularity.def```

Two libraries are used from the /shared/lib area,
VarScan.v2.4.6.jar and picard.jar (and the picard isn't actually used),
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
