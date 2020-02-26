# MutSig2CV_NC (PCAWG)

MutSig2CV ([Lawrence *et al.*, 2014](https://www.nature.com/articles/nature12912)) adapted for noncoding 
significance analysis, as run for the PCAWG drivers paper ([Rheinbay *et al.*, 2020](https://doi.org/10.1038/s41586-020-1965-x)).

## Installing

MutSig is implemented in MATLAB. If you have a MATLAB installation and wish to run MutSig interactively on the MATLAB console, skip to the [Running](#running) section below. If you do not have MATLAB installed, or do not wish to run interactively, MutSig can be run as a standalone executable. The standalone executable is available for 64 bit Linux systems only, and requires that the MATLAB R2016a runtime be installed.
You can download and install the runtime environment from [here](https://ssd.mathworks.com/supportfiles/downloads/R2016a/deployment_files/R2016a/installers/glnxa64/MCR_R2016a_glnxa64_installer.zip). Runtime installation instructions can be found [here](http://www.mathworks.com/help/compiler/install-the-matlab-runtime.html).

Once the runtime is successfully installed, you must add it to your `LD_LIBRARY_PATH`.

```bash
MCRROOT=<path to runtime you specified when installing>
export LD_LIBRARY_PATH=.:${MCRROOT}/runtime/glnxa64
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/opengl/lib/glnxa64
```

To run on the cohort in the manuscript, MutSig requires ~100 GB of reference files. This is due to the fact that each cohort has a unique set of covariates and coverage models (see "[Cohort Parameters](#cohort_params)" below). Since these files are too large to include in a GitHub repository, they are hosted in cloud storage; please see the [README in the `ref/` folder](blob/master/ref/README.md) for instructions on how to obtain them.

## Running <a name="running"></a>

To run on the MATLAB console, start MATLAB in this directory, and run:

```matlab
MutSig2CV_NC(<path to mutations>, <path to output directory>, <path to cohort parameters>)
```

To run the standalone application, `cd` to this directory, and run:

```bash
bin/MutSig2CV_NC <path to mutations> <path to output directory> <path to cohort parameters>
```

MutSig looks for its reference files relative to this directory, so it is essential it is run here.

Each input is decribed below.

### Description of inputs

* **Mutations file**: Absolute path to the set of mutations to analyze.  Format is specified in the
following section, [Mutation Input Format](#mutation_inputs)

* **Output directory**: Absolute path to the directory where MutSig will save its output.  Will be
created if necessary.  NB: Any previous MutSig results in this directory will
be overwritten!  <!-- A description of the output files and formats are in the
following section, [Output Format] -->

* **Cohort parameters**: <a name="cohort_params"></a> Each cohort and noncoding feature (e.g., promoters, UTRs, enhancers, etc.) analyzed in the manuscript has its own set of reference files. This is due to the fact that genomic coverage in noncoding regions can vary substantially from cohort to cohort, and MutSig requires an accurate estimate of sequencing coverage to properly compute mutation densities. MutSig must be explicitly told which cohort/feature-specific reference files to use, which are located in `run/params`. For example, to analyze transcription factor binding sites in the glioblastoma cohort, use `run/params/CNS-GBM_TFBS.params.txt`. Cohort parameter files' names follow the convention `<cohort>_<feature>.txt`, where `<cohort>` corresponds to the "Tier 3 Abbreviations" in the [PCAWG tumor subtype table](https://dcc.icgc.org/api/v1/download?fn=/PCAWG/clinical_and_histology/tumour_subtype_consolidation_map.xlsx), and `<feature>` corresponds to the following noncoding features: <!--  * `alt_prom` -->
  * `enhancers`: Enhancers
  * `gc_pc.3utr`: 3' UTRs
  * `gc_pc.5utr`: 5' UTRs
  * `lncrna`: lncRNA genes
  * `lncrna.prom`: lncRNA promoters
  * `mirna.mat`: Mature miRNAs
  * `mirna.pre`: Pre-miRNAs
  * `mirna.prom`: miRNA promoters <!-- * promcore:  -->
  * `promoters`: Coding gene promoters <!-- * smallrna.ncrna -->
  * `TFBS`: Transcription Factor Binding Sites

### Mutation Input Format <a name="mutation_inputs"></a>

As input, MutSig takes a tab-delimited file with each line annotating a single
mutation in a single patient.  Columns can be in any order, with names and
formats as follows.  To provide maximal input flexibility, MutSig accepts
synonyms for each column name.  Column names are case sensitive.

* `chr`: Chromosome of the mutation.  MutSig only analyzes mutations on autosomal or sex chromosomes, and does not consider the mitochondrial chromosome or unplaced/alternate contigs.
  * Range: `(chr)?[1..24XY]`
  * Synonyms: `Chromosome`

* `pos`: hg19 position of the mutation, 1-indexed.
  * Regex: `[0-9]+`
  * Synonyms: `Position`, `start`, `Start_position`

* `patient`: Unique identifier for the patient.
  * Regex: `[A-Za-z0-9]+`
  * Synonyms: `Tumor_Sample_Barcode`, `Patient_name`

* `ref_allele`: hg19 reference base(s) for the position.  In the case of insertions, must be "-".
  * Regex: `(-|[ACGT]+)`
  * Synonyms: `Reference_Allele`

* `newbase`: Observed variant allele at the position.  In the case of deletions, must be "-".
  * Regex: `(-|[ACGT]+)`
  * Synonyms: `Tumor_Allele`, `Tum_allele`, `Alt_allele`, `Alternate_allele`, `Tumor_Seq_Allele2`
  
Note that MutSig does not require any other mutation annotations; it infers everything else on its own.
