# GWAS / PRS Nextflow Pipeline

This repository contains a modular Nextflow pipeline for:

- genotype QC
- ancestry PCA against a local 1000 Genomes reference
- GWAS association testing
- PRS calculation with PRSice-2
- result plotting and summary

This README includes:

1. what the pipeline does
2. repository layout
3. environment setup
4. data requirements
5. how to build/download reference data
6. how to run the validated simulated chr22 test
7. recommended resource allocation
8. how to run a real full study
9. what must change when moving from test to full run
10. optional simulation validation
11. common errors and fixes
12. version-control guidance

---

# 1. Current status

This pipeline has been validated end-to-end on a **simulated chr22 test run** using:

- a real chr22 1000 Genomes genotype backbone
- simulated phenotype/covariates
- FinnGen R12 lung adenocarcinoma summary statistics

Validated components:

- `QC`
- `ANCESTRY_PCA`
- `ASSOCIATION`
- `PLOTS`
- `PRS`

A successful run completed with:

- QC completed
- ancestry PCA completed
- association completed
- association plots completed
- PRS completed

This is a **workflow validation**, not a biological validation of lung cancer results.

---

# 2. Repository layout

Main repository location:

```bash
/data/rezaa2/gwas-pipeline
```

Typical structure:

```text
gwas-pipeline/
├── conf/
├── data/
├── modules/
├── scripts/
├── main.nf
├── nextflow.config
├── run_download_1kg_ref.sbatch
├── run_download_finngen_r12_lung.sbatch
├── run_full_verify.sbatch
├── simulate_and_run.sh
├── simulated_ld_data/           # generated test data
├── results_sim_ld_chr22/        # generated test outputs
└── README.md
```

Important tracked files:

- `main.nf`
- `nextflow.config`
- `conf/`
- `modules/`
- `scripts/`
- `run_*.sbatch`
- small example files under `data/`

Generated outputs that should usually **not** be tracked:

- `.nextflow/`
- `work/`
- `results*/`
- `prs_results/`
- `simulated_ld_data/`
- logs
- zip files
- generated PLINK outputs

---

# 3. Environment setup

This pipeline was developed on NIH Biowulf.

Typical setup:

```bash
cd /data/rezaa2/gwas-pipeline
source ~/bin/myconda
conda activate gwas_plink19
module load nextflow
module load R
```

Sanity checks:

```bash
which plink
which nextflow
which Rscript
python3 --version
```

Expected:
- `plink` from conda env
- `nextflow` from module
- `Rscript` from module
- `python3` available

Notes:
- `PRSice` is loaded inside the PRS process with `module load prsice`
- `R` is required for helper scripts and plots

---

# 4. Data requirements

## 4.1 Minimum required for a real full run

You need:

1. a study cohort in PLINK format:
   - `cohort.bed`
   - `cohort.bim`
   - `cohort.fam`

2. a phenotype CSV with at least:
   - `FID`
   - `IID`
   - `PHE`

3. a base GWAS summary statistics file for PRS

4. an ancestry / LD reference panel, typically local 1000 Genomes

## 4.2 Phenotype file format

Expected schema:

```text
FID,IID,PHE,age,sex,batch
sample1,sample1,2,64,2,batch1
sample2,sample2,1,58,1,batch2
```

Rules:

- `PHE`: `1=control`, `2=case`
- `sex`: `1=male`, `2=female`
- `batch`: optional
- `age`: optional but recommended

## 4.3 Optional phenotype columns

Supported if present:
- `age`
- `sex`
- `batch`

Behavior:
- `batch` is **used if present**
- `batch` is **ignored if absent**
- text batch labels are better converted to dummy covariates before PLINK association
- the pipeline should never fail just because `batch` is absent

---

# 5. Reference data

## 5.1 Full reusable 1000G reference

Expected path:

```bash
/data/${USER}/gwas_refs/1000g_phase3_allchr/1kg_all
```

Files:
- `1kg_all.bed`
- `1kg_all.bim`
- `1kg_all.fam`

This is the correct ancestry/LD reference for a real full run.

## 5.2 chr22 reference subset for fast testing

Expected path:

```bash
/data/${USER}/gwas_refs/1000g_phase3_allchr/1kg_chr22_test
```

Files:
- `1kg_chr22_test.bed`
- `1kg_chr22_test.bim`
- `1kg_chr22_test.fam`

Use this for fast validation and debugging, not for a full real analysis.

## 5.3 FinnGen base summary statistics

Expected path:

```bash
/data/${USER}/gwas_refs/finngen_r12_lung/finngen_R12_C3_NSCLC_ADENO_EXALLC.gz
```

---

# 6. Building/downloading references

## 6.1 1000G reference build script

The pipeline work used a custom script like:

- `download_1kg_ref.sh`

This script:
- downloads chr1-22 VCFs
- converts to PLINK
- merges chromosomes
- records checksums
- supports resume/skipping
- cleans up or can be extended to clean up unnecessary intermediate files

Important lessons learned during 1000G build:
- PLINK 1.9 cannot use `--set-all-var-ids`
- use `--set-missing-var-ids`
- deduplication options differ by `bcftools` version
- per-chromosome files must be validated before merge
- merge failures can happen from:
  - missing files
  - duplicate IDs
  - multiallelic sites
  - disk-space exhaustion

## 6.2 FinnGen download

The FinnGen file was downloaded separately and used as PRS base input.

---

# 7. Validated test strategy

The best test strategy used here was:

1. take a real chr22 1000G genotype backbone
2. subset individuals
3. rename samples to synthetic IDs like `SIM0001`, `SIM0002`
4. simulate phenotypes on top of real genotype LD structure
5. run the full pipeline end-to-end

Why this is good:
- preserves real LD
- preserves real allele frequencies
- uses real PLINK structure
- avoids needing real study samples during pipeline development
- allows truth-aware validation via causal SNPs

This is better than:
- purely random independent-SNP simulation
- fake phenotypes on top of full 1000G without truth structure

---

# 8. Simulated chr22 validation run

## 8.1 Simulation script

Use:

```bash
scripts/simulate_ld_gwas.sh
```

This script:
- uses `/data/${USER}/gwas_refs/1000g_phase3_allchr/1kg_chr22_test`
- samples a subset of individuals
- renames IDs to `SIM####`
- selects causal SNPs
- simulates balanced case/control phenotypes
- writes:
  - `simulated_ld_data/sim_ld_gwas.{bed,bim,fam}`
  - `simulated_ld_data/sim_phenotypes.csv`
  - `simulated_ld_data/causal_snps.txt`

## 8.2 Run the simulation

```bash
cd /data/rezaa2/gwas-pipeline
source ~/bin/myconda
conda activate gwas_plink19
module load R
bash scripts/simulate_ld_gwas.sh
```

## 8.3 Check simulation outputs

```bash
ls -lh simulated_ld_data
head -5 simulated_ld_data/sim_phenotypes.csv
wc -l simulated_ld_data/sim_ld_gwas.fam
wc -l simulated_ld_data/sim_ld_gwas.bim
wc -l simulated_ld_data/causal_snps.txt
head -5 simulated_ld_data/sim_ld_gwas.fam
```

Expected:
- sample IDs like `SIM0001`
- phenotype file with matching `SIM` IDs
- balanced case/control counts
- causal SNP truth file present

---

# 9. Validated full test run command

Run from repo root:

```bash
nextflow run main.nf \
  -profile biowulflocal,conda \
  --mode full \
  --input_bed /data/rezaa2/gwas-pipeline/simulated_ld_data/sim_ld_gwas \
  --phenotype /data/rezaa2/gwas-pipeline/simulated_ld_data/sim_phenotypes.csv \
  --sumstats /data/${USER}/gwas_refs/finngen_r12_lung/finngen_R12_C3_NSCLC_ADENO_EXALLC.gz \
  --local_1kg_prefix /data/${USER}/gwas_refs/1000g_phase3_allchr/1kg_chr22_test \
  --ancestry_ref_prefix /data/${USER}/gwas_refs/1000g_phase3_allchr/1kg_chr22_test \
  --assoc_test logistic \
  --pop EUR \
  --outdir results_sim_ld_chr22
```

If resuming:

```bash
nextflow run main.nf \
  -profile biowulflocal,conda \
  --mode full \
  --input_bed /data/rezaa2/gwas-pipeline/simulated_ld_data/sim_ld_gwas \
  --phenotype /data/rezaa2/gwas-pipeline/simulated_ld_data/sim_phenotypes.csv \
  --sumstats /data/${USER}/gwas_refs/finngen_r12_lung/finngen_R12_C3_NSCLC_ADENO_EXALLC.gz \
  --local_1kg_prefix /data/${USER}/gwas_refs/1000g_phase3_allchr/1kg_chr22_test \
  --ancestry_ref_prefix /data/${USER}/gwas_refs/1000g_phase3_allchr/1kg_chr22_test \
  --assoc_test logistic \
  --pop EUR \
  --outdir results_sim_ld_chr22 \
  -resume
```

---

# 10. Resource allocation

## 10.1 Recommended VS Code / interactive test session

For chr22 simulation testing:

- CPUs: `8`
- Memory: `32 GB`
- Local scratch: `100 GB`
- Hours: `4`

This was enough for:
- editing
- local Nextflow execution
- simulation generation
- avoiding earlier scratch-space errors

## 10.2 Recommended sbatch chr22 test resources

For a Slurm test run:

```bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=32g
#SBATCH --time=12:00:00
#SBATCH --gres=lscratch:100
```

## 10.3 Full `1kg_all` as target test

This is much heavier and not recommended for routine debugging.

Typical needs:
- CPUs: `16-32`
- Memory: `64-128 GB`
- Local scratch: `600-800 GB`
- Time: `24-48 hours`

## 10.4 Real full-run starting recommendation

Reasonable initial allocation:

- CPUs: `32`
- Memory: `128 GB`
- Time: `48 hours`

Per-process settings should remain moderate.
Do not give every subtask all CPUs.

A good starting pattern:
- `plink`: `8 CPUs`, `32 GB`
- `r_gwas`: `8 CPUs`, `48 GB`
- `prsice`: `8-16 CPUs`, `32-64 GB`

---

# 11. Output structure

For the validated test run:

```bash
results_sim_ld_chr22/
```

Key subdirectories:
- `02_ancestry`
- `03_association`
- `04_plots`
- `05_prs`
- `pipeline_info`

Useful inspection commands:

```bash
find results_sim_ld_chr22 -maxdepth 2 -type f | sort
ls -lh results_sim_ld_chr22/03_association
ls -lh results_sim_ld_chr22/04_plots
ls -lh results_sim_ld_chr22/05_prs
```

---

# 12. What the test run proves

The successful simulated chr22 run proves that the pipeline works end-to-end for:

- PLINK input handling
- QC
- ancestry PCA
- phenotype/covariate handoff
- GWAS association
- association plotting
- PRSice execution
- PRS plotting
- FinnGen summary-stat harmonization

It does **not** prove:
- biological validity of the findings
- transferability of FinnGen PRS to a real lung cancer cohort
- meaningful real-world effect estimation

It is a workflow validation.

---

# 13. What to change for a real full run

To move from the validated simulated chr22 test to a real full run, change these things.

## 13.1 Change the target input cohort

### Test:
```bash
--input_bed /data/rezaa2/gwas-pipeline/simulated_ld_data/sim_ld_gwas
```

### Real:
```bash
--input_bed /path/to/real_cohort
```

Where `/path/to/real_cohort` corresponds to:
- `/path/to/real_cohort.bed`
- `/path/to/real_cohort.bim`
- `/path/to/real_cohort.fam`

## 13.2 Change the phenotype file

### Test:
```bash
--phenotype /data/rezaa2/gwas-pipeline/simulated_ld_data/sim_phenotypes.csv
```

### Real:
```bash
--phenotype /path/to/real_pheno.csv
```

## 13.3 Change ancestry / LD reference from chr22 test to full 1kg_all

### Test:
```bash
--local_1kg_prefix /data/${USER}/gwas_refs/1000g_phase3_allchr/1kg_chr22_test
--ancestry_ref_prefix /data/${USER}/gwas_refs/1000g_phase3_allchr/1kg_chr22_test
```

### Real:
```bash
--local_1kg_prefix /data/${USER}/gwas_refs/1000g_phase3_allchr/1kg_all
--ancestry_ref_prefix /data/${USER}/gwas_refs/1000g_phase3_allchr/1kg_all
```

## 13.4 Keep `--sumstats` as the chosen PRS base GWAS

Example:
```bash
--sumstats /path/to/base_sumstats.gz
```

## 13.5 Change execution profile

### Test:
```bash
-profile biowulflocal,conda
```

### Real:
```bash
-profile biowulf,conda
```

Use local profile for small debugging only.
Use Slurm profile for real full runs.

---

# 14. Real full-run command template

```bash
nextflow run main.nf \
  -profile biowulf,conda \
  --mode full \
  --input_bed /path/to/real_cohort \
  --phenotype /path/to/real_pheno.csv \
  --sumstats /path/to/base_sumstats.gz \
  --local_1kg_prefix /data/${USER}/gwas_refs/1000g_phase3_allchr/1kg_all \
  --ancestry_ref_prefix /data/${USER}/gwas_refs/1000g_phase3_allchr/1kg_all \
  --assoc_test logistic \
  --pop EUR \
  --outdir results_real_full
```

---

# 15. Handling `batch` in a real full run

## 15.1 Is `batch` necessary?

No. `batch` is not inherently required.

Use it only if it is a real technical covariate, such as:
- array version
- lab batch
- processing wave
- recruitment center
- technical study subgroup

## 15.2 Desired behavior

The pipeline should:
- include `batch` if present
- ignore it if absent
- never fail just because `batch` is missing

## 15.3 Recommended implementation

In `scripts/run_association.R`:
- always include:
  - `FID`
  - `IID`
  - available PCs
  - `age` if present
  - `sex` if present
- include `batch` only if present
- if `batch` is text, convert to numeric dummy variables before passing to PLINK
- if absent, do nothing

For the test run, `batch` was dropped from association covariates because text batch values caused PLINK to treat all covariates as missing.

---

# 16. Preparing the pipeline for real sample input

Before running a real cohort, do these checks.

## 16.1 Check genotype files exist

```bash
ls -lh /path/to/real_cohort.bed /path/to/real_cohort.bim /path/to/real_cohort.fam
```

## 16.2 Check phenotype file

```bash
head -5 /path/to/real_pheno.csv
```

Verify:
- `FID`
- `IID`
- `PHE`
- `PHE` values are 1/2
- optional `age`, `sex`, `batch` look reasonable

## 16.3 Check ancestry reference

```bash
ls -lh /data/${USER}/gwas_refs/1000g_phase3_allchr/1kg_all.bed
ls -lh /data/${USER}/gwas_refs/1000g_phase3_allchr/1kg_all.bim
ls -lh /data/${USER}/gwas_refs/1000g_phase3_allchr/1kg_all.fam
```

## 16.4 Check base sumstats

```bash
zcat -f /path/to/base_sumstats.gz | head
```

## 16.5 Check phenotype / FAM overlap

At minimum, confirm study IDs in the phenotype file overlap IDs in `.fam`.

---

# 17. Optional simulation validation

For simulated-data runs, it is useful to compare:

- causal SNP truth file
- recovered GWAS hits
- PRS separation between cases and controls

## 17.1 Add optional parameter

In `nextflow.config`:

```groovy
params.causal_snps = null
```

## 17.2 Expected simulation truth file

Example:
```bash
simulated_ld_data/causal_snps.txt
```

## 17.3 Validation goals

### GWAS recovery
Compare:
- causal SNP list
- association results

Metrics:
- number of causal SNPs present in GWAS results
- number passing:
  - `p < 5e-8`
  - `p < 1e-5`
  - `p < 0.05`

### PRS separation
Compare:
- `prs.best`
- phenotype file

Metrics:
- mean PRS in cases vs controls
- logistic model pseudo-R²
- quartile enrichment
- case fraction by PRS quartile

## 17.4 Recommended optional files

A future validation module can write:
- `causal_recovery.txt`
- `prs_validation.txt`
- `prs_distribution.png`
- `prs_quartile_summary.txt`

---

# 18. Important implementation lessons from debugging

These issues were encountered and fixed during development.

## 18.1 PLINK 1.9 output filenames
PLINK 1.9 writes:
- `*.assoc.logistic`
- `*.assoc.linear`

not:
- `*.logistic`
- `*.linear`

## 18.2 Missing-sex handling
If `.fam` sex is `0`, PLINK may ignore phenotypes unless:

```bash
--allow-no-sex
```

is included.

## 18.3 Phenotype file for PLINK
`plink_pheno.txt` must be:
- no header
- no quotes
- whitespace/tab delimited
- `FID IID PHE`

## 18.4 `glm(..., family=binomial())` in R
R expects `0/1` or a two-level factor, not PLINK `1/2`.

So in PRS summary scripts:
- recode `1 -> 0`
- recode `2 -> 1`

## 18.5 PRS summary option parsing
The PRS summary script must parse named options correctly with `optparse`, otherwise arguments like `--best` can be misread as literal file paths.

## 18.6 Simulation IDs
Simulated target samples must **not** keep original 1000G IDs when ancestry PCA uses 1000G as reference.
Rename to `SIM####` to avoid confusing target and reference.

## 18.7 Batch covariates
Text `batch` columns can break PLINK association if passed directly.
They should be:
- omitted
- or converted to numeric dummy variables

---

# 19. Common errors and fixes

## `nextflow: command not found`
Load module:

```bash
module load nextflow
```

## `Rscript: command not found`
Load R:

```bash
module load R
```

## `plink: command not found`
Activate conda:

```bash
source ~/bin/myconda
conda activate gwas_plink19
```

## `No space left on device`
Increase local scratch.
For local chr22 testing, `100 GB` was a good working value.

## `Process requirement exceeds available CPUs`
Lower per-process CPU requests in `nextflow.config`, especially `prsice`.

## `y values must be 0 <= y <= 1`
Recode `PHE` from PLINK `1/2` to `0/1` before binomial GLM in R scripts.

## `Skipping --logistic since less than two phenotypes are present`
Usually caused by:
- phenotype file mismatch
- invalid covariates
- missing-sex phenotype suppression
- no overlapping IDs

## `gwas_raw.logistic missing`
Use the correct PLINK 1.9 filename:
- `gwas_raw.assoc.logistic`

## `File '--best' does not exist`
Fix option parsing in `summarise_prs.R` with `optparse`.

---

# 20. Git and .gitignore

Recommended `.gitignore` should include:

- `.nextflow/`
- `work/`
- `.nextflow.log*`
- `results/`
- `results_*/`
- `prs_results/`
- `simulated_ld_data/`
- `*.out`
- `*.err`
- `*.zip`
- generated PLINK files if produced inside repo

Keep tracked:
- `main.nf`
- `nextflow.config`
- `conf/`
- `modules/`
- `scripts/`
- `run_*.sbatch`
- small example files in `data/`

---

# 21. Recommended workflow for future development

## Fast debugging
Use:
- `1kg_chr22_test`
- simulated phenotype data
- `-profile biowulflocal,conda`

## Branch-specific testing
To isolate failures:
- `--mode assoc`
- `--mode prs`

## Full real analysis
Use:
- real cohort PLINK files
- real phenotype CSV
- `1kg_all`
- `-profile biowulf,conda`

---

# 22. Example commands

## 22.1 Simulate chr22 data

```bash
cd /data/rezaa2/gwas-pipeline
source ~/bin/myconda
conda activate gwas_plink19
module load R
bash scripts/simulate_ld_gwas.sh
```

## 22.2 Full validated test run

```bash
nextflow run main.nf \
  -profile biowulflocal,conda \
  --mode full \
  --input_bed /data/rezaa2/gwas-pipeline/simulated_ld_data/sim_ld_gwas \
  --phenotype /data/rezaa2/gwas-pipeline/simulated_ld_data/sim_phenotypes.csv \
  --sumstats /data/${USER}/gwas_refs/finngen_r12_lung/finngen_R12_C3_NSCLC_ADENO_EXALLC.gz \
  --local_1kg_prefix /data/${USER}/gwas_refs/1000g_phase3_allchr/1kg_chr22_test \
  --ancestry_ref_prefix /data/${USER}/gwas_refs/1000g_phase3_allchr/1kg_chr22_test \
  --assoc_test logistic \
  --pop EUR \
  --outdir results_sim_ld_chr22
```

## 22.3 Resume a partial run

```bash
nextflow run main.nf \
  -profile biowulflocal,conda \
  --mode full \
  --input_bed /data/rezaa2/gwas-pipeline/simulated_ld_data/sim_ld_gwas \
  --phenotype /data/rezaa2/gwas-pipeline/simulated_ld_data/sim_phenotypes.csv \
  --sumstats /data/${USER}/gwas_refs/finngen_r12_lung/finngen_R12_C3_NSCLC_ADENO_EXALLC.gz \
  --local_1kg_prefix /data/${USER}/gwas_refs/1000g_phase3_allchr/1kg_chr22_test \
  --ancestry_ref_prefix /data/${USER}/gwas_refs/1000g_phase3_allchr/1kg_chr22_test \
  --assoc_test logistic \
  --pop EUR \
  --outdir results_sim_ld_chr22 \
  -resume
```

## 22.4 Real full run

```bash
nextflow run main.nf \
  -profile biowulf,conda \
  --mode full \
  --input_bed /path/to/real_cohort \
  --phenotype /path/to/real_pheno.csv \
  --sumstats /path/to/base_sumstats.gz \
  --local_1kg_prefix /data/${USER}/gwas_refs/1000g_phase3_allchr/1kg_all \
  --ancestry_ref_prefix /data/${USER}/gwas_refs/1000g_phase3_allchr/1kg_all \
  --assoc_test logistic \
  --pop EUR \
  --outdir results_real_full
```

---

# 23. Final note

This repository is now in a good state for:

- reproducible simulated validation
- continued pipeline development
- transition to real study cohorts

The next major step is to plug in real sample input and perform a full run using:
- a real PLINK cohort
- a real phenotype file
- the full local `1kg_all` reference
- an appropriate base GWAS summary-stat file
