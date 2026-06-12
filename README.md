# vivaxReC

## What is vivaxReC ?
**vivaxReC** (**v**ivax **R**ecurrence **C**onsensus) is a consensus framework for processing MicroHaplotype data (genotyping data following protocols and design workflows created by **VivaxGEN**) to use towards classifying *Plasmodium vivax* recurrent infections - this will be done by integrating outputs from established genotyping-based recurrence tools (e.g. *Pv3Rs*, with *Aster* being added at a later date), allele-sharing, relatedness estimates (including *Dcifer* and background relatedness).

**vivaxReC** does not introduce a new recurrence classification algorithm, nor does it reinvent existing methods. All primary methodological credit remains with the original tool authors.

## Methods and attribution

**vivaxReC** integrates outputs from several existing analysis tools. Users of this repository should cite the original methods accordingly:

- **Pv3Rs:** [Pv3Rs GitHub](https://github.com/aimeertaylor/Pv3Rs) – probabilistic classification of recrudescence, relapse, and reinfection using within-host genotype series.
- **DCiFER:** [dcifer GitHub](https://github.com/EPPIcenter/dcifer) – identity-by-descent–based relatedness estimation across loci.
- **MOIRE:** [moire GitHub](https://github.com/EPPIcenter/moire) – an MCMC based approach to estimating complexity of infection (COI), also sometimes referred to as multiplicity of infection (MOI), population allele frequencies, and within-host relatedness from polyallelic genomics data.

Coming soon:
- **Aster** [asterTES GitHub](https://github.com/EPPIcenter/asterTES)  – likelihood-based recurrence inference incorporating complexity of infection, detection probability, allele frequencies, and background relatedness.

**vivaxReC** does not modify or reimplement the underlying algorithms from these tools. All likelihoods, relatedness estimates, and probabilistic classifications are computed using the original software, with **vivaxReC** providing post hoc integration, quality control, and consensus-calling logic.

Please refer to the original publications and documentation for for full methodological details and citation information.


## Pipeline overview - ADDING MORE QC AND DOWNSTREAM ANALYSES

```
  metadata.tsv                      ──┐
  mhap_effective_COI_summary.tsv    ──┤  Step 1: vivaxrec standardize
  mhap_between_relatedness_estimate.tsv ┘
                                        ↓
                           standardized_inputs/
                             analysis_ready/
                               metadata.cleaned.tsv
                               coi.cleaned.tsv
                               relatedness.cleaned.tsv
                             dropped/   (samples not in metadata)
                             qc/        (duplicate IDs, missing dates, etc.)
                                        ↓
                                   Step 2: vivaxrec build
                                        ↓
                           patient_level_tables/
                             patient_infection_sequence.tsv
                             baseline_pairs_with_relatedness.tsv
                             consecutive_pairs_with_relatedness.tsv
                             baseline_samples.tsv
                             summary_by_trial_site.tsv
                                        ↓
                                   Step 3: vivaxrec prep
                                        ↓
                           recurrence-classification/<trial>/
                             master_samples.tsv
                             per_sample_locus_alleles.tsv
                             d0_dr_pairs.tsv
                             allele_freq_from_D0.tsv
                                        ↓
                                   Step 4: vivaxrec classify
                                        ↓
                           recurrence-classification/<trial>/
                             baseline-to-recurrence/
                               per_infection_classifications.tsv
                             consecutive-pairs/
                               per_infection_classifications.tsv
                                        ↓
                                   Step 5: vivaxrec plot
                                        ↓
                           recurrence-classification/<trial>/
                             figures/
```

## Installation
### Step 1: Install micromamba (if you don't already have it).
```
curl -L micro.mamba.pm/install.sh | bash
exec bash  # restart your shell so micromamba is available
```
### Step 2: Clone GitHub repo in appropriate location on device or HPC, such as a tools directory.
```
git clone https://github.com/aosborne13/vivaxReC
cd vivaxReC
```
### Step 3: Create *vivaxReC* envrionment and install all dependencies.
*Note - if conda is also installed on your device, run* ```conda deactivate``` *before proceeding with this next step*
```
bash scripts/setup_env.sh
```
### Step 4: Test the installation.
```
vivaxrec --help
```
If installed correctly, you should see:
```
usage: vivaxrec [-h] {vivaxrec prep,vivaxrec classify,vivaxrec plot,vivaxrec run-all} ...
```

### Updating
To update, activate the *vivaxReC* enironment and then run the following:
```
cd vivaxReC
git pull
micromamba activate vivaxReC
python -m pip install -e .
```
If only environment dependencies have changed:
```
micromamba activate vivaxReC
micromamba update -f env/environment.yml
Rscript env/install_R_packages.R
```
## Input data
**vivaxReC** was designed to analyse *P. vivax* microhaplotpe genotyping data following protocols and design workflows created by **VivaxGEN**.

- Sequencing assay and panel design methodology - https://www.nature.com/articles/s41467-025-62357-x
- Bioinformatic processing pipeline -  https://github.com/vivaxgen/MicroHaps

*Note - vivaxReC should work with other forms of multi-allelic sequencing data but accuracy will vary depending on the number of sequencing markers used. **Do not use WGS data**.*

The pipeline takes three raw input files — the `trial` column in the metadata file determines how samples are grouped, and any trial name(s) can be used. Multiple trials, or categories, can coexist in the same files.

### `metadata.tsv`

One row per sample. Required columns:

| Column | Description |
|--------|-------------|
| `sample_id` | Unique sample identifier |
| `external_id` | Individual identifier — authoritative grouping key |
| `trial` | Trial name (any string) |
| `treatment_arm` | Treatment arm label |
| `site` | Study site |
| `day` | Day relative to enrolment |
| `time_point` | Time point label (e.g. `D0` for Day 0, `DR1` for first recurrence,`DR2` for second recurrence, etc.) |
| `samp_date` | Sample collection date (dd/mm/yyyy) |

### Output from VivaxGEN pipeline

```
outputCIGAR.tsv
		↓ (prepare_data.r --infile outputCIGAR.tsv --rmindel --outdir analysis-filtered-no-indel)
	mhap_filtered.tsv
	↓				↓
MOIRE				Dcifer
(run_Moire.r)			(run_Dcifer.r)
↓									↓
mhap_effective_COI_summary.tsv		mhap_between_relatedness_estimate.tsv
```

#### `mhap_effective_COI_summary.tsv`

MOIRE complexity of infection estimates. Main columns: `sample_id`, `post_effective_coi_mean`. Additional columns: `post_effective_coi_lower`, `post_effective_coi_med`, `post_effective_coi_upper`.

#### `mhap_between_relatedness_estimate.tsv`

Pairwise DCiFER relatedness estimates. Main columns: `sample_id1`, `sample_id2`, `relatedness`.

---

## Usage

### Step 1 — Standardise inputs

Cleans sample IDs, filters to the metadata universe, and produces QC reports.

```
vivaxrec standardize \
  --rel  /path/to/mhap_between_relatedness_estimate.tsv \
  --coi  /path/to/mhap_effective_COI_summary.tsv \
  --meta /path/to/metadata.tsv \
  --out_dir /path/to/standardized_inputs
```

Outputs:

- `standardized_inputs/analysis_ready/` — cleaned metadata, COI, and relatedness files restricted to samples present in metadata
- `standardized_inputs/dropped/` — rows excluded because the sample was not in metadata
- `standardized_inputs/qc/` — duplicate sample IDs, missing dates, day/time_point mismatches

### Step 2 — Build patient-level tables

Constructs infection sequences and recurrence pairs by individual, merges in relatedness estimates.

```
vivaxrec build \
  --analysis_dir /path/to/standardized_inputs/analysis_ready \
  --out_dir /path/to/patient_level_tables
```

Outputs written to `patient_level_tables/`:

- `patient_infection_sequence.tsv` — one row per sample, with infection index and COI
- `baseline_pairs_with_relatedness.tsv` — D0-to-each-recurrence pairs with DCiFER relatedness
- `consecutive_pairs_with_relatedness.tsv` — consecutive recurrence pairs
- `baseline_samples.tsv` — baseline (D0) samples only
- `summary_by_trial_site.tsv` — patient and pair counts by trial and site

### Step 3 — Prepare Pv3Rs inputs

Formats the patient-level tables into the structures required by the classification step.

```
vivaxrec prep \
  --pis  /path/to/patient_level_tables/patient_infection_sequence.tsv \
  --geno /path/to/mhap_filtered.tsv \
  --out_dir /path/to/recurrence-classification
```

Outputs are written to `<out_dir>/` and `<out_dir>/<trial>/` for each trial found in the data.

### Step 4 — Classify recurrences

Run once per trial and pair mode. The `trial` flag must match a value in the `trial` column of `patient_infection_sequence.tsv` (case-insensitive).

```
BASE=/path/to/your/analysis
TRIAL=mytrial

# Baseline pairs (D0 → recurrence)
vivaxrec classify \
  trial $TRIAL pair_mode baseline max_recrud_days 45 \
  input_dir  $BASE/patient_level_tables \
  geno       $BASE/mhap_filtered.tsv \
  master_dir $BASE/recurrence-classification \
  out_dir    $BASE/recurrence-classification/$TRIAL/baseline-to-recurrence

# Consecutive pairs (recurrence → next recurrence)
vivaxrec classify \
  trial $TRIAL pair_mode consecutive max_recrud_days 45 \
  input_dir  $BASE/patient_level_tables \
  geno       $BASE/mhap_filtered.tsv \
  master_dir $BASE/recurrence-classification \
  out_dir    $BASE/recurrence-classification/$TRIAL/consecutive-pairs
```

Repeat for each trial in your dataset.

**Option flags - modify depending on modelling parameters logical to your dataset or analysis requirements:**

| Flag | Default | Description |
|------|---------|-------------|
| `max_recrud_days` | `45` | Interval threshold for post-hoc C→L/I reassignment |
| `af_model` | `mix` | Allele frequency model: `d0`, `pair`, or `mix` |
| `window_days` | `35 120` | Time window boundaries (days) |
| `pC` | `0.08`, `0.01`, `0.005` | Per-window recrudescence prior |
| `pL` | `0.25`, `0.70`, `0.40` | Per-window relapse prior (non-PQ arms) |
| `lambda_foi_w` | `0.002` | Daily reinfection rate (force of infection) |
| `pq_discount` | `0.19` | Relapse prior multiplier for PQ arms in windows 2–3 |
| `pq_arms` | *(none)* | Treatment arm names to receive PQ discount (e.g. ` "7-day"`, `"14-day"`, `ALPQ`, `CQPQ`) |
| `first_rec_only` | `off` | `on` to additionally re-run Pv3Rs restricted to D0 + first recurrence per patient |
| `mono_only` | `off` | `on` to additionally re-run Pv3Rs restricted to patients with monoclonal D0 (raw eMOI ≤ threshold) |
| `mono_emoi_th` | `1.0` | eMOI threshold for monoclonal D0 definition (raw `post_effective_coi_mean`, not rounded) |

The default prior values are derived from published *P. vivax* radical cure trial data. Users working in different transmission settings should review and adjust these to reflect local epidemiology and treatment efficacy. See the [Prior system](#prior-system) section below.

> **Note on Recrudesence:** Pv3Rs does not include a time-to-event model and may assign recrudesence, p(C), to individuals with a recurrence that is an exact clone (Dcifer = 1; scale 0 to 1.0) of a previous infection, regardless of time point. High allele sharing between an infection episode and a subsequent recurrent infection greater than 45 days after the previous infection is most parsimoniously attributable to a hypnozoite-derived relapse over a recrudescence from a persisting blood-stage clone.

> **Note on subset analyses:** `first_rec_only` and `mono_only` each trigger a completely independent Pv3Rs run on a restricted dataset — they do not just filter the full output. This is essential because Pv3Rs uses all presented episodes as joint context when computing posteriors; subsetting the output after a full run would give *incorrect* results.

> **Note on pairs if the total genotype count exceeds 8:** When the total number of genotypes within an individual exceeds 8, Pv3Rs will fail and produce an *NA* value. To circumvent this, pairs within an individual are subset and run independently to produce posterior estimates. This must be acknowledged as a limitation of full-picture results when analysing outputs. When the number of genotypes still exceeds 8 within individual pairs, the output will remain *NA*.

### Step 5 — Plot

```
vivaxrec plot \
  trial $TRIAL pair_mode both \
  bl_dir   $BASE/recurrence-classification/$TRIAL/baseline-to-recurrence \
  con_dir  $BASE/recurrence-classification/$TRIAL/consecutive-pairs \
  figs_dir $BASE/recurrence-classification/$TRIAL/figures \
  pis      $BASE/patient_level_tables/patient_infection_sequence.tsv
```

### Full pipeline shortcut

Run all five steps in sequence with `run-all`. Pass `--trials` to specify which trial(s) to process through steps 4 and 5; all other flags are forwarded to the relevant step.

```
vivaxrec run-all \
  --trials mytrial \
  --rel     /data/mhap_between_relatedness_estimate.tsv \
  --coi     /data/mhap_effective_COI_summary.tsv \
  --meta    /data/metadata.tsv \
  --out_dir /data/standardized_inputs
```

For multiple trials:

```
vivaxrec run-all --trials trial1 trial2 trial3 ...
```

## Output files

### `per_infection_classifications.tsv`

The primary output — one row per classified pair, all patients included. When `first_rec_only on` or `mono_only on` are passed, additional files are written to the same directory:

| File | Description |
|------|-------------|
| `per_infection_classifications.tsv` | All patients, all recurrences |
| `per_infection_classifications_first_rec_only.tsv` | D0 + first recurrence only; Pv3Rs re-run on this two-episode context |
| `per_infection_classifications_mono_only.tsv` | Patients with monoclonal D0 (eMOI ≤ threshold) only; Pv3Rs re-run on subset |

Each file has a corresponding `sanity_check.txt` / `sanity_check_first_rec_only.txt` / `sanity_check_mono_only.txt` clearly documenting which subset was used and how many patients and pairs were included.

Key columns in every classifications file:

| Column | Description |
|--------|-------------|
| `id` | Patient identifier (`external_id`) |
| `treatment_arm` | Treatment arm |
| `site` | Study site |
| `dcifer_r` | DCiFER pairwise IBD estimate |
| `pv3rs_C` / `pv3rs_L` / `pv3rs_I` | Pv3Rs posterior probabilities for C, L, I |
| `pv3rs_top` | Top Pv3Rs call (post recrudescence reassignment) |
| `pv3rs_reassigned` | `TRUE` if recrudescence call was time-reassigned |
| `share_frac` | Allele-sharing fraction |
| `pv3rs_note` | Notes (e.g. `stepwise_fallback`) |

### `sanity_check.txt`

Summary of key parameters and data quality metrics for the classification run.

### `summary_by_trial_site.tsv`

Patient and pair counts by trial and site, produced by the `build` step.

---

## Prior system

The pipeline uses a biologically-informed three-window time prior. Defaults are derived from published *P. vivax* radical cure trial data and should be reviewed for each new study setting. These can be adjusted by any user for their specific datasets.

| Window | Days | Description |
|--------|------|-------------|
| W1 | ≤ 35 | Recrudescence possible; short-latency relapse possible |
| W2 | 36–120 | Recrudescence very unlikely; primary relapse window |
| W3 | > 120 | Recrudescence implausible; reinfection accumulates |

Window boundaries and prior values are all configurable via CLI flags. If primaquine (PQ) treatment arms are present, the relapse prior in Windows 2 and 3 is discounted by `pq_discount` (default 0.19, reflecting ~5-fold reduction in relapse from published trial data). Adjust this multiplier to match the efficacy observed in your study.

The `lambda_foi_w = 0.002/day` default reinfection rate is conservative; relapse proportions derived under this prior should be treated as lower bounds in high-transmission settings.
