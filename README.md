# vivaxReC

## What is vivaxReC ?
**vivaxReC** (**v**ivax **R**ecurrence **C**onsensus) is a consensus framework for processing MicroHaplotype data (genotyping data following protocols and design workflows created by **VivaxGEN**) to use towards classifying *Plasmodium vivax* recurrent infections - this will be done by integrating outputs from established genotyping-based recurrence tools (e.g. *Pv3Rs*, with *Aster* being added at a later date), allele-sharing, relatedness estimates (including *Dcifer* and background relatedness).

**vivaxReC** does not introduce a new recurrence classification algorithm, nor does it reinvent existing methods. All primary methodological credit remains with the original tool authors.

## Methods and attribution

**vivaxReC** integrates outputs from several existing analysis tools. Users of this repository should cite the original methods accordingly:

- **Pv3Rs** – probabilistic classification of recrudescence, relapse, and reinfection using within-host genotype series. Available at: https://github.com/aimeertaylor/Pv3Rs
- **Dcifer** – identity-by-descent–based relatedness estimation across loci. Available at: https://github.com/EPPIcenter/dcifer

Coming soon:
- **Aster** – likelihood-based recurrence inference incorporating complexity of infection, detection probability, allele frequencies, and background relatedness. Available at: https://github.com/EPPIcenter/asterTES

**vivaxReC** does not modify or reimplement the underlying algorithms from these tools. All likelihoods, relatedness estimates, and probabilistic classifications are computed using the original software, with **vivaxReC** providing post hoc integration, quality control, and consensus-calling logic.

Please refer to the original publications and documentation for for full methodological details and citation information.


## Pipeline overview - ADDING MORE QC AND DOWNSTREAM ANALYSES

```
patient_infection_sequence.tsv  ──┐
mhap_filtered.tsv               ──┤  Step 1: vivaxrec prep
                                   ↓
                        master_samples.tsv
                        per_sample_locus_alleles.tsv
                        d0_dr_pairs.tsv
                        allele_freq_from_D0.tsv
                                   ↓
baseline_pairs_with_relatedness.tsv  ──┐
consecutive_pairs_with_relatedness.tsv ┤  Step 2: vivaxrec classify
                                       ↓
                  per_infection_classifications.tsv  (baseline)
                  per_infection_classifications.tsv  (consecutive)
                                   ↓
                               Step 3: vivaxrec plot
                                   ↓
                          figures/  (PNG outputs)
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
usage: vivaxrec [-h] {prepare-input-data,calculate-priors,run-consensus-pipeline,run-aster-dcifer,plot-outputs} ...
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

*Note - vivaxReC should work with other forms of multi-allelic sequencing data but accuracy will vary depending on the number of sequencing markers used.*
