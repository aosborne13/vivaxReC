# vivaxReC

## What is vivaxReC ?
**vivaxReC** (**v**ivax **R**ecurrence **C**onsensus) is a consensus framework for classifying *Plasmodium vivax* recurrent infections by integrating outputs from established genotyping-based recurrence tools (e.g. *Pv3Rs* and *Aster*), allele-sharing, relatedness estimates (including *Dcifer* and *Aster*-derived background relatedness), and time-to-event information.

**vivaxReC** does not introduce a new recurrence classification algorithm, nor does it reimplement existing methods. Instead, it acts as an adjudication and evidence-integration layer, synthesizing complementary probabilistic, genetic,
and temporal signals to generate robust, biologically informed, recurrence classifications, while explicitly retaining uncertainty in ambiguous cases. All primary methodological credit remains with the original tool authors.

## Methods and attribution

**vivaxReC** integrates outputs from several existing recurrence and relatedness inference tools. Users of this repository should cite the original methods accordingly:

- **Pv3Rs** – probabilistic classification of recrudescence, relapse, and reinfection using within-host genotype series. Available at: https://github.com/aimeertaylor/Pv3Rs
- **Aster** – likelihood-based recurrence inference incorporating complexity of infection, detection probability, allele frequencies, and background relatedness. Available at: https://github.com/EPPIcenter/asterTES
- **Dcifer** – identity-by-descent–based relatedness estimation across loci. Available at: https://github.com/EPPIcenter/dcifer

**vivaxReC** does not modify or reimplement the underlying algorithms from these tools. All likelihoods, relatedness estimates, and probabilistic classifications are computed using the original software, with **vivaxReC** providing post hoc integration, quality control, and consensus-calling logic.

Please refer to the original publications and documentation for *Pv3Rs*, *Aster*, and *Dcifer* for full methodological details and citation information.

## Installation
Clone GitHub repo in appropriate location on device or HPC (e.g. ~/tools/ folder).
```
git clone https://github.com/YOUR_GITHUB_USERNAME/vivaxReC.git
cd vivaxReC
```
Create conda environment using the ```.yml``` file in the GitHub repo install.
*Note - is conda is active, first run* ```conda deactivate```
```
mamba env create -f env/environment.yml
```
Activate environment and install remaining R-based software packages (inclides *Pv3Rs*,*asterTES*, and *Dcifer*).
```
conda activate vivaxReC
Rscript env/install_R_packages.R
```
Finally, install vivaxReC CLI and test install.
```
pip install -e .
vivaxrec --help
```
If installed correctly, you should see:
```
usage: vivaxrec [-h] {prepare-input-data,calculate-priors,run-consensus-pipeline,run-aster-dcifer,plot-outputs} ...
```

## Input data
**vivaxReC** was designed to analyse *P. vivax* microhaplotpe genotyping data following protocols and design workflows created by **VivaxGEN**.

- Sequencing assay and panel design methodology - https://www.nature.com/articles/s41467-025-62357-x
- Bioinformatic processing pipeline -  https://github.com/vivaxgen/MicroHaps

*Note - vivaxReC should work with other forms of multi-allelic sequencing data but accuracy will vary depending on the number of sequencing markers used.*

### Formatting input data to use vivaxReC pipeline

