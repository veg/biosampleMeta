# biosampleMeta

Some biosample metadata scraping efforts from the Pond lab for the ArgosDB project.

## Installation and setup

Recommend [miniconda](https://docs.conda.io/en/latest/miniconda.html) and [mamba](https://github.com/mamba-org/mamba). Choose your environment (currently `m1`), then

```
mamba env create -f environment-$ENVIRONMENT.yml
conda activate fdaargos
```

It is also recommended that you set environment variables `ENTREZ_EMAIL` and `ENTREZ_API_KEY` via

```
export ENTREZ_EMAIL=$ENTREZ_EMAIL
export ENTREZ_API_KEY=$ENTREZ_API_KEY
```

## Usage

Fetch all biosamples from a given bioproject accession:

```
snakemake -j 1 bioprojects/$BP_ACCESSION/biosample_links.txt
snakemake -j 1 bioprojects/$BP_ACCESSION/biosamples/all.txt
```

For instance, try [PRJNA231221](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA231221), the original ArgosDB bioproject.
