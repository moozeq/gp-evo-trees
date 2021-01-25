# evopipe

A phylogenetic pipeline for inferring a species/genome tree from a set of
genomes by clustering, inferring gene families and their trees.

## Requirements

- [RAxML](https://github.com/stamatak/standard-RAxML) (must be in `$PATH` as `raxml`)
- [ninja](http://nimbletwist.com/software/ninja/index.html) (must be in `$PATH` as `ninja`)
- [clann](https://github.com/ChrisCreevey/clann) (must be in `$PATH` as `clann`)
- [muscle](https://anaconda.org/bioconda/muscle) (`conda install -c bioconda muscle`)
- [mmseqs2](https://github.com/soedinglab/MMseqs2) (`conda install -c conda-forge -c bioconda mmseqs2`)
- [biopython](https://biopython.org/) (`conda install -c conda-forge biopython`)
- [ete3](http://etetoolkit.org/) (`conda install -c etetoolkit ete3`)
- [joblib](https://joblib.readthedocs.io/) (`conda install -c anaconda joblib`)
- [requests](https://requests.readthedocs.io/en/master/) (`conda install -c anaconda requests`)

## Run

```bash
usage: pipe.py [-h] [-n NUM] [--cluster-min CLUSTER_MIN] [--cluster-highest CLUSTER_HIGHEST] [--filter-min FILTER_MIN] [--filter-max FILTER_MAX] [--fastas-dir FASTAS_DIR] [-l LOG] [-d] [-o OUTPUT]
               {family,file} input

Phylogenetic pipeline to infer a species/genome tree from a set of genomes

positional arguments:
  {family,file}         pipeline mode
  input                 family name or .json file with species names list which will be inferred

optional arguments:
  -h, --help            show this help message and exit
  -n NUM, --num NUM     limit downloading species to specific number
  --cluster-min CLUSTER_MIN
                        filter cluster proteomes minimum, by default: 4
  --cluster-highest CLUSTER_HIGHEST
                        get only "n" most populated clusters
  --filter-min FILTER_MIN
                        filter proteomes minimum
  --filter-max FILTER_MAX
                        filter proteomes maximum
  --fastas-dir FASTAS_DIR
                        directory name with fasta files, by default: "fastas/"
  -l LOG, --log LOG     logger file
  -d, --duplications    allow duplications (paralogs)
  -o OUTPUT, --output OUTPUT
                        output directory, by default: name of family if "family" mode, otherwise "results"
```

## Example

### File with species

Providing `.json` file with list of species to be inferred:
```bash
$ head species.json
[
    "SARS coronavirus civet020",
    "Bat coronavirus",
    "Dromedary camel coronavirus HKU23",
    "Hipposideros bat coronavirus HKU10",
    "Human betacoronavirus 2c Jordan-N3/2012",
    "Murine coronavirus SA59/RJHM",
    "Feline coronavirus UU20",
    "SARS coronavirus Sino3-11",
    "Bat SARS-like coronavirus YNLF_34C",
```

Run pipeline with:
```bash
./pipe.py file species.json -o coronaviruses
```

All proteomes will be stored under `fastas` directory. All trees will be
available at `coronaviruses/trees`, i.e.: `coronaviruses/trees/nj_super_tree_species.nwk`

### Family name

Providing family name, proteomes (each for one organism) will be downloaded
(sorted by a score - so from best to worst) and then used to build trees:

Run pipeline with:
```bash
./pipe.py family Coronaviridae -n 100
```

With above command, we'll obtain maximum 100 proteomes from _Coronaviridae_ which
will be then stored under `fastas` directory. All trees will be available at `Coronaviridae/trees`, i.e.:
`Coronaviridae/trees/nj_super_tree_species.nwk`

## Pipeline

Following steps are performed:
1) Download proteomes for provided tax family, or for species names from `.json` file
2) Filter fasta files with `min` and `max` sequences within, specified with `args`
3) Change sequences IDs to inner species IDs for building trees purposes
4) Merge all proteomes into one, big fasta file
5) Cluster merged fasta file to obtain protein families from all species and filter
   them getting only first `n` most populated clusters; with `min` sequences; with or
   without duplications (`dup`) - all options specified in `args`
6) Save clusters to separate files, each per protein family
7) Align all sequences within each protein family fasta file
8) Build NJ, ML and MP trees from provided aligned protein families fasta files
9) Retrieve species names for consensus and super trees from their inner IDs