# evopipe

A phylogenetic pipeline for inferring a species/genome tree from a set of
genomes by clustering, inferring gene families and their trees.

It can be used for inferring all selected species, below are results obtained
for _Coronaviridae_ viruses family:

![Coronaviridae](demo/ML_supertree_organisms.png)

Report for _Coronaviridae_ project is available [here](docs/report_EN.pdf).

## Run (Docker)

```basha
# build
docker build -t evopipe github.com/moozeq/gp-evo-trees#pipeline

# run for species from file analysis (see example below)
docker run --name evopipe_file \
           -v $PWD/results:/evopipe/results \
           -v $PWD/species.json:/evopipe/species.json \
           -t evopipe file species.json

# run for family analysis (see example below)
docker run --name evopipe_family \
           -v $PWD/Coronaviridae:/evopipe/Coronaviridae \
           -t evopipe family Coronaviridae
```

All analyzed data will be stored under mounted point:
- `$PWD/results` - default directory for results when read species from file
- `$PWD/Coronaviridae` - directory corresponds to family name if specified instead of file

If output directory specified as argument, i.e. `-o coronaviruses`, then mounting
point must be the same - `-v $PWD/coronaviruses:/evopipe/coronaviruses`.

Logs from run will be available at `$PWD/<output dir>/info.log` (file may also be specified
using `--log` option).

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
docker run --name evopipe_file \
           -v $PWD/coronaviruses:/evopipe/coronaviruses \
           -v $PWD/species.json:/evopipe/species.json \
           -t evopipe file species.json -o coronaviruses
```

All proteomes will be stored under `fastas` directory inside docker, if you want
to perform multiple analysis without downloading proteomes each time, add mounting
`fastas` directory to `docker run` command:

```bash
docker run --name evopipe_file \
           -v $PWD/coronaviruses:/evopipe/coronaviruses \
           -v $PWD/species.json:/evopipe/species.json \
           -v $PWD/fastas:/evopipe/fastas \
           -t evopipe file species.json -o coronaviruses
```

All trees will be available at `coronaviruses/filter_trees`
an `coronaviruses/corr_trees`, i.e.: `coronaviruses/filter_trees/nj_super_tree_species.nwk`

### Family name

Providing family name, proteomes (each for one organism) will be downloaded
(sorted by a score - so from best to worst) and then used to build trees:

Run pipeline with:
```bash
# run for family analysis
docker run --name evopipe_family \
           -v $PWD/Coronaviridae:/evopipe/Coronaviridae \
           -t evopipe family Coronaviridae -n 100
```

With above command, we'll obtain maximum 100 proteomes from _Coronaviridae_.
All trees will be available at `Coronaviridae/filter_trees` and `Coronaviridae/corr_trees`,
i.e.: `Coronaviridae/filter_trees/nj_super_tree_species.nwk`

## Commands

```
usage: pipe.py [-h] [-n NUM] [--cluster-min CLUSTER_MIN] [--cluster-highest CLUSTER_HIGHEST] [--cluster-min-species-part CLUSTER_MIN_SPECIES_PART] [--filter-min FILTER_MIN] [--filter-max FILTER_MAX]
               [--fastas-dir FASTAS_DIR] [--duplications] [--super-search] [--cpu CPU] [-l LOG] [-o OUTPUT]
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
  --cluster-min-species-part CLUSTER_MIN_SPECIES_PART
                        what part of all species should be guaranteed one-to-one correspondence clusters, by default 5, so 1/5 of all species
  --filter-min FILTER_MIN
                        filter proteomes minimum
  --filter-max FILTER_MAX
                        filter proteomes maximum
  --fastas-dir FASTAS_DIR
                        directory name with fasta files, by default: "fastas/"
  --duplications        allow duplications (paralogs)
  --super-search        use more exhaustive search for super trees
  --cpu CPU             specify how many cores use for parallel computations
  -l LOG, --log LOG     logger file
  -o OUTPUT, --output OUTPUT
                        output directory, by default: name of family if "family" mode, otherwise "results"
```

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

## Development

It's highly recommended running pipeline in Docker, but if not, following packages are required with
[Conda](https://docs.conda.io/) environment:

- [RAxML](https://github.com/stamatak/standard-RAxML) (must be in `$PATH` as `raxml`)
- [ninja](http://nimbletwist.com/software/ninja/index.html) (must be in `$PATH` as `ninja`)
- [clann](https://github.com/ChrisCreevey/clann) (must be in `$PATH` as `clann`)
- [muscle](https://anaconda.org/bioconda/muscle) (`conda install -c bioconda muscle`)
- [mmseqs2](https://github.com/soedinglab/MMseqs2) (`conda install -c conda-forge -c bioconda mmseqs2`)
- [biopython](https://biopython.org/) (`conda install -c conda-forge biopython`)
- [ete3](http://etetoolkit.org/) (`conda install -c etetoolkit ete3`)
- [joblib](https://joblib.readthedocs.io/) (`conda install -c anaconda joblib`)
- [requests](https://requests.readthedocs.io/en/master/) (`conda install -c anaconda requests`)
