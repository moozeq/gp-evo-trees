#!/usr/bin/env python3
import argparse
import itertools
import json
import logging
import shutil
import subprocess
from string import ascii_lowercase

import requests
from collections import defaultdict
from pathlib import Path
from typing import List, Dict, Tuple, Set

from Bio import Phylo
from Bio.Phylo.Consensus import majority_consensus
from Bio.Align.Applications import MuscleCommandline
from Bio.SeqRecord import SeqRecord
from joblib import Parallel, delayed


class RecUtils:
    """Utilities for operations on FASTA records, species names etc."""

    @staticmethod
    def read_records(filename: str):
        """Get records from file as Seqs objects."""
        from Bio import SeqIO
        seqs = [record for record in SeqIO.parse(filename, 'fasta')]
        return seqs

    @staticmethod
    def save_fasta_records(recs, filename):
        """Save records as single fasta."""
        from Bio import SeqIO
        SeqIO.write(recs, filename, 'fasta')
        return filename

    @staticmethod
    def count_records(filename: str):
        """Count how many records in fasta file."""
        from Bio import SeqIO
        return sum(True for _ in SeqIO.parse(filename, 'fasta'))

    @staticmethod
    def load_species(filename: str) -> List[str]:
        """Load species from json file."""
        with open(filename) as f:
            return json.load(f)

    @staticmethod
    def normalize_species(sp: str) -> str:
        """Change species name, removing characters which may cause issues in pipeline."""
        not_valid = [' ', '-', '/', '(', ')', '#', ':', ',', ';', '[', ']', '\'', '"', '___', '__']
        for ch in not_valid:
            sp = sp.replace(ch, '_')
        return sp

    @staticmethod
    def generate_ids():
        """Generate inner IDs using for species while building trees.

        For each species generate its inner IDs used when building trees to omit errors connected
        with species names. Generated IDs are in format as presented below:

            'A', 'B', 'C', ..., 'Z', 'AA', 'AB', ...
        """
        for size in itertools.count(1):
            for s in itertools.product(ascii_lowercase, repeat=size):
                yield ''.join(s).upper()

    @staticmethod
    def retrieve_species_names(tree_file, sp_map: Dict[str, str], out: str, rm_zero_lengths: bool = False) -> str:
        """Replace inner IDs for species with their real names."""
        def remove_zero_lengths(tree_f):
            """For super trees remove zero-length branches to avoid issues with displaying trees."""
            with open(tree_f) as f:
                tree_str = f.read()
            tree_str = tree_str.replace(':0.00000', '')
            with open(tree_f, 'w') as f:
                f.write(tree_str)

        from Bio.Phylo.Newick import Tree
        try:
            tree: Tree = Phylo.read(tree_file, 'newick')
            terms = tree.get_terminals()
            for clade in terms:
                clade.name = sp_map[clade.name]
            Phylo.write(tree, out, 'newick')
            # remove lengths which cause issues when viewing tree
            if rm_zero_lengths:
                remove_zero_lengths(out)
            return out
        except Exception as e:
            logging.info(f'Could not retrieve species names for: {tree_file}, error = {str(e)}')
            return ''


class Tools:
    """External tools, mostly using subprocess.run() to be executed."""

    @staticmethod
    def make_RAxML_trees(aligned_fasta: str, output_dir: str, sub_model: str = 'PROTGAMMAGTR') -> (str, str):
        """Calculate trees using RAxML, returns filenames for ML and parsimony trees."""
        name = Path(aligned_fasta).name[:-len('.fasta')]
        output_dir = f'{output_dir}/{name}'
        (path := Path(output_dir)).mkdir(exist_ok=True)
        ml_tree_nwk, mp_tree_nwk = f'{output_dir}/RAxML_bestTree.results', f'{output_dir}/RAxML_parsimonyTree.results'

        # already exist
        if Path(ml_tree_nwk).exists() and Path(mp_tree_nwk).exists():
            return ml_tree_nwk, mp_tree_nwk, name

        cline = f'raxml -s {aligned_fasta} -w {path.absolute()} -n results -m {sub_model} -p 12345'
        try:
            subprocess.run(str(cline).split(), check=True)
            return ml_tree_nwk, mp_tree_nwk, name
        except subprocess.CalledProcessError:
            logging.error(f'Could not build ML an MP trees for: {aligned_fasta}')
            shutil.rmtree(path.absolute())
            return ''

    @staticmethod
    def make_ninja_tree(aligned_fasta: str, output_dir: str) -> str:
        """Calculate neighbour-joining tree using ninja, returns filename for NJ tree."""
        name = Path(aligned_fasta).name[:-len('.fasta')]
        nj_tree_nwk = f'{output_dir}/{name}'

        # already exist
        if Path(nj_tree_nwk).exists():
            return nj_tree_nwk

        cline = f'ninja --in {aligned_fasta} --out {nj_tree_nwk}'
        try:
            subprocess.run(str(cline).split(), check=True)
            return nj_tree_nwk
        except subprocess.CalledProcessError:
            logging.error(f'Could not build NJ tree for: {aligned_fasta}')
            return ''

    @staticmethod
    def make_clann_super_tree(trees_dir: str, out_tree_nwk: str, super_search: bool = False) -> str:
        """Make supertree using clann, if super_search - perform more exhaustive search."""
        def get_config_str(ss: bool):
            """Get config as single string which will be saved to file."""
            hs_config_params = {
                'swap': 'nni',
                'maxswaps': 100000,
                'nreps': 5,
                'weight': 'equal',
            }
            if ss:
                hs_config_params['swap'] = 'spr'
                hs_config_params['maxswaps'] = 500000
                hs_config_params['nreps'] = 3
            hs_config_params_str = ' '.join(f'{p}={v}' for p, v in hs_config_params.items())
            config_str = f'execute; hs {hs_config_params_str} savetrees={out_tree_nwk}'
            return config_str

        if Path(out_tree_nwk).exists():
            return out_tree_nwk

        config = get_config_str(super_search)

        merged_trees = f'{trees_dir}/_alltrees.ph'
        cmds_file = f'{trees_dir}/_clanncmds'

        try:
            with open(merged_trees, 'w') as f:
                for tree in Path(trees_dir).glob('*.nwk'):
                    with open(tree) as ft:
                        tree_w_only_one_nl = ft.read().replace('\n', '')
                        tree_w_only_one_nl = f'{tree_w_only_one_nl}\n'
                        f.write(tree_w_only_one_nl)

            with open(cmds_file, 'w') as cmds:
                cmds.write(config)

            cline = f'clann -n -c {cmds_file} {merged_trees}'
            subprocess.run(str(cline).split(), check=True)
            return out_tree_nwk
        except subprocess.CalledProcessError as e:
            logging.error(f'Could not build super-tree for: {trees_dir}, err = {str(e)}')
            return ''

    @staticmethod
    def make_phylo_consensus_tree(trees_dir: str, out_tree_nwk: str) -> str:
        """Make consensus tree using biopython phylo package."""
        if Path(out_tree_nwk).exists():
            return out_tree_nwk

        trees_files = Path(trees_dir).glob(f'*.nwk')
        try:
            trees = [Phylo.read(tf, 'newick') for tf in trees_files]
            majority_tree = majority_consensus(trees)
            Phylo.write(majority_tree, out_tree_nwk, 'newick')
        except Exception as e:
            logging.error(f'Could not build consensus tree for: {trees_dir}, err = {str(e)}')
            return ''
        return majority_tree

    @staticmethod
    def mmseqs2(merged_fasta: str, out: str):
        """Cluster sequences in one, merged FASTA file using mmseqs2."""
        if not Path((cluster_file := f'{out}/_all_seqs.fasta')).exists():
            subprocess.run(f'mmseqs easy-cluster {merged_fasta} mmseqs2 {out}'.split())
            shutil.move('mmseqs2_all_seqs.fasta', cluster_file)
            shutil.move('mmseqs2_cluster.tsv', f'{out}/_cluster.tsv')
            shutil.move('mmseqs2_rep_seq.fasta', f'{out}/_rep_seq.fasta')

        return cluster_file


class Uniprot:
    """Functions calling Uniprot API to retrieve proteomes or their IDs."""

    @staticmethod
    def get_proteome_id_by_organism(organism: str) -> str:
        """Get ID of best proteome for specified organism."""
        query = f'query=organism:{organism}&format=list&sort=score'
        url = f'https://www.uniprot.org/proteomes/?{query}'
        try:
            ids = requests.get(url)
            ids = ids.content.decode()
            if not ids:
                raise Exception('empty list')
            pid = ids.splitlines()[0]
            if not pid.startswith('UP'):
                raise Exception(f'wrong pid = {pid}')
            logging.info(f'Get proteome ID: {organism} -> {pid}')
            return pid
        except Exception as e:
            logging.error(f'Could not download proteome IDs list for: {organism}, error = {str(e)}')
            return ''

    @staticmethod
    def download_proteomes_ids(tax: str, out: str) -> str:
        """Download all proteomes IDs for specific tax family."""
        query = f'query=taxonomy:{tax}&format=tab&sort=score'
        url = f'https://www.uniprot.org/proteomes/?{query}'
        try:
            if Path(out).exists():
                return out
            ids = requests.get(url)
            ids = ids.content.decode()
            if not ids:
                raise Exception('empty list')
            logging.info(f'Downloaded proteomes IDs list: {len(ids.splitlines()) - 1}')
            with open(out, 'w') as fp:
                fp.write(ids)
            return out
        except Exception as e:
            logging.error(f'Could not download proteomes IDs list for: {tax}, error = {str(e)}')
            return ''

    @staticmethod
    def download_proteome(pid: str, org: str, o_dir: str):
        """Download proteome using proteome Uniprot ID."""
        query = f'query=proteome:{pid}&format=fasta&compress=no'
        url = f'https://www.uniprot.org/uniprot/?{query}'
        try:
            if Path(pfile := f'{o_dir}/{RecUtils.normalize_species(org)}.fasta').exists():
                return pfile
            ids = requests.get(url)
            ids = ids.content.decode()
            if not ids:
                raise Exception('empty proteome')
            logging.info(f'Downloaded proteome for: {org}')
            with open(pfile, 'w') as fp:
                fp.write(ids)
            return pfile
        except Exception as e:
            logging.error(f'Could not download proteome for: {org}, error = {str(e)}')
            return ''


def download_proteomes_by_names(names: List[str], fastas_out: str, limit: int = 100000) -> List[str]:
    """Providing list of organisms names, try to download proteomes with max limit."""
    pids = {
        pid: org
        for org in names
        if (
            not Path(f'{fastas_out}/{RecUtils.normalize_species(org)}.fasta').exists() and
            (pid := Uniprot.get_proteome_id_by_organism(org))
        )
    }

    proteomes_files = {
        org: str(prot_file)
        for org in names
        if (prot_file := Path(f'{fastas_out}/{RecUtils.normalize_species(org)}.fasta')).exists()
    }

    if not pids and not proteomes_files:
        raise Exception('No proteome IDs loaded')

    if len(proteomes_files) >= limit:
        logging.info(f'Used only local proteomes (limit = {limit}): {len(proteomes_files)}/{len(names)}')
        return list(proteomes_files.values())[:limit]

    logging.info(f'Translated organisms names to proteomes IDs: {len(pids) + len(proteomes_files)}/{len(names)}')
    for i, (pid, org) in enumerate(pids.items()):
        if (
            len(proteomes_files) < limit and
            org not in proteomes_files and
            (prot_file := Uniprot.download_proteome(pid, org, fastas_out))
        ):
            proteomes_files[org] = prot_file

    logging.info(f'Downloaded proteomes for: {len(proteomes_files)}/{len(names)}')
    return list(proteomes_files.values())


def download_proteomes_by_family(family: str, fastas_out: str, limit: int = 100000) -> List[str]:
    """Providing taxonomy family name, try to download proteomes from it with max limit."""
    ids_file = Uniprot.download_proteomes_ids(family, f'{fastas_out}/_ids.tsv')
    proteomes_files = {}
    with open(ids_file) as ifp:
        reader = iter(ifp)
        header = next(reader)
        for i, entry in enumerate(reader):
            pid, org, *_ = entry.split('\t')
            if (
                len(proteomes_files) < limit and
                org not in proteomes_files and
                (prot_file := Uniprot.download_proteome(pid, org, fastas_out))
            ):
                proteomes_files[org] = prot_file

    logging.info(f'Downloaded proteomes for: {len(proteomes_files)}/{i}')
    Path(ids_file).unlink(missing_ok=True)  # remove downloaded file with IDs after work
    return list(proteomes_files.values())


def download_proteomes(mode: str, input: str, fastas_out: str, limit: int = 100000) -> List[str]:
    """Base on mode, set proper proteomes downloading option: by family or by file."""
    Path(fastas_out).mkdir(exist_ok=True)
    if mode == 'family':
        return download_proteomes_by_family(input, fastas_out, limit)
    elif mode == 'file':
        with open(input) as f:
            organisms = json.load(f)
        return download_proteomes_by_names(organisms, fastas_out, limit)
    else:
        raise Exception(f'Wrong mode: mode = {mode}, input = {input}')


def filter_fastas(fastas: List[str], min_seqs: int = 0, max_seqs: int = 100000) -> List[str]:
    """Simple fasta filtering, all fastas without meeting requirements will be excluded."""
    if min_seqs == 0 and max_seqs == -1:
        return fastas
    filtered_fastas = [
        file
        for file in fastas
        if max_seqs > RecUtils.count_records(file) >= min_seqs
    ]
    logging.info(f'Filtered fastas with min = {min_seqs}, max = {max_seqs}: {len(filtered_fastas)}/{len(fastas)}')
    return filtered_fastas


def map_recs_to_species(fastas: List[str], out: str) -> (Dict[str, str], Dict[str, str]):
    """Create mapping for record IDs to species names and their inner IDs created with `RecUtils.generate_ids()`."""
    if Path(out).exists():
        with open(out) as f:
            maps = json.load(f)
            return maps['recs'], maps['orgs']

    # map unique IDs 'A', 'B', ... 'aa', 'ab' to organisms names
    org_ids = RecUtils.generate_ids()
    orgs_map = {
        Path(file).name[:-len('.fasta')]: next(org_ids)
        for file in fastas
    }
    rev_orgs_map = {v: k for k, v in orgs_map.items()}

    recs_map = {}
    for file in fastas:
        seqs = RecUtils.read_records(file)
        org_name = Path(file).name[:-len('.fasta')]
        seqs = {
            seq.id: orgs_map[org_name]
            for seq in seqs
        }
        recs_map.update(seqs)
        logging.info(f'Mapped records ({len(seqs)}) IDs for: {org_name}')

    logging.info(f'Mapped records ({len(recs_map)}) to species ({len(fastas)})')
    with open(out, 'w') as f:
        json.dump({'recs': recs_map, 'orgs': rev_orgs_map}, f, indent=4)
    return recs_map, rev_orgs_map


def merge_fastas(fastas: List[str], out: str) -> str:
    """Simply merging all fasta files into one big FASTA."""
    if Path(out).exists():
        return out
    recs = [
        RecUtils.read_records(file)
        for file in fastas
    ]
    merged = list(itertools.chain.from_iterable(recs))
    RecUtils.save_fasta_records(merged, out)
    logging.info(f'Saved all records ({len(merged)}) to file: {out}')
    return out


def clustering(merged_fasta: str,
               out: str,
               recs_map: Dict[str, str],
               min_len: int = 2,
               min_species_part: int = 5,
               highest: int = 0,
               duplications: bool = False) -> (Dict[str, list], Dict[str, list]):
    """Cluster merged fasta file (proteomes) to proteins families."""
    def rename_fasta_record(seq: SeqRecord, name: str):
        """Set inner ID as sequence ID."""
        seq.description = ''
        seq.id = name
        return seq

    def get_clusters(f: str) -> Dict[str, List[SeqRecord]]:
        """Based on output from mmseqs2, get clusters as dict."""
        recs = iter(RecUtils.read_records(f))
        # order is important!
        unfiltered_clusters = defaultdict(list)
        cluster_name = next(recs).id
        for rec in recs:
            if len(rec.seq) == 0:  # cluster sequence (id only)
                cluster_name = rec.id
            else:  # real sequence after cluster sequence
                rename_fasta_record(rec, recs_map[rec.id])
                unfiltered_clusters[cluster_name].append(rec)

        logging.info(f'Loaded clusters: {len(unfiltered_clusters)}')

        return unfiltered_clusters

    def filter_duplications(cls_recs, dup):
        """Filter correspondence, if dup == True, it means allow duplications (paralogs) in clusters."""
        if dup:
            return cls_recs
        # remove duplicates, one-to-one correspondence
        ids = set()
        corr_cls_recs = []
        for cls_rec in cls_recs:
            if cls_rec.id in ids:
                continue
            corr_cls_recs.append(cls_rec)
            ids.add(cls_rec.id)
        return corr_cls_recs

    def filter_clusters(cls: Dict[str, List[SeqRecord]], lim: int, hg: int, dup: bool) -> Dict[str, list]:
        """Filter out clusters with min length (lim), first highest (hg), duplications allowance (dup)."""
        def filter_highest(cls_recs, hg_filter: int):
            """Get only first n most populated clusters."""
            if hg_filter:
                clusters_from_highest = sorted(cls_recs, key=lambda k: len(cls_recs[k]), reverse=True)
                hg_lim = hg if hg < len(clusters_from_highest) else len(clusters_from_highest)
                cls_recs = {k: cls_recs[k] for k in clusters_from_highest[:hg_lim]}
            return cls_recs

        unfiltered_records_cnt = sum(len(fc) for fc in cls.values())

        # filter out clusters without meeting minimum len
        filtered_clusters = {
            cluster_name: f_cluster_recs
            for cluster_name, cluster_recs in cls.items()
            if len(f_cluster_recs := filter_duplications(cluster_recs, duplications)) >= lim
        }
        # get only first n highest clusters (most populated) if specified
        filtered_clusters = filter_highest(filtered_clusters, hg)

        filtered_records_cnt = sum(len(fc) for fc in filtered_clusters.values())
        logging.info(
            f'Filtered clusters with duplication = {dup}, min_len = {lim}, highest = {highest}: '
            f'clusters {len(filtered_clusters)}/{len(cls)}, records {filtered_records_cnt}/{unfiltered_records_cnt}')
        return filtered_clusters

    def filter_one_to_one(cls: Dict[str, List[SeqRecord]], lim: int, min_sp_part: int) -> Dict[str, list]:
        """Guarantee one-to-one correspondence for minimum species specified"""
        def prune_cluster(cls_recs_to_prune: List[SeqRecord], only_species: Set[str]):
            pruned_cls_recs = [
                rec
                for rec in cls_recs_to_prune
                if rec.id in only_species
            ]
            return pruned_cls_recs

        # at least part of all species must remain after filtering
        species_min = (all_species_cnt := len(set(recs_map.values()))) // min_sp_part

        unfiltered_records_cnt = sum(len(fc) for fc in cls.values())

        # filter out clusters without meeting minimum len
        filtered_clusters = {
            cluster_name: f_cluster_recs
            for cluster_name, cluster_recs in cls.items()
            if len(f_cluster_recs := filter_duplications(cluster_recs, False)) >= lim
        }

        # guarantee correspondence one-to-one
        clusters_from_highest = sorted(filtered_clusters, key=lambda k: len(filtered_clusters[k]), reverse=True)
        species = set()
        one_to_one_clusters = {}
        # start from most populated clusters and add clusters till species_min included
        for cls_name in clusters_from_highest:
            cur_cls: List[SeqRecord] = filtered_clusters[cls_name]
            cls_species = {sp.id for sp in cur_cls}
            # first iteration, biggest family comes first
            if not species:
                species = cls_species
                one_to_one_clusters[cls_name] = cur_cls
                continue
            # still can add clusters since min species correspondence guaranteed
            if len(cls_species_intersect := species.intersection(cls_species)) >= species_min:
                species = cls_species_intersect
                one_to_one_clusters[cls_name] = cur_cls
            else:
                break

        # remove species which occurs in filtered clusters but without one-to-one correspondence
        filtered_one_to_one_clusters = {
            cls_name: pruned_cls_recs
            for cls_name, cls_recs in one_to_one_clusters.items()
            if (pruned_cls_recs := prune_cluster(cls_recs, species))
        }

        corr_records_cnt = sum(len(fc) for fc in filtered_one_to_one_clusters.values())
        logging.info(
            f'Filtered clusters one-to-one correspondence without duplication, min_len = {lim}, '
            f'min species = {species_min}: '
            f'clusters {len(filtered_one_to_one_clusters)}/{len(cls)}, '
            f'species {len(species)}/{all_species_cnt}, '
            f'records {corr_records_cnt}/{unfiltered_records_cnt}'
        )
        return filtered_one_to_one_clusters

    clusters_file = Tools.mmseqs2(merged_fasta, out)
    clusters = get_clusters(clusters_file)
    filter_clusters = filter_clusters(clusters, min_len, highest, duplications)
    corr_clusters = filter_one_to_one(clusters, min_len, min_species_part)

    logging.info(
        f'Clustered records, filtered clusters: {len(filter_clusters)}/{len(clusters)}, '
        f'correspondence guaranteed clusters: {len(corr_clusters)}/{len(clusters)}'
    )
    return filter_clusters, corr_clusters


def make_genes_families(clusters: Dict[str, list], out: str) -> List[str]:
    """Based od clusters dict, save fasta files with protein families."""
    Path(out).mkdir(exist_ok=True)

    families = []
    for cls_name, cls_recs in clusters.items():
        if not Path(family_filename := f'{out}/{cls_name}.fasta').exists():
            RecUtils.save_fasta_records(cls_recs, family_filename)
        families.append(family_filename)

    logging.info(f'Created protein families from clusters: {len(families)}/{len(clusters)}')
    return families


def align_families(families: List[str], out: str) -> List[str]:
    """Align protein families fastas within them, not with each other."""
    Path(out).mkdir(exist_ok=True)

    def get_output_filename(fasta: str, output: str) -> str:
        return f'{output}/{Path(fasta).name}'

    def align_fasta_file(fasta_in: str, fasta_out: str):
        if Path(fasta_out).exists():
            return fasta_out
        cline = MuscleCommandline(input=fasta_in, out=fasta_out)
        try:
            subprocess.run(str(cline).split(), check=True)
            return fasta_out
        except subprocess.CalledProcessError:
            logging.error(f'Could not align fasta file: {fasta_in}')
            return ''

    aligned = Parallel(n_jobs=args.cpu)(delayed(align_fasta_file)(
        fasta, get_output_filename(fasta, out)
    ) for fasta in families)

    aligned = [fasta_file for fasta_file in aligned if fasta_file]
    logging.info(f'Aligned families: {len(aligned)}/{len(families)}')
    return aligned


def build_trees(aligned_fastas: List[str], out: str, super_search: bool = False) -> (List[str], List[str]):
    """Core of pipeline, building NJ, ML and MP trees.

    Trees are built using:
        - [NJ] neighbor-joining:    using ninja
        - [ML] maximum-likelihood:  using RAxML
        - [MP] maximum-parsimony:   using RAxML

    Args:
        aligned_fastas: filenames for aligned fasta files with protein families
        out:            output directory for trees
        super_search:   if True - perform exhaustive search for super trees

    Returns:
        (list, list):   two lists, corresponding to filenames for:
                            - built consensus trees (NJ, ML, MP)
                            - built super trees (NJ, ML, MP)

    """
    Path(out).mkdir(exist_ok=True)
    Path(nj_trees_dir := f'{out}/nj-trees').mkdir(exist_ok=True)
    Path(ml_trees_dir := f'{out}/ml-trees').mkdir(exist_ok=True)
    Path(mp_trees_dir := f'{out}/mp-trees').mkdir(exist_ok=True)

    nj_cons = f'{out}/nj_consensus_tree.nwk'
    ml_cons = f'{out}/ml_consensus_tree.nwk'
    mp_cons = f'{out}/mp_consensus_tree.nwk'
    nj_super = f'{out}/nj_super_tree.nwk'
    ml_super = f'{out}/ml_super_tree.nwk'
    mp_super = f'{out}/mp_super_tree.nwk'

    def unroot_tree(tree_file: str) -> str:
        """Need to unroot trees from ninja, since they are rooted, while from RAxML are unrooted."""
        try:
            cline = f'ete3 mod --unroot -t {tree_file}'
            proc_out = subprocess.run(cline.split(), check=True, capture_output=True)
            proc_out = proc_out.stdout.decode()
            if not proc_out:
                raise Exception(f'invalid oputput')
            # overwrite tree with unrooted version
            with open(tree_file, 'w') as f:
                f.write(proc_out)
            return tree_file
        except subprocess.CalledProcessError as e:
            logging.error(f'Unrooting failed for tree: {tree_file}, error = {str(e)}')
            return ''

    def move_ninja_trees(nj_trees: List[str]):
        """Move output trees from ninja to proper directory (NJ trees)."""
        Path(nj_trees_dir).mkdir(exist_ok=True)
        for nj_tree in nj_trees:
            family = Path(nj_tree).name
            shutil.move(nj_tree, f'{nj_trees_dir}/{family}.nwk')

    def move_raxml_trees(raxml_trees: List[Tuple[str, str, str]]):
        """Move output trees from RAxML to proper directories (MP and ML trees)."""
        Path(ml_trees_dir).mkdir(exist_ok=True)
        Path(mp_trees_dir).mkdir(exist_ok=True)
        for ml_tree, mp_tree, family in raxml_trees:
            shutil.move(ml_tree, f'{ml_trees_dir}/{family}.nwk')
            shutil.move(mp_tree, f'{mp_trees_dir}/{family}.nwk')

    def make_ninja_trees():
        """Make NJ trees using ninja.

        Make neighbor-joining trees using ninja, then unroot them and move to proper directory.
        Unroot trees for consistency with RAxML ML and MP trees.

        """
        Path(f'{out}/ninja').mkdir(exist_ok=True, parents=True)
        ninja_trees = Parallel(n_jobs=args.cpu)(
            delayed(Tools.make_ninja_tree)
            (a_fasta, f'{out}/ninja')
            for a_fasta in aligned_fastas
        )
        ninja_trees = [tree for tree in ninja_trees if tree]
        unrooted_ninja_trees = [
            unrooted_tree
            for tree in ninja_trees
            if (unrooted_tree := unroot_tree(tree))
        ]
        logging.info(f'Unrooted ninja NJ trees: {len(unrooted_ninja_trees)}/{len(ninja_trees)}')
        move_ninja_trees(unrooted_ninja_trees)
        logging.info(f'Built unrooted NJ trees using ninja '
                     f'(store at {nj_trees_dir}): {len(ninja_trees)}')

    def make_raxml_trees():
        """Make ML and MP trees using RAxML.

        Make unrooted maximum-likelihood and maximum-parsimony trees using RAxML,
        then move them to proper directory (no need to unroot like from ninja).

        """
        Path(f'{out}/raxml').mkdir(exist_ok=True, parents=True)
        raxml_trees = Parallel(n_jobs=args.cpu)(
            delayed(Tools.make_RAxML_trees)
            (a_fasta, f'{out}/raxml')
            for a_fasta in aligned_fastas
        )
        raxml_trees = [tree for tree in raxml_trees if tree]
        move_raxml_trees(raxml_trees)
        logging.info(f'Built unrooted ML and MP trees using RAxML '
                     f'(store at {ml_trees_dir}, {mp_trees_dir}): {len(raxml_trees)}')

    def make_consensus_trees(o_nj: str, o_ml: str, o_mp: str):
        """Make NJ, ML and MP consensus trees."""
        parameters = [
            (nj_trees_dir, o_nj),
            (ml_trees_dir, o_ml),
            (mp_trees_dir, o_mp),
        ]
        Parallel(n_jobs=args.cpu)(
            delayed(Tools.make_phylo_consensus_tree)(  # make_phylo_consensus_tree
                t_dir, t_out
            ) for t_dir, t_out in parameters
        )

    def make_super_trees(o_nj: str, o_ml: str, o_mp: str, ss: bool):
        """Make NJ, ML and MP super trees."""
        parameters = [
            (nj_trees_dir, o_nj),
            (ml_trees_dir, o_ml),
            (mp_trees_dir, o_mp),
        ]
        Parallel(n_jobs=args.cpu)(
            delayed(Tools.make_clann_super_tree)(  # make_clann_super_tree
                t_dir, t_out, ss
            ) for t_dir, t_out in parameters
        )

    # if trees already under directory, don't make them again
    if not list(Path(nj_trees_dir).glob('*.nwk')):
        make_ninja_trees()
    if not list(Path(ml_trees_dir).glob('*.nwk')) and not list(Path(mp_trees_dir).glob('*.nwk')):
        make_raxml_trees()

    make_consensus_trees(nj_cons, ml_cons, mp_cons)
    make_super_trees(nj_super, ml_super, mp_super, super_search)

    return [nj_cons, ml_cons, mp_cons], [nj_super, ml_super, mp_super]


def retrieve_species_names(
        trees_files: List[str],
        orgs_map: Dict[str, str],
        rm_zero_lengths: bool = False) -> List[str]:
    """Retrieve real species names instead of inner IDs.

    Since consensus and super trees are built using inner IDs for species ('A', 'B', 'AA', ...)
    we need to retrieve real species names using `orgs_map` created before.

    Also we remove zero-length branches from supertree, since Phylo while writes trees, writes
    also branches lengths sometimes causing issues with displaying them.

    Returns:
        list: trees filenames with successfully retrieved species names

    """
    def get_tree_only(tree: str) -> str:
        """Supertrees are written with score after ';', remove it to be able to read newick tree."""
        tree_only = tree.split(';')[0]
        tree_only = f'{tree_only};'
        return tree_only

    def prune_trees(trees: List[str]):
        """Remove scores from all trees, making them readable as newick."""
        for tree_f in trees:
            try:
                with open(tree_f) as f:
                    pruned_tree = get_tree_only(f.read())
                with open(tree_f, 'w') as f:
                    f.write(pruned_tree)
            except Exception as e:
                logging.error(f'Could not prune tree from score: {tree_f}, error = {str(e)}')

    prune_trees(trees_files)
    successfully_retrieved = []
    for tree_file in trees_files:
        tree_file_ret = f'{tree_file[:-len(".nwk")]}_species.nwk'
        ret = RecUtils.retrieve_species_names(tree_file, orgs_map, tree_file_ret, rm_zero_lengths)
        if ret:
            successfully_retrieved.append(ret)
            Path(tree_file).unlink()  # remove tree without species translation
    logging.info(f'Retrieved species names for: {len(successfully_retrieved)}/{len(trees_files)}')
    return successfully_retrieved


def set_logger(log_file: str):
    """Log info to stdout and file specified in argparse."""
    logging.root.handlers = []
    # noinspection PyArgumentList
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )


def pipeline(input_args) -> List[str]:
    """Core part of script, performing pipeline steps with provided arguments.

    Arguments are read from user using argparse, run script with `-h` to see available options.
    Following steps are performed in pipeline:
        1) setting logger to log info to stdout to file (by default under output directory)
        2) download proteomes for provided tax family, or for species names from .json file
        3) filter fasta files with min and max sequences within, specified in args
        4) change sequences IDs to inner species IDs for building trees purposes
        5) merge all proteomes into one, big fasta file
        6) cluster merged fasta file to obtain protein families from all species
        7) save clusters to separate files, each per protein family
        8) align all sequences within each protein family fasta file
        9) build NJ, ML and MP trees from provided aligned protein families fasta files
        10) retrieve species names for consensus and super trees from their inner IDs

    Returns:
        list: successfully built consensus and super NJ, ML and MP trees

    """
    Path(input_args.output).mkdir(exist_ok=True)

    set_logger(input_args.log)

    prots = download_proteomes(
        input_args.mode,
        input_args.input,
        fastas_out=input_args.fastas_dir,
        limit=input_args.num
    )
    prots = filter_fastas(
        prots,
        min_seqs=input_args.filter_min,
        max_seqs=input_args.filter_max
    )
    recs_map, orgs_map = map_recs_to_species(
        prots,
        f'{input_args.output}/_recsmap.json'
    )
    all_prots = merge_fastas(
        prots,
        f'{input_args.output}/_merged.fasta'
    )
    filtered_clusters, corr_clusters = clustering(
        all_prots,
        f'{input_args.output}/mmseqs2',
        recs_map=recs_map,
        min_len=input_args.cluster_min,
        min_species_part=input_args.cluster_min_species_part,
        highest=input_args.cluster_highest,
        duplications=input_args.duplications
    )

    def second_stage(clusters: Dict[str, List[SeqRecord]], prefix: str):
        """Second stage for pipeline to perform after clustering"""
        logging.info(f'[2nd STAGE] Started pipeline 2nd stage for clusters count: {len(clusters)}, prefix: {prefix}')
        families = make_genes_families(
            clusters,
            f'{input_args.output}/{prefix}clusters'
        )
        aligned_families = align_families(
            families,
            f'{input_args.output}/{prefix}align'
        )
        consensus_trees, super_trees = build_trees(
            aligned_families,
            f'{input_args.output}/{prefix}trees',
            super_search=input_args.super_search
        )
        consensus_trees = retrieve_species_names(
            consensus_trees,
            orgs_map=orgs_map
        )
        super_trees = retrieve_species_names(
            super_trees,
            orgs_map=orgs_map,
            rm_zero_lengths=True
        )
        logging.info(f'[2nd STAGE] Ended pipeline 2nd stage for clusters count: {len(clusters)}, prefix: {prefix}, '
                     f'built consensus trees: {consensus_trees}, '
                     f'built super trees: {super_trees}')
        return consensus_trees, super_trees

    corr_consensus_trees, corr_super_trees = second_stage(corr_clusters, 'corr_')
    filtered_consensus_trees, filtered_super_trees = second_stage(filtered_clusters, 'filter_')

    return [
        *corr_super_trees,
        *corr_consensus_trees,
        *filtered_super_trees,
        *filtered_consensus_trees,
    ]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Phylogenetic pipeline to infer a species/genome tree from a set of genomes')
    parser.add_argument('mode', type=str, choices=['family', 'file'],
                        help='pipeline mode')
    parser.add_argument('input', type=str,
                        help='family name or .json file with species names list which will be inferred')
    parser.add_argument('-n', '--num', type=int, default=100000,
                        help='limit downloading species to specific number')
    parser.add_argument('--cluster-min', type=int, default=4,
                        help='filter cluster proteomes minimum, by default: 4')
    parser.add_argument('--cluster-highest', type=int,
                        help='get only "n" most populated clusters')
    parser.add_argument('--cluster-min-species-part', type=int, default=5,
                        help='what part of all species should be guaranteed one-to-one correspondence '
                             'clusters, by default 5, so 1/5 of all species')
    parser.add_argument('--filter-min', type=int, default=0,
                        help='filter proteomes minimum')
    parser.add_argument('--filter-max', type=int, default=100000,
                        help='filter proteomes maximum')
    parser.add_argument('--fastas-dir', type=str, default='fastas',
                        help='directory name with fasta files, by default: "fastas/"')
    parser.add_argument('--duplications', action='store_true', default=False,
                        help='allow duplications (paralogs)')
    parser.add_argument('--super-search', action='store_true', default=False,
                        help='use more exhaustive search for super trees')
    parser.add_argument('--cpu', type=int, default=4,
                        help='specify how many cores use for parallel computations')
    parser.add_argument('-l', '--log', type=str, default='info.log',
                        help='logger file')
    parser.add_argument('-o', '--output', type=str,
                        help='output directory, by default: name of family if "family" mode, otherwise "results"')
    args = parser.parse_args()

    if not args.output:
        if args.mode == 'family':
            args.output = args.input
        else:
            args.output = 'results'

    args.log = f'{args.output}/{args.log}'

    built_trees = pipeline(args)
    logging.info(f'Successfully built trees: {built_trees}')
