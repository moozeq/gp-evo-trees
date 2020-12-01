#!/usr/bin/env python3
import argparse
import json
from shutil import copy
from pathlib import Path
from typing import Optional, Tuple

from seq_utils import (read_records, align_records, save_fasta_records, download_sequence,
                       get_16S, get_gene_homologous, make_RAxML_trees, make_ninja_tree, )


def get_HBA1_genes(args):
    import requests
    hba1 = requests.get('https://www.uniprot.org/uniprot/P69905.fasta')
    with open('P69905.fasta', 'w') as f:
        f.write(hba1.text)
    return get_gene_homologous('P69905.fasta', args.output, args.limit)


def get_16S_genes(args):
    if not args.file:
        raise Exception('No file provided')
    with open(args.file, 'r') as f:
        organisms = f.read().splitlines()

    # create directory to store results
    Path(args.output).mkdir(exist_ok=True)

    if (file := Path(f'{args.output}/genes.json')).exists():
        with file.open('r') as f:
            genes_names = json.load(f)
    else:
        genes_names = {get_16S(org): org for org in organisms}
        with file.open('w') as f:
            json.dump(genes_names, f, indent=4)

    def get_loc(name: str) -> Optional[Tuple[int, int]]:
        seq_id, locs = name.split(':')
        if not locs:
            return None
        start, end = locs.split('-')
        return int(start) - 1, int(end)  # starts with 1

    def get_id(name: str) -> str:
        seq_id, locs = name.split(':')
        return seq_id

    seqs_ids = {get_id(gene_name): get_loc(gene_name) for gene_name in genes_names}

    seqs_files = [
        download_sequence(seq_id, loc, directory=f'{args.output}/fastas', known_orgs=genes_names)
        for seq_id, loc in seqs_ids.items()
    ]
    return seqs_files


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Building evolutionary trees based on 16S rRNA or HBA1')
    parser.add_argument('mode', type=str, choices=['16S', 'HBA1'], help='basis for building evolutionary tree')
    parser.add_argument('-f', '--file', type=str, help='file with organisms names in case of 16S rRNA mode')
    parser.add_argument('-l', '--limit', type=int, default=100, help='limit organisms count for HBA1-based trees building')
    parser.add_argument('-o', '--output', type=str, default='results', help='output directory, default "results"')
    args = parser.parse_args()

    if args.mode == '16S':
        sequences_files = get_16S_genes(args)
        sub_model = 'GTRGAMMA'
    elif args.mode == 'HBA1':
        sequences_files = get_HBA1_genes(args)
        sub_model = 'PROTGAMMAGTR'
    else:
        raise Exception('Wrong mode')

    recs = []
    for seq_f in sequences_files:
        recs += read_records(seq_f)

    # save merged records
    merged_fasta = f'{args.output}/merged.fasta'
    if not Path(merged_fasta).exists():
        save_fasta_records(recs, merged_fasta)

    aligned_fasta = f'{args.output}/aligned.fasta'
    if not Path(aligned_fasta).exists():
        align_records(merged_fasta, aligned_fasta)

    tree_ml, tree_mp = make_RAxML_trees(aligned_fasta, args.output, sub_model)
    tree_nj = make_ninja_tree(aligned_fasta, args.output)

    # collect trees
    trees_dir = f'{args.output}/trees'
    Path(trees_dir).mkdir(exist_ok=True)
    copy(Path(tree_ml), Path(f'{trees_dir}/tree_ml.nwk'))
    copy(Path(tree_mp), Path(f'{trees_dir}/tree_mp.nwk'))
    copy(Path(tree_nj), Path(f'{trees_dir}/tree_nj.nwk'))

    print(f'[+] analysis results saved under directory: {args.output}, trees planted at: {args.output}/trees')
