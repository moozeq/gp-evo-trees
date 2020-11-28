#!/usr/bin/env python3
import argparse
import json
from pathlib import Path
from typing import Optional, Tuple

from seq_utils import (read_records, align_records, save_fasta_records, download_sequence, make_trees,
                       get_16SrRNA_gene_name, )

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Building evolutionary trees based on 16S rRNA')
    parser.add_argument('file', type=str, help='file with organisms names')
    parser.add_argument('-o', '--output', type=str, default='results', help='output directory, default "results"')
    args = parser.parse_args()

    with open(args.file, 'r') as f:
        organisms = f.read().splitlines()

    # create directory to store results
    Path(args.output).mkdir(exist_ok=True)

    if (file := Path(f'{args.output}/genes.json')).exists():
        with file.open('r') as f:
            genes_names = json.load(f)
    else:
        genes_names = {get_16SrRNA_gene_name(org): org for org in organisms}
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

    seqs_files = [download_sequence(seq_id, loc, directory=f'{args.output}/fastas') for seq_id, loc in seqs_ids.items()]
    recs = []
    for seq_f in seqs_files:
        recs += read_records(seq_f)

    # save merged records
    merged_fasta = f'{args.output}/merged.fasta'
    save_fasta_records(recs, merged_fasta)

    aligned_fasta = f'{args.output}/aligned.fasta'
    align_records(merged_fasta, aligned_fasta)

    make_trees(aligned_fasta, args.output)
    print(f'[+] analysis results saved under directory: {args.output}')
