#!/usr/bin/env python3
import random
import time
import subprocess
from pathlib import Path

import prody


def download_family(name: str) -> str:
    return prody.fetchPfamMSA(name, format='fasta', outname=f'fs/{name}.fasta', timeout=10)


def read_records(filename: str):
    """Get records from file as Seqs objects"""
    from Bio import SeqIO
    seqs = [record for record in SeqIO.parse(filename, 'fasta')]
    return seqs


def save_fasta_records(recs, filename):
    """Save records as single fasta"""
    from Bio import SeqIO
    SeqIO.write(recs, filename, 'fasta')


def get_pfam_100_families():
    if Path('../lab6/merged.fasta').exists():
        return get_clusters('../lab6/merged.fasta', 30)
    seqs = []
    fams = set()
    Path('../lab6/fs').mkdir(exist_ok=True)
    while len(seqs) != 100 * 30:
        i = random.randrange(101, 999, 1)
        try:
            fam = f'PF00{i}'
            if fam in fams:
                continue
            o = f'PF00{i}.fasta_full.fasta'
            if not Path(o).exists():
                o = download_family(f'PF00{i}')
            fams.add(fam)
            recs = read_records(o)
            if len(recs) >= 30:
                seqs += recs[:30]
            time.sleep(3)
            print(f'>>>> >>>> {len(seqs)} / {100 * 30} <<<< <<<<')
        except:
            continue
    print(len(seqs))
    save_fasta_records(seqs, 'temp_merged.fasta')
    subprocess.run('cat temp_merged.fasta | tr -d \'-\' > merged.fasta'.split())
    return get_clusters('../lab6/merged.fasta', 30)


def mmseqs2(file):
    if not Path('../lab6/mmseqs2_cluster.tsv').exists():
        subprocess.run('mmseqs easy-cluster merged.fasta mmseqs2 tmp'.split())  # mmseqs2_cluster.tsv

    clusters = {}
    with open(file) as f:
        for line in f:
            cat, g = line.strip().split('\t')
            if cat not in clusters:
                clusters[cat] = []
            clusters[cat].append(g)

    return clusters


def cdhit(file):
    if not Path('../lab6/cdhit.clstr').exists():
        subprocess.run('cd-hit -i merged.fasta -o cdhit -c 0.66'.split())  # cdhit.clstr

    def find_between(s, start, end):
        return (s.split(start))[1].split(end)[0]

    def get_gene(l):
        return find_between(l, '>', '...')

    clusters = {}
    with open(file) as f:
        for line in f:
            if 'Cluster' in line:
                clusters[(c := line.strip())] = []
                continue
            g = get_gene(line)
            clusters[c].append(g)

    return clusters


def get_clusters(file, size):
    def chunks(lst, n):
        """Yield successive n-sized chunks from lst."""
        for i in range(0, len(lst), n):
            yield lst[i:i + n]

    def get_names(seqs):
        return [seq.id for seq in seqs]

    recs = read_records(file)
    clusters = {
        i: get_names(chunk)
        for i, chunk in enumerate(chunks(recs, size))
    }
    return clusters


def jaccard_similarity(list1, list2):
    s1 = set(list1)
    s2 = set(list2)
    return float(len(s1.intersection(s2)) / len(s1.union(s2)))


def calc_j_sim(rclus: dict, clus: dict):
    done = set()
    js_all = 0.0
    for c_name, c_list in clus.items():
        for rc_name, rc_list in rclus.items():
            if (js := jaccard_similarity(rc_list, c_list)) > 0.0:
                done.add(c_name)
                js_all += js
                break
    return js_all / len(clus)


if __name__ == '__main__':
    real_clusters = get_pfam_100_families()
    cdhit_clusters = cdhit('../lab6/cdhit.clstr')
    mmseqs2_clustesrs = mmseqs2('../lab6/mmseqs2_cluster.tsv')

    cdhit_js = calc_j_sim(real_clusters, cdhit_clusters)
    mmseqs2_js = calc_j_sim(real_clusters, mmseqs2_clustesrs)
    print(f'Clustering results:\n'
          f'\t[Jaccard similarity] mmseqs2 = {mmseqs2_js:6f}\n'
          f'\t[Jaccard similarity] cd-hit = {cdhit_js:6f}')
