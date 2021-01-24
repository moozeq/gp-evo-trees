import re
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import List, Tuple

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

transitions = {
    'A': ['G'],
    'C': ['T'],
    'T': ['C'],
    'G': ['A'],
}


transversions = {
    'A': ['C', 'T'],
    'C': ['A', 'G'],
    'T': ['A', 'G'],
    'G': ['C', 'T'],
}


def is_transition(n1: str, n2: str) -> bool:
    return n1 in transitions[n2]


def is_transversion(n1: str, n2: str) -> bool:
    return n1 in transversions[n2]


def download_sequence(seq_id: str, seq_range: (int, int) = None, directory: str = '', known_orgs: dict = {}) -> str:
    """Download sequence by id"""
    directory = directory if directory.endswith('/') else f'{directory}/'
    seq_filename = f'{directory}{seq_id}.fasta'

    if Path(seq_filename).exists():
        return seq_filename

    from Bio import Entrez
    print(f'[*] Downloading sequence {seq_id}')
    Entrez.email = 'A.N.Other@example.com'
    handle = Entrez.efetch(db='nucleotide', id=seq_id, rettype='fasta', retmode='text')
    record = handle.read()

    if not record:
        raise Exception(f'[-] Downloading sequence {seq_id} FAILED')

    Path(seq_filename).parent.mkdir(exist_ok=True)
    # save full fasta
    with open(seq_filename, 'w') as f:
        f.write(record)

    # change names
    from Bio import SeqIO
    rec = next(SeqIO.parse(seq_filename, 'fasta'))
    full_id = rec.id
    # truncate seq if needed
    if seq_range:
        start, stop = seq_range
        rec = rec[start:stop]
        full_id = f'{rec.id}:{start + 1}-{stop}'
    # replace spaces and get only first part of description (before ',')
    rec.id = rec.description.replace(' ', '_').split(',')[0]
    rec.description = ''
    # check if can swap all id to organism name
    rec.id = known_orgs.get(full_id, rec.id).replace(' ', '_')

    SeqIO.write(rec, seq_filename, 'fasta')
    print(f'[+] Downloaded sequence stored at {seq_filename}')
    return seq_filename


def retrieve_genes_data(id_list):
    """Annotates Entrez Gene IDs using Bio.Entrez, in particular epost (to
    submit the data to NCBI) and esummary to retrieve the information.
    Returns a list of dictionaries with the annotations."""

    from Bio import Entrez
    Entrez.email = "A.N.Other@example.com"
    request = Entrez.epost("gene", id=",".join(id_list))
    try:
        result = Entrez.read(request)
    except RuntimeError as e:
        # FIXME: How generate NAs instead of causing an error with invalid IDs?
        print("An error occurred while retrieving the annotations.")
        print("The error returned was %s" % e)
        exit(-1)

    webEnv = result["WebEnv"]
    queryKey = result["QueryKey"]
    data = Entrez.esummary(db="gene", webenv=webEnv, query_key=queryKey)
    annotations = Entrez.read(data)

    print("Retrieved %d annotations for %d genes" % (len(annotations), len(id_list)))

    def get_gene_data(genomic_inf):
        return (
            genomic_inf['ChrAccVer'],              # gene id
            (int(genomic_inf['ChrStart']) + 1,     # gene start loc
             int(genomic_inf['ChrStop']) + 1)      # gene stop loc
        )

    genes_anns = []
    for gene_data in annotations['DocumentSummarySet']['DocumentSummary']:
        genomic_info = gene_data['GenomicInfo'][0]
        gene_info = get_gene_data(genomic_info)
        genes_anns.append(gene_info)

    return genes_anns


def gene_to_accession(query: str):
    """Based on provided search string, get gene number and translate it to NBCI accession number"""
    from Bio import Entrez
    Entrez.email = "A.N.Other@example.com"
    handle = Entrez.esearch(db="gene", term=query, idtype="acc")
    records = Entrez.read(handle)
    ids = records['IdList']
    anns = retrieve_genes_data(ids)
    print(f'"{query}" = {ids}')
    if len(ids) == 1 and len(anns) == 1:
        gene_id, (loc_start, loc_stop) = anns[0]
        gene_name = f'{gene_id}:{loc_start}-{loc_stop}'
        return gene_name
    else:
        raise Exception('[-] Gene ID must be single and be annotated')


def get_16S(species: str):
    cmd = f'"16S ribosomal RNA"[All Fields] AND "{species}"[porgn] AND ("source mitochondrion"[property] AND alive[prop])'
    return gene_to_accession(cmd)


def get_gene_homologous(gene: str, output_dir: str, limit: int = 100):
    from Bio import SeqIO
    from Bio.Blast import NCBIWWW

    path = Path(output_dir)
    path.mkdir(exist_ok=True)

    blast_results_filename = f'{output_dir}/results.xml'
    if not Path(blast_results_filename).exists():
        rec = next(SeqIO.parse(gene, 'fasta'))
        print('[*] Blastp-ing HBA1, please wait...')
        result_handle = NCBIWWW.qblast('blastp', 'nr', rec.seq, hitlist_size=1000)
        with open(blast_results_filename, 'w') as save_file:
            blast_results = result_handle.read()
            save_file.write(blast_results)
            print(f'[*] Blastp-ed, results at {blast_results_filename}')

    from Bio import SearchIO
    results = SearchIO.read(blast_results_filename, 'blast-xml')

    Path(f'{output_dir}/fastas').mkdir(exist_ok=True)
    organisms = set()
    seqs_files = []

    def good_org(org_str: str):
        return (org_str and
                org_str not in organisms and
                org_str != 'synthetic construct' and
                # 'bacter' not in org_str and  ### uncomment to exclude bacterias
                'unclassified' not in org_str)

    for hit in results.hits:
        orgs = re.findall(r"\[([A-Za-z ]+)\]", hit.description)
        org = orgs[0] if orgs else None
        if not good_org(org):
            continue
        organisms.add(org)
        rec = hit.hsps[0].hit
        rec.description = ''
        rec.id = org.replace(' ', '_')
        print(f'{org} = {rec.id}')
        seq_file = f'{output_dir}/fastas/{rec.id}.fasta'
        SeqIO.write(rec, seq_file, 'fasta')
        seqs_files.append(seq_file)

        # limit species
        if len(organisms) == limit:
            break

    with open(f'{output_dir}/species.txt', 'w') as f:
        f.write('\n'.join(organisms))

    return seqs_files


def read_records(filename: str):
    """Get records from file as Seqs objects"""
    from Bio import SeqIO
    seqs = [record for record in SeqIO.parse(filename, 'fasta')]
    return seqs


def count_records(filename: str):
    """Count how many records in fasta file"""
    from Bio import SeqIO
    return sum(True for _ in SeqIO.parse(filename, 'fasta'))


def save_fasta_records(recs, filename):
    """Save records as single fasta"""
    from Bio import SeqIO
    SeqIO.write(recs, filename, 'fasta')
    return filename


def read_sequences(filename: str) -> List[str]:
    """Get sequences from file as list of strings"""
    seqs = [str(record.seq) for record in read_records(filename)]
    return seqs


def align_records(filename: str, output: str):
    """Align records using muscle"""
    from Bio.Align.Applications import MuscleCommandline
    cline = MuscleCommandline(input=filename, out=output)
    subprocess.run(str(cline).split())


def mmseqs2(merged_fasta: str, out: str):
    if not Path((cluster_file := f'{out}/_all_seqs.fasta')).exists():
        subprocess.run(f'mmseqs easy-cluster {merged_fasta} mmseqs2 {out}'.split())
        Path('mmseqs2_all_seqs.fasta').replace(cluster_file)
        Path('mmseqs2_cluster.tsv').replace(f'{out}/_cluster.tsv')
        Path('mmseqs2_rep_seq.fasta').replace(f'{out}/_rep_seq.fasta')

    # def get_clusters():
    #     clusters = defaultdict(list)
    #     with open(f'{out}/_cluster.tsv') as f:
    #         for line in f:
    #             cat, g = line.strip().split('\t')
    #             clusters[cat].append(g)

    return cluster_file


def make_RAxML_trees(filename: str, output_dir: str, sub_model: str = '') -> (str, str):
    """Calculate trees using RAxML, returns filenames for ML and parsimony trees"""
    path = Path(output_dir)
    path.mkdir(exist_ok=True)
    cline = f'raxml -s {filename} -w {path.absolute()} -n results -m {sub_model} -p 12345'
    subprocess.run(str(cline).split())
    return f'{output_dir}/RAxML_bestTree.results', f'{output_dir}/RAxML_parsimonyTree.results'


def make_ninja_tree(filename: str, output_dir: str) -> str:
    """Calculate neighbour-joining tree using ninja, returns filename for NJ tree"""
    path = Path(output_dir)
    path.mkdir(exist_ok=True)
    cline = f'ninja --in {filename} --out {output_dir}/ninja_nj_tree.nwk'
    subprocess.run(str(cline).split())
    return f'{output_dir}/ninja_nj_tree.nwk'


def make_clann_supertree(trees_dir: str, output: str) -> str:
    """Make supertree using clann"""
    merged_trees = 'all_trees.ph'
    merge_trees_cline = f'for f in {trees_dir}/*; do sed \'/^$/d\' $f >> {merged_trees}; done'
    subprocess.run(str(merge_trees_cline).split())

    with open('clann_cmds', 'w') as cmds:
        cmds.write(f'execute; alltrees create weight=equal savetrees={output}')
    cline = f'clann -ln -c clann_cmds {merged_trees}'
    subprocess.run(str(cline).split())
    return output


def plot(data_lists: List[List[int]], data_labels: List[str], /, *,
         title: str = '', xlabel: str = '', ylabel: str = '',
         output_file: str = ''):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    for data, label in zip(data_lists, data_labels):
        ax.plot(data, label=label)
    ax.legend()
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if output_file:
        fig.savefig(output_file, dpi=fig.dpi)
    plt.show()
