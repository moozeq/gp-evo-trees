import subprocess
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


def download_sequence(seq_id: str, seq_range: (int, int) = None, directory: str = '', known_ids: dict = {}) -> str:
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
        full_id = f'{rec.id}:{start}-{stop}'
    # replace spaces and get only first part of description (before ',')
    rec.id = rec.description.replace(' ', '_').split(',')[0]
    rec.description = ''
    # check if can swap all id to organism name
    rec.id = known_ids.get(full_id, rec.id)

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


def get_16SrRNA_gene_name(species: str):
    cmd = f'"16S ribosomal RNA"[All Fields] AND "{species}"[porgn] AND ("source mitochondrion"[property] AND alive[prop])'
    from Bio import Entrez
    Entrez.email = "A.N.Other@example.com"
    handle = Entrez.esearch(db="gene", term=cmd, idtype="acc")
    records = Entrez.read(handle)
    ids = records['IdList']
    anns = retrieve_genes_data(ids)
    print(f'{species} = {ids}')
    if len(ids) == 1 and len(anns) == 1:
        gene_id, (loc_start, loc_stop) = anns[0]
        gene_name = f'{gene_id}:{loc_start}-{loc_stop}'
        return gene_name
    else:
        raise Exception('[-] Gene ID must be single and be annotated')


def read_records(filename: str):
    """Get records from file as Seqs objects"""
    from Bio import SeqIO
    print(f'[*] Reading sequences from {filename}')
    seqs = [record for record in SeqIO.parse(filename, 'fasta')]
    print(f'[+] Read {len(seqs)} sequences')
    return seqs


def save_fasta_records(recs, filename):
    """Save records as single fasta"""
    from Bio import SeqIO
    SeqIO.write(recs, filename, 'fasta')


def read_sequences(filename: str) -> List[str]:
    """Get sequences from file as list of strings"""
    seqs = [str(record.seq) for record in read_records(filename)]
    return seqs


def align_records(filename: str, output: str):
    """Align records using muscle"""
    from Bio.Align.Applications import MuscleCommandline
    cline = MuscleCommandline(input=filename, out=output)
    subprocess.run(str(cline).split())


def make_trees(filename: str, output_dir: str):
    """Calculate trees using RAxML, returns filenames for ML and parsimony trees"""
    path = Path(output_dir)
    path.parent.mkdir(exist_ok=True)
    cline = f'raxml -s {filename} -w {path.absolute()} -n results -m GTRCATIX -p 12345'
    subprocess.run(str(cline).split())
    return f'{output_dir}/RAxML_bestTree', f'{output_dir}/RAxML_parsimonyTree'


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
