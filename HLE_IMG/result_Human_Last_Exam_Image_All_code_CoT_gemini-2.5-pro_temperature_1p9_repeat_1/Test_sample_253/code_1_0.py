# First, ensure you have biopython installed:
# pip install biopython

from Bio import Entrez
from Bio import SeqIO
import io

def get_protein_length(organism, protein_name):
    """
    Searches the NCBI protein database for a specific protein in a given organism
    and returns its length.
    """
    Entrez.email = "your.email@example.com"  # Please replace with your email address

    # Construct the search query
    search_term = f'"{organism}"[Organism] AND "{protein_name}"[Protein Name]'

    # Search the protein database
    handle = Entrez.esearch(db="protein", term=search_term, retmax="1")
    record = Entrez.read(handle)
    handle.close()

    if not record["IdList"]:
        return f"Protein '{protein_name}' not found for organism '{organism}'."

    protein_id = record["IdList"][0]

    # Fetch the full protein record in FASTA format
    handle = Entrez.efetch(db="protein", id=protein_id, rettype="fasta", retmode="text")
    fasta_record_str = handle.read()
    handle.close()

    # Parse the FASTA record to get the sequence length
    # Using io.StringIO to treat the string as a file
    with io.StringIO(fasta_record_str) as fasta_file:
      seq_record = SeqIO.read(fasta_file, "fasta")
      protein_length = len(seq_record.seq)

    return protein_length

# The organism is Xanthoria parietina and the protein is XPH1
organism_name = "Xanthoria parietina"
protein_name = "XPH1"

# Get the length
num_amino_acids = get_protein_length(organism_name, protein_name)

# Print the result
if isinstance(num_amino_acids, int):
    print(f"The XPH1 protein of {organism_name} has {num_amino_acids} amino acids.")
else:
    print(num_amino_acids)
