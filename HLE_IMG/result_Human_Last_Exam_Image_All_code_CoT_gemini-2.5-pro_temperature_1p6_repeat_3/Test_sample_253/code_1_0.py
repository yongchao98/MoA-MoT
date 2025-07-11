# First, ensure you have biopython installed:
# pip install biopython

from Bio import Entrez
from Bio import SeqIO
import io

def get_protein_length(organism, protein_name):
    """
    Fetches a protein sequence from NCBI and returns its length.
    """
    try:
        # Provide your email to NCBI (required)
        Entrez.email = "your.email@example.com"

        # Search for the protein in the NCBI protein database
        search_term = f'"{organism}"[Organism] AND "{protein_name}"[Protein Name]'
        handle = Entrez.esearch(db="protein", term=search_term, retmax="1")
        record = Entrez.read(handle)
        handle.close()

        if not record["IdList"]:
            print(f"Protein '{protein_name}' not found for organism '{organism}'.")
            return None

        protein_id = record["IdList"][0]

        # Fetch the protein record in FASTA format
        fetch_handle = Entrez.efetch(db="protein", id=protein_id, rettype="fasta", retmode="text")
        fasta_record_str = fetch_handle.read()
        fetch_handle.close()

        # Use a string IO object to parse the record from the string
        fasta_record = SeqIO.read(io.StringIO(fasta_record_str), "fasta")
        
        # The length of the sequence is the number of amino acids
        num_amino_acids = len(fasta_record.seq)
        
        return num_amino_acids
        
    except Exception as e:
        print(f"An error occurred: {e}")
        return None

# The organism in the image is likely Xanthoria parietina
organism = "Xanthoria parietina"
protein_name = "XPH1"

# Get the number of amino acids
amino_acid_count = get_protein_length(organism, protein_name)

if amino_acid_count is not None:
    print(f"The XPH1 protein of the organism {organism} has the following number of amino acids:")
    # The prompt requires printing the number in the final equation.
    # In this context, the final "equation" is the statement of the result.
    print(f"Number of amino acids = {amino_acid_count}")
