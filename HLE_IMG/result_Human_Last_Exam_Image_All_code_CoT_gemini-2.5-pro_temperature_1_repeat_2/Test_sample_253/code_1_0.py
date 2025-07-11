# The BioPython library is required to run this code.
# If you don't have it installed, you can install it using pip:
# pip install biopython

from Bio import Entrez
from Bio import SeqIO

def get_protein_length(accession_id):
    """
    Fetches a protein record from NCBI and returns its length.
    """
    # NCBI requires an email address for identification
    Entrez.email = "example@example.com"
    
    try:
        # Fetch the protein record from the NCBI 'protein' database in FASTA format
        handle = Entrez.efetch(db="protein", id=accession_id, rettype="fasta", retmode="text")
        
        # Parse the FASTA record
        record = SeqIO.read(handle, "fasta")
        handle.close()
        
        # The length of the sequence is the number of amino acids
        protein_length = len(record.seq)
        
        return protein_length, record.description
        
    except Exception as e:
        return f"An error occurred: {e}", None

# The accession ID for the XPH1 protein from Xanthoria parietina
accession_number = "CAA72420.1"

# Get the protein length and description
length, description = get_protein_length(accession_number)

if isinstance(length, int):
    print(f"Organism and Protein identified: {description}")
    print(f"Using NCBI Accession Number: {accession_number}")
    print(f"The number of amino acids in the protein is the length of its sequence.")
    print(f"Length of protein (accession {accession_number}) = {length} amino acids")
else:
    # Print the error message
    print(length)

# Final answer in the required format
if isinstance(length, int):
    print(f'<<<{length}>>>')
