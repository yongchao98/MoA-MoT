# First, ensure you have Biopython installed:
# pip install biopython

from Bio import Entrez
from Bio import SeqIO
import io

def get_protein_length():
    """
    Searches the NCBI database for the XPH1 protein in Xanthoria parietina
    and returns the number of amino acids in its sequence.
    """
    # Always provide an email to NCBI
    Entrez.email = "user@example.com"

    # Define the search terms
    organism = "Xanthoria parietina"
    protein_name = "XPH1"
    
    print(f"Searching NCBI for the protein '{protein_name}' from the organism '{organism}'...")

    # Construct the search query
    search_term = f'"{organism}"[Organism] AND "{protein_name}"[Protein Name]'

    # Search the protein database
    try:
        handle = Entrez.esearch(db="protein", term=search_term, retmax="1")
        record = Entrez.read(handle)
        handle.close()

        if not record["IdList"]:
            print("Could not find the specified protein.")
            return

        # Get the GI number for the protein
        protein_id = record["IdList"][0]
        print(f"Found protein with NCBI ID: {protein_id}")

        # Fetch the full protein record in FASTA format
        fetch_handle = Entrez.efetch(db="protein", id=protein_id, rettype="fasta", retmode="text")
        fasta_string = fetch_handle.read()
        fetch_handle.close()
        
        # Parse the FASTA record from the string
        seq_record = SeqIO.read(io.StringIO(fasta_string), "fasta")
        
        # The length of the sequence is the number of amino acids
        amino_acid_count = len(seq_record.seq)

        # There's no equation, so we print the final number directly.
        print(f"The number of amino acids in the encoded XPH1 protein is: {amino_acid_count}")

    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please check your internet connection and ensure Biopython is installed.")

if __name__ == "__main__":
    get_protein_length()
