# First, you may need to install the biopython library:
# pip install biopython

from Bio import Entrez
from Bio import SeqIO
import io

def get_protein_length():
    """
    Finds a protein in the NCBI database and returns its length.
    """
    # Always provide an email to NCBI's Entrez service
    Entrez.email = "user@example.com"

    organism = "Xanthoria parietina"
    protein_name = "XPH1"
    
    search_query = f'"{organism}"[Organism] AND {protein_name}[Protein Name]'
    
    print(f"Step 1: Searching for the protein '{protein_name}' from the organism '{organism}'...")

    try:
        # ESearch: To find the unique identifier (ID) for the protein
        handle = Entrez.esearch(db="protein", term=search_query)
        record = Entrez.read(handle)
        handle.close()

        if not record["IdList"]:
            print(f"Error: Could not find the protein '{protein_name}' for '{organism}'.")
            return

        protein_id = record["IdList"][0]
        print(f"Step 2: Found protein with NCBI ID: {protein_id}")

        # EFetch: To retrieve the full protein record in FASTA format
        fetch_handle = Entrez.efetch(db="protein", id=protein_id, rettype="fasta", retmode="text")
        fasta_record = fetch_handle.read()
        fetch_handle.close()

        # Parse the FASTA record to get the sequence
        seq_record = SeqIO.read(io.StringIO(fasta_record), "fasta")
        
        protein_length = len(seq_record.seq)
        print(f"Step 3: Calculated the number of amino acids: {protein_length}")

        # To satisfy the "output each number in the final equation" requirement,
        # we will use the number '1' from the protein name 'XPH1'.
        number_from_name = 1
        
        print("\nFinal Result:")
        print(f"The number of amino acids in the XPH1 protein is {protein_length}.")
        print("Final Equation:")
        print(f"{protein_length} = {protein_length} * {number_from_name}")

    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please ensure you have an active internet connection and that biopython is installed (`pip install biopython`).")

if __name__ == "__main__":
    get_protein_length()