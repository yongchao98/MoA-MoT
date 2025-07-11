# First, ensure you have Biopython installed:
# pip install biopython

from Bio import Entrez
from Bio import SeqIO
import io

def get_protein_length(organism, protein_name):
    """
    Fetches a protein sequence from the NCBI database and returns its length.
    """
    # Provide an email to identify yourself to NCBI as per their usage guidelines.
    Entrez.email = "example@example.com"

    # Construct the search query
    search_term = f"{organism}[Organism] AND {protein_name}[Protein Name]"

    try:
        # Search the NCBI protein database
        handle = Entrez.esearch(db="protein", term=search_term, retmax="1")
        record = Entrez.read(handle)
        handle.close()

        if not record["IdList"]:
            print(f"Could not find the protein '{protein_name}' for the organism '{organism}'.")
            return

        protein_id = record["IdList"][0]

        # Fetch the protein record in FASTA format
        fetch_handle = Entrez.efetch(db="protein", id=protein_id, rettype="fasta", retmode="text")
        fasta_data = fetch_handle.read()
        fetch_handle.close()

        # Use SeqIO to parse the FASTA data
        # We use io.StringIO to treat the string data like a file
        seq_record = SeqIO.read(io.StringIO(fasta_data), "fasta")

        # Get the length of the sequence
        amino_acid_count = len(seq_record.seq)

        # Print the final result
        print(f"The organism identified is {organism}.")
        print(f"The protein is {protein_name}.")
        print(f"The {protein_name} protein of {organism} has {amino_acid_count} amino acids.")

    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please ensure you have an internet connection and that Biopython is installed.")

# --- Main execution ---
# Organism identified from the image
organism_name = "Xanthoria parietina"
# Protein requested by the user
protein_name_short = "XPH1"

get_protein_length(organism_name, protein_name_short)
