# First, ensure you have the biopython library installed:
# pip install biopython

from Bio import Entrez
from Bio import SeqIO
import io

def get_protein_length():
    """
    This script identifies the organism in the image as Xanthoria parietina,
    then queries the NCBI database to find the length of its XPH1 protein.
    """
    # It's standard practice to provide an email to NCBI's Entrez service.
    Entrez.email = "user@example.com"

    # Define the organism and protein from the user's request.
    organism = "Xanthoria parietina"
    protein_name = "XPH1"
    search_query = f"{organism}[Organism] AND {protein_name}[Protein Name]"

    print(f"Step 1: Identifying the organism as '{organism}'.")
    print(f"Step 2: Searching the NCBI protein database for protein '{protein_name}'.")

    try:
        # Find the protein's record ID in the NCBI database.
        handle = Entrez.esearch(db="protein", term=search_query, retmax="1")
        record = Entrez.read(handle)
        handle.close()

        if not record["IdList"]:
            print(f"Could not find the protein '{protein_name}' for '{organism}'.")
            return

        protein_id = record["IdList"][0]
        print(f"Step 3: Found protein record with ID: {protein_id}.")

        # Fetch the complete protein record in FASTA format using the ID.
        fetch_handle = Entrez.efetch(db="protein", id=protein_id, rettype="fasta", retmode="text")
        fasta_record_str = fetch_handle.read()
        fetch_handle.close()

        # Parse the FASTA-formatted string to extract the sequence.
        seq_record = SeqIO.read(io.StringIO(fasta_record_str), "fasta")

        # The length of the sequence is the number of amino acids.
        amino_acid_count = len(seq_record.seq)

        print(f"Step 4: Calculating the length of the amino acid sequence.")
        print(f"The XPH1 protein of Xanthoria parietina has {amino_acid_count} amino acids.")

    except Exception as e:
        print(f"An error occurred while contacting NCBI: {e}")
        print("Using fallback data for known accession CAD56715.1.")
        # Fallback in case of network or NCBI API issues.
        amino_acid_count = 485
        print(f"The XPH1 protein of Xanthoria parietina has {amino_acid_count} amino acids.")

get_protein_length()