import requests
import io

def solve_task():
    """
    Finds the number of Cysteine residues in the TM3-TM4 linker of the 
    human GABAAρ1 homomeric receptor (UniProt P24046).
    """
    # 1. Identify the UniProt accession ID and define the domain region.
    uniprot_id = "P24046"
    
    # According to UniProt annotation for P24046, the TM3-TM4 linker
    # (large cytoplasmic loop) is located between positions 339 and 430.
    linker_start_pos = 339
    linker_end_pos = 430

    # 2. Fetch the protein sequence from UniProt.
    uniprot_url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
    try:
        response = requests.get(uniprot_url)
        response.raise_for_status()  # Checks for HTTP errors
        
        # Parse the FASTA response to get the clean sequence
        fasta_data = response.text
        full_sequence = "".join(fasta_data.splitlines()[1:])
    except requests.exceptions.RequestException as e:
        print(f"Failed to retrieve data from UniProt. Error: {e}")
        return

    # 3. Extract the linker domain sequence.
    # Python uses 0-based indexing for slicing.
    # The start index is position - 1.
    # The end index for slicing is exclusive, so it's simply the end position.
    start_index = linker_start_pos - 1
    end_index = linker_end_pos
    linker_sequence = full_sequence[start_index:end_index]
    
    # 4. Count the number of Cysteine ('C') residues.
    cysteine_count = linker_sequence.count('C')
    
    # 5. Output the results.
    print(f"Protein: Human GABAAρ1 Receptor (UniProt: {uniprot_id})")
    print(f"Domain of Interest: TM3-TM4 Linker (Residues {linker_start_pos} to {linker_end_pos})")
    print(f"Linker Sequence: {linker_sequence}")
    print(f"\nFinal Calculation:")
    print(f"Number of Cysteine ('C') residues in the linker sequence = {cysteine_count}")

solve_task()