import requests
import re

def count_cysteines_in_linker():
    """
    Finds the number of Cysteine residues in the TM3-TM4 linker of the
    human GABAAρ1 homomeric receptor (UniProt P24046).
    """
    uniprot_id = "P24046"
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.txt"

    print(f"Fetching protein data for UniProt ID: {uniprot_id}...")
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()  # Raise an exception for bad status codes
        data = response.text
    except requests.exceptions.RequestException as e:
        print(f"Error fetching data from UniProt: {e}")
        return

    tm_helices = []
    full_sequence = ""
    in_sequence_block = False

    # Parse the UniProt text file for TRANSMEM annotations and sequence
    for line in data.splitlines():
        if line.startswith("FT   TRANSMEM"):
            # Use regex to find start and end positions of the domain
            match = re.search(r'(\d+)\.\.(\d+)', line)
            if match:
                start, end = map(int, match.groups())
                tm_helices.append({'start': start, 'end': end})

        # The sequence block starts with 'SQ'
        if line.startswith("SQ   "):
            in_sequence_block = True
            continue
        # The sequence block ends with '//'
        if line.startswith("//"):
            in_sequence_block = False
            continue
        # If we are inside the sequence block, build the sequence string
        if in_sequence_block:
            # Remove spaces and line breaks
            full_sequence += "".join(line.strip().split())

    if len(tm_helices) < 4:
        print("Error: Could not find at least 4 transmembrane domains in the UniProt data.")
        return
    
    if not full_sequence:
        print("Error: Could not extract the protein sequence.")
        return

    # UniProt lists domains sequentially, so TM3 is the 3rd and TM4 is the 4th
    tm3_end = tm_helices[2]['end']
    tm4_start = tm_helices[3]['start']

    # The linker is the region between the end of TM3 and the start of TM4.
    # We use 1-based indexing for reporting, consistent with UniProt.
    linker_start_pos = tm3_end + 1
    linker_end_pos = tm4_start - 1

    # Python string slicing is 0-based.
    # To get from linker_start_pos to linker_end_pos (inclusive),
    # we slice from index (linker_start_pos - 1) up to (linker_end_pos).
    linker_sequence = full_sequence[linker_start_pos - 1 : linker_end_pos]

    # Count the number of Cysteine ('C') residues
    cysteine_count = linker_sequence.count('C')

    print("\n--- Analysis Results ---")
    print(f"TM3 Location: Residues {tm_helices[2]['start']} to {tm_helices[2]['end']}")
    print(f"TM4 Location: Residues {tm_helices[3]['start']} to {tm_helices[3]['end']}")
    print(f"TM3-TM4 Linker is therefore from residue {linker_start_pos} to {linker_end_pos}.")
    print(f"The length of the linker sequence is {len(linker_sequence)} amino acids.")
    print("\nCounting the Cysteine ('C') residues in the linker...")
    print(f"Number of Cysteine residues found = {cysteine_count}")


if __name__ == "__main__":
    count_cysteines_in_linker()