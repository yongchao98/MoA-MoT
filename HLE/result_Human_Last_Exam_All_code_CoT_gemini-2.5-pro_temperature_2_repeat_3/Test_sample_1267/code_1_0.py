def get_rna_complement(base):
    """Returns the complementary RNA base."""
    if base == 'A':
        return 'U'
    elif base == 'U':
        return 'A'
    elif base == 'C':
        return 'G'
    elif base == 'G':
        return 'C'
    else:
        return 'N' # For unknown/modified bases

def get_codon_from_anticodon(anticodon_seq):
    """
    Calculates the corresponding mRNA codon from a 5'-3' tRNA anticodon sequence.
    The anticodon is reversed to simulate 3'-5' reading and then complemented.
    """
    # Reverse the anticodon sequence to get 3'-5' orientation
    reversed_anticodon = anticodon_seq[::-1]
    
    # Get the complement of each base to find the 5'-3' codon
    codon = "".join([get_rna_complement(base) for base in reversed_anticodon])
    return codon

def main():
    """
    Analyzes the effect of a tRNA mutation and explains the implications.
    """
    # Define the core sequences of the anticodons
    original_anticodon_seq = "UAA"
    mutated_anticodon_seq = "UUG"

    # Define a simplified genetic code dictionary
    genetic_code = {
        'UUA': 'Leucine (Leu)', 'UUG': 'Leucine (Leu)',
        'CUU': 'Leucine (Leu)', 'CUC': 'Leucine (Leu)', 'CUA': 'Leucine (Leu)', 'CUG': 'Leucine (Leu)',
        'CAA': 'Glutamine (Gln)', 'CAG': 'Glutamine (Gln)'
    }

    # --- Step 1: Analyze the original tRNA ---
    print("--- Analysis of the Original tRNA ---")
    original_codon = get_codon_from_anticodon(original_anticodon_seq)
    original_aa = genetic_code.get(original_codon, "Unknown")
    
    print(f"Original Anticodon (5'-3'): {original_anticodon_seq}")
    print(f"Pairs with mRNA Codon (5'-3'): {original_codon}")
    print(f"The original tRNA carries and inserts: {original_aa}\n")
    
    # --- Step 2: Analyze the mutated tRNA ---
    print("--- Analysis of the Mutated tRNA ---")
    mutated_codon_recognized = get_codon_from_anticodon(mutated_anticodon_seq)
    intended_aa_at_codon = genetic_code.get(mutated_codon_recognized, "Unknown")
    
    print(f"Mutated Anticodon (5'-3'): {mutated_anticodon_seq}")
    print(f"Now recognizes mRNA Codon (5'-3'): {mutated_codon_recognized}")
    print(f"This codon normally codes for: {intended_aa_at_codon}\n")
    
    # --- Step 3: Explain the implication ---
    print("--- Implication of the Mutation ---")
    print("The mutation only changed the anticodon, not the tRNA's identity.")
    print(f"Therefore, the mutated tRNA is still charged with its original amino acid, {original_aa}.")
    print(f"During translation, this mutated tRNA now incorrectly binds to the '{mutated_codon_recognized}' codon.")
    print(f"This causes {original_aa} to be inserted into the protein where {intended_aa_at_codon} was supposed to be.")
    print("This event is rare (1 in 1000) because the mutated tRNA must compete with the correct, and more common, tRNA for Glutamine.\n")

    # --- Step 4: Final Conclusion ---
    print("Conclusion: The mutation allows a tRNA carrying Leucine to misread a Glutamine codon,")
    print("causing the insertion of an amino acid (Leucine) that is normally inserted by a different,")
    print("more common tRNA at a different codon. This matches answer choice C.")

if __name__ == '__main__':
    main()
<<<C>>>