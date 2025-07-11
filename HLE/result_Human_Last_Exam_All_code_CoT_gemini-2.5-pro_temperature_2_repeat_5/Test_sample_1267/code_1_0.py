def get_mrna_codon(anticodon):
    """
    Calculates the corresponding mRNA codon for a given tRNA anticodon.
    Anticodons are read 3'->5' on the mRNA, so we need the reverse complement.
    """
    complement = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
    # Assumes input is a standard 5'->3' anticodon sequence
    # No need to reverse the input string as Python slicing handles it
    complementary_strand = ""
    for base in anticodon:
        complementary_strand += complement.get(base, 'N')
    # Reverse the complementary strand to get the 5'->3' mRNA codon
    mrna_codon = complementary_strand[::-1]
    return mrna_codon

def main():
    """
    Analyzes the effect of the tRNA mutation.
    """
    # Define a simplified genetic code for the relevant amino acids
    genetic_code = {
        'UUA': 'Leucine (Leu)',
        'CAA': 'Glutamine (Gln)',
    }

    # Define the anticodons (ignoring modifications for sequence pairing)
    original_anticodon_seq = "UAA"
    mutated_anticodon_seq = "UUG"

    # Step 1: Analyze the original tRNA
    original_mrna_codon = get_mrna_codon(original_anticodon_seq)
    original_amino_acid = genetic_code.get(original_mrna_codon, 'Unknown')
    print("--- Original tRNA Function ---")
    print(f"Original Anticodon (5'->3'): {original_anticodon_seq}")
    print(f"Recognized mRNA Codon (5'->3'): {original_mrna_codon}")
    print(f"Inserts Amino Acid: {original_amino_acid}")
    print("-" * 30)

    # Step 2: Analyze the mutated tRNA
    mutated_mrna_target = get_mrna_codon(mutated_anticodon_seq)
    print("--- Mutated tRNA Function ---")
    print(f"Mutated Anticodon (5'->3'): {mutated_anticodon_seq}")
    print(f"New Recognized mRNA Codon (5'->3'): {mutated_mrna_target}")
    print("-" * 30)
    
    # Step 3: Analyze the consequence of the mutation
    normal_amino_acid_for_target = genetic_code.get(mutated_mrna_target, 'Unknown')
    print("--- Consequence of Mutation ---")
    print(f"The mutated tRNA now reads the '{mutated_mrna_target}' codon, which normally codes for {normal_amino_acid_for_target}.")
    print(f"The tRNA synthetase likely still charges the mutated tRNA with {original_amino_acid}.")
    print(f"Therefore, the mutated tRNA competes with the normal {normal_amino_acid_for_target}-tRNA at the '{mutated_mrna_target}' codon site.")
    print("\n--- Final Implication ---")
    print(f"This leads to the occasional (1 in 1000) mis-insertion of {original_amino_acid} where {normal_amino_acid_for_target} should be.")
    print("This perfectly matches answer choice C: 'It allows insertion of an amino acid usually inserted by another, more common anticodon.'")


if __name__ == "__main__":
    main()

<<<C>>>