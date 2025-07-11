def get_codon_from_anticodon(anticodon):
    """Calculates the corresponding mRNA codon from a tRNA anticodon."""
    # Reverse the anticodon sequence
    reversed_anticodon = anticodon[::-1]
    # Get the complementary base for each nucleotide
    complement_map = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
    codon = "".join([complement_map[base] for base in reversed_anticodon])
    return codon

def main():
    """
    Analyzes the effect of a tRNA mutation and explains the implications.
    """
    # Genetic code mapping for relevant amino acids
    genetic_code = {
        'UUA': 'Leucine',
        'CAA': 'Glutamine'
    }

    # Define the core sequences of the anticodons
    original_anticodon_core = "UAA"
    mutated_anticodon_core = "UUG"

    # --- Analysis of the original tRNA ---
    original_codon = get_codon_from_anticodon(original_anticodon_core)
    original_amino_acid = genetic_code[original_codon]

    print("--- Original tRNA Function ---")
    print(f"Original Anticodon (5'->3'): {original_anticodon_core}")
    print(f"Recognized mRNA Codon (5'->3'): {original_codon}")
    print(f"Amino Acid Inserted: {original_amino_acid}")
    print("-" * 30)

    # --- Analysis of the mutated tRNA ---
    targeted_codon = get_codon_from_anticodon(mutated_anticodon_core)
    correct_amino_acid_for_target = genetic_code[targeted_codon]
    misincorporated_amino_acid = original_amino_acid # The tRNA is still charged with Leucine

    print("--- Mutated tRNA Function ---")
    print(f"Mutated Anticodon (5'->3'): {mutated_anticodon_core}")
    print(f"Recognized mRNA Codon (5'->3'): {targeted_codon}")
    print(f"This codon should code for: {correct_amino_acid_for_target}")
    print(f"However, the mutated tRNA misincorporates: {misincorporated_amino_acid}")
    print("-" * 30)

    # --- Conclusion ---
    print("\nImplication:")
    print("The mutation causes a tRNA that should carry Leucine to now recognize a codon for Glutamine.")
    print(f"This leads to the occasional misincorporation of {misincorporated_amino_acid} at a protein site where {correct_amino_acid_for_target} should be.")
    print("This means the cell inserts an amino acid (Leucine) that is normally coded for by other codons (like UUA) at a new site (CAA).")
    print("This corresponds to answer choice C.")

if __name__ == "__main__":
    main()