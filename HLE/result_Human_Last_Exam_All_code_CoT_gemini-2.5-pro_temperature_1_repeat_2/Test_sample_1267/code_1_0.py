def get_codon_from_anticodon(anticodon):
    """
    Calculates the corresponding mRNA codon for a given tRNA anticodon.
    Note: This is a simplified model and does not account for modified bases or wobble pairing variations.
    """
    # Reverse the anticodon to align 3' -> 5' for pairing
    reversed_anticodon = anticodon[::-1]
    
    # Define the base pairing rules
    pairing_rules = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
    
    # Determine the complementary codon
    codon = ""
    for base in reversed_anticodon:
        codon += pairing_rules.get(base, '?')
        
    return codon

def main():
    """
    Analyzes the effect of a tRNA anticodon mutation.
    """
    genetic_code = {
        'UUA': 'Leucine (Leu)', 'UUG': 'Leucine (Leu)',
        'CUU': 'Leucine (Leu)', 'CUC': 'Leucine (Leu)', 'CUA': 'Leucine (Leu)', 'CUG': 'Leucine (Leu)',
        'CAA': 'Glutamine (Gln)', 'CAG': 'Glutamine (Gln)',
        # Add other codons as needed for a full map
    }

    # Define the anticodons (using standard bases for simplicity)
    original_anticodon_core = "UAA"
    mutated_anticodon_core = "UUG"

    # --- Analysis of the Original tRNA ---
    print("--- Original tRNA Analysis ---")
    original_codon = get_codon_from_anticodon(original_anticodon_core)
    original_amino_acid = genetic_code.get(original_codon, "Unknown")
    print(f"Original Anticodon (5'->3'): {original_anticodon_core}")
    print(f"Recognized mRNA Codon (5'->3'): {original_codon}")
    print(f"Amino Acid Normally Inserted: {original_amino_acid}")
    print(f"Conclusion: The original tRNA is a tRNA-{original_amino_acid.split(' ')[0]}.")
    print("-" * 30)

    # --- Analysis of the Mutated tRNA ---
    print("--- Mutated tRNA Analysis ---")
    mutated_codon = get_codon_from_anticodon(mutated_anticodon_core)
    intended_amino_acid = genetic_code.get(mutated_codon, "Unknown")
    print(f"Mutated Anticodon (5'->3'): {mutated_anticodon_core}")
    print(f"NEW Recognized mRNA Codon (5'->3'): {mutated_codon}")
    print(f"This codon normally codes for: {intended_amino_acid}")
    print("\n--- Effect on Translation ---")
    print("The tRNA's charging identity (which amino acid it carries) is determined by its structure, not the anticodon.")
    print(f"Therefore, the mutated tRNA still carries: {original_amino_acid}")
    print(f"During translation, when the ribosome encounters a {mutated_codon} codon:")
    print(f"1. The correct tRNA for {intended_amino_acid} should bind.")
    print(f"2. However, the mutated tRNA can now also bind to the {mutated_codon} codon.")
    print(f"3. If the mutated tRNA binds, it incorrectly inserts {original_amino_acid} instead of {intended_amino_acid}.")
    print("\nThis means the mutation allows the insertion of Leucine (an amino acid normally inserted by tRNA-Leu) at a site for Glutamine.")
    print("This perfectly matches answer choice C.")


if __name__ == "__main__":
    main()
