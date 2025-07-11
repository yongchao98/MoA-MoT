import sys
from io import StringIO

def analyze_trna_mutation():
    """
    Analyzes the effect of a tRNA anticodon mutation.
    """

    # Simplified genetic code dictionary {codon: amino_acid}
    genetic_code = {
        'UUA': 'Leucine', 'UUG': 'Leucine',
        'CAA': 'Glutamine', 'CAG': 'Glutamine'
        # Add other codons as needed for a full simulation
    }

    # Helper to get the corresponding codon from an anticodon
    def get_codon_from_anticodon(anticodon):
        """Calculates the complementary mRNA codon for a given tRNA anticodon."""
        # Reverse the anticodon to align 5'-3' with mRNA's 3'-5'
        reversed_anticodon = anticodon[::-1]
        codon = ""
        for base in reversed_anticodon:
            if base == 'A':
                codon += 'U'
            elif base == 'U':
                codon += 'A'
            elif base == 'C':
                codon += 'G'
            elif base == 'G':
                codon += 'C'
        return codon

    # Original tRNA analysis
    original_anticodon_raw = "5'-xm5s2UAA-3'"
    # Simplify to standard bases for lookup
    original_anticodon = "UAA"
    original_codon = get_codon_from_anticodon(original_anticodon)
    original_amino_acid = genetic_code.get(original_codon, 'Unknown')

    print("--- Original tRNA ---")
    print(f"Anticodon: {original_anticodon_raw}")
    print(f"Recognizes mRNA Codon: 5'-{original_codon}-3'")
    print(f"This codon codes for: {original_amino_acid}")
    print(f"Conclusion: The original tRNA is a {original_amino_acid}-tRNA.")
    print("-" * 25)

    # Mutated tRNA analysis
    mutated_anticodon_raw = "5'-xm5s2UUG-3'"
    # Simplify to standard bases for lookup
    mutated_anticodon = "UUG"
    mutated_codon = get_codon_from_anticodon(mutated_anticodon)
    coded_amino_acid = genetic_code.get(mutated_codon, 'Unknown')

    print("--- Mutated tRNA ---")
    print(f"Mutated Anticodon: {mutated_anticodon_raw}")
    print(f"Now Recognizes mRNA Codon: 5'-{mutated_codon}-3'")
    print(f"This codon normally codes for: {coded_amino_acid}")
    print(f"However, the mutated tRNA is likely still charged with: {original_amino_acid}")
    print("\n--- Implication ---")
    print(f"During translation, at the {mutated_codon} codon, the mutated tRNA will compete with the normal {coded_amino_acid}-tRNA.")
    print(f"This leads to the occasional misincorporation of {original_amino_acid} where {coded_amino_acid} should be.")
    print("This is a missense mutation at the protein level, caused by a misreading tRNA.")
    print("\nThis matches description C: It allows insertion of an amino acid usually inserted by another, more common anticodon.")

# Execute the analysis
analyze_trna_mutation()
