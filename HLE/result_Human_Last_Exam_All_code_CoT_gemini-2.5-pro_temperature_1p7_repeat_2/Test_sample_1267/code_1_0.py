def get_mrna_codon(anticodon):
    """Calculates the corresponding mRNA codon for a tRNA anticodon."""
    # Reverse the anticodon to get 3' -> 5' orientation for pairing
    reversed_anticodon = anticodon[::-1]
    
    # Define complementary bases
    complements = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
    
    # Build the mRNA codon
    mrna_codon = ""
    for base in reversed_anticodon:
        mrna_codon += complements.get(base, 'N') # Use 'N' for unknown/modified bases
        
    return mrna_codon

def main():
    """Analyzes the tRNA mutation and its effects."""
    
    # Genetic code mapping (relevant codons)
    genetic_code = {
        'UUA': 'Leucine',
        'CAA': 'Glutamine'
    }

    # Original and mutated anticodons (simplified)
    original_anticodon = "UAA" # Simplified from 5'-xm5s2UAA-3'
    mutated_anticodon = "UUG"   # Simplified from 5'-xm5s2UUG-3'

    # Step 1: Analyze the original tRNA
    original_mrna_codon = get_mrna_codon(original_anticodon)
    original_amino_acid = genetic_code.get(original_mrna_codon)
    print("--- Analysis of Original tRNA ---")
    print(f"Original anticodon: 5'-{original_anticodon}-3'")
    print(f"Binds to mRNA codon: 5'-{original_mrna_codon}-3'")
    print(f"This tRNA is charged with and inserts: {original_amino_acid}\n")

    # Step 2: Analyze the mutated tRNA
    mutated_mrna_codon = get_mrna_codon(mutated_anticodon)
    intended_amino_acid = genetic_code.get(mutated_mrna_codon)
    print("--- Analysis of Mutated tRNA ---")
    print(f"Mutated anticodon: 5'-{mutated_anticodon}-3'")
    print(f"Now binds to mRNA codon: 5'-{mutated_mrna_codon}-3'")
    print(f"This codon normally codes for: {intended_amino_acid}")
    
    # Step 3: Conclude the effect
    print("\n--- Conclusion ---")
    print(f"The mutation does not change the amino acid the tRNA carries, which is still {original_amino_acid}.")
    print("However, it now incorrectly binds to the mRNA codon for a different amino acid.")
    print(f"Result: The mutated tRNA causes {original_amino_acid} to be inserted where {intended_amino_acid} should be.")
    print("This corresponds to a missense mutation in the final protein.")
    print("\nThis allows the insertion of Leucine, an amino acid that is usually brought in by tRNAs with different anticodons, supporting answer choice C.")


if __name__ == "__main__":
    main()
