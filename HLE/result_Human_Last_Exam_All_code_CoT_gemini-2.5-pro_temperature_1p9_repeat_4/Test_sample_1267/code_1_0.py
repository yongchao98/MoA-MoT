def get_mrna_codon(anticodon):
    """
    Calculates the complementary mRNA codon for a given tRNA anticodon.
    Anticodons are typically written 5'-3', but they pair with codons 3'-5'.
    """
    # Reverse the anticodon to get the 3'-5' orientation for pairing
    reversed_anticodon = anticodon[::-1]
    
    # Define base pairing rules
    complementary_bases = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
    
    # Build the mRNA codon
    mrna_codon_list = [complementary_bases.get(base, 'N') for base in reversed_anticodon]
    mrna_codon = "".join(mrna_codon_list)
    
    return mrna_codon

def main():
    """
    Analyzes the effect of the tRNA mutation.
    """
    # Standard genetic code dictionary
    genetic_code = {
        'UUA': 'Leucine', 'UUG': 'Leucine',
        'CAA': 'Glutamine', 'CAG': 'Glutamine',
        # ... and all other codons
    }

    original_anticodon = 'UAA'
    mutated_anticodon = 'UUG'
    
    # 1. Analyze the original tRNA
    original_mrna_codon = get_mrna_codon(original_anticodon)
    original_amino_acid = genetic_code.get(original_mrna_codon, 'Unknown')
    
    print("--- Analysis of the Original tRNA ---")
    print(f"Original anticodon (5'-3'): {original_anticodon}")
    print(f"Binds to mRNA codon (5'-3'): {original_mrna_codon}")
    print(f"Therefore, the tRNA is charged with: {original_amino_acid}\n")
    
    # 2. Analyze the mutated tRNA
    recognized_mrna_codon = get_mrna_codon(mutated_anticodon)
    coded_amino_acid = genetic_code.get(recognized_mrna_codon, 'Unknown')
    
    print("--- Analysis of the Mutated tRNA ---")
    print(f"Mutated anticodon (5'-3'): {mutated_anticodon}")
    print(f"Now binds to mRNA codon (5'-3'): {recognized_mrna_codon}")
    print(f"This codon normally codes for: {coded_amino_acid}\n")
    
    # 3. Determine the consequence
    print("--- Consequence of the Mutation ---")
    print(f"The mutation does not change the amino acid the tRNA carries. It still carries '{original_amino_acid}'.")
    print(f"However, it now delivers '{original_amino_acid}' to the '{recognized_mrna_codon}' codon site.")
    print(f"This results in the misincorporation of '{original_amino_acid}' at protein positions that should contain '{coded_amino_acid}'.")
    print("\nThis means the mutation allows the insertion of an amino acid (Leucine) usually specified by a different codon, which aligns with option C.")

if __name__ == '__main__':
    main()

<<<C>>>