def explain_trna_mutation():
    """
    Analyzes a tRNA mutation and explains its impact on protein synthesis.
    """
    
    # A simplified dictionary for the genetic code
    genetic_code = {
        'UUA': 'Leucine (Leu)', 
        'UUG': 'Leucine (Leu)',
        'CAA': 'Glutamine (Gln)',
        'CAG': 'Glutamine (Gln)',
    }

    # Function to find the corresponding mRNA codon for a tRNA anticodon
    def find_mrna_codon(anticodon):
        """
        Calculates the complementary mRNA codon for a 5'-3' tRNA anticodon.
        Pairing is antiparallel, so we reverse and complement.
        """
        # Remove the 5'- and -3' for processing
        clean_anticodon = anticodon.strip("5'-3'")
        
        # Reverse to simulate antiparallel pairing (3'->5' reading)
        reversed_anticodon = clean_anticodon[::-1]
        
        # Complement the bases
        codon = ""
        base_complement = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
        for base in reversed_anticodon:
            codon += base_complement.get(base, 'N') # 'N' for unknown base
            
        return f"5'-{codon}-3'"

    # 1. Define the original and mutated tRNA properties
    original_anticodon = "5'-UAA-3'"
    mutated_anticodon = "5'-UUG-3'"
    carried_amino_acid = "Leucine (Leu)"

    # 2. Analyze the original tRNA function
    original_codon_recognized = find_mrna_codon(original_anticodon)
    original_codon_aa = genetic_code.get(original_codon_recognized.strip("5'-3'"), "Unknown")
    
    print("--- Analysis of the Original tRNA ---")
    print(f"The original tRNA has an anticodon of {original_anticodon}.")
    print(f"This anticodon pairs with the mRNA codon {original_codon_recognized}.")
    print(f"The codon {original_codon_recognized} normally codes for {original_codon_aa}.")
    print(f"Conclusion: The tRNA correctly inserts Leucine at Leucine codons.\n")

    # 3. Analyze the mutated tRNA function
    mutated_codon_recognized = find_mrna_codon(mutated_anticodon)
    mutated_codon_aa = genetic_code.get(mutated_codon_recognized.strip("5'-3'"), "Unknown")
    
    print("--- Analysis of the Mutated tRNA ---")
    print(f"The mutated tRNA has an anticodon of {mutated_anticodon}.")
    print(f"This new anticodon now pairs with the mRNA codon {mutated_codon_recognized}.")
    print(f"The codon {mutated_codon_recognized} is supposed to code for {mutated_codon_aa}.")
    print(f"However, the mutated tRNA is still charged with {carried_amino_acid}.\n")
    
    # 4. Explain the final implication
    print("--- Implication of the Mutation ---")
    print(f"During translation, the mutated tRNA will cause {carried_amino_acid} to be inserted")
    print(f"whenever the ribosome encounters a {mutated_codon_recognized} codon.")
    print(f"This leads to the incorrect substitution of {mutated_codon_aa} with {carried_amino_acid}.")
    print("This matches option C: It allows insertion of an amino acid (Leucine) usually inserted by another, more common anticodon.")

explain_trna_mutation()