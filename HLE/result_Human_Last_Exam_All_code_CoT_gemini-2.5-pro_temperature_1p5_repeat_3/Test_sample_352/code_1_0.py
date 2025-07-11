def suggest_mutagenesis():
    """
    Suggests a site-directed mutagenesis experiment to neutralize a negatively charged
    patch in a protein.
    """
    
    # Define the target residues for mutagenesis
    original_patch = {
        47: {"name": "Serine", "one_letter_code": "S"},
        48: {"name": "Glutamate", "one_letter_code": "E"},
        49: {"name": "Glutamate", "one_letter_code": "E"},
        50: {"name": "Aspartate", "one_letter_code": "D"}
    }
    
    # Define the replacement amino acid
    replacement_amino_acid = {"name": "Alanine", "one_letter_code": "A"}

    print("--- Analysis of the Target Patch ---")
    print("The goal is to eliminate the negative charge from the S-E-E-D patch at positions 47-50.")
    print(f"Position {47} ({original_patch[47]['name']}): A phosphorylation site which, when modified, becomes highly negative.")
    print(f"Position {48} ({original_patch[48]['name']}): Negatively charged.")
    print(f"Position {49} ({original_patch[49]['name']}): Negatively charged.")
    print(f"Position {50} ({original_patch[50]['name']}): Negatively charged.")
    
    print("\n--- Proposed Mutagenesis Strategy ---")
    print(f"The best replacement strategy is to mutate all four residues to {replacement_amino_acid['name']} ({replacement_amino_acid['one_letter_code']}).")
    print("This choice is based on the following reasons:")
    print("1. Neutralizes Charge: Alanine is a neutral amino acid, which will eliminate the negative charge from the patch.")
    print("2. Minimally Disruptive: Alanine is small and chemically inert, reducing the risk of unintended structural changes or interactions.")
    print("3. Standard Practice: Replacing residues with Alanine is a common and effective technique ('Alanine Scanning') to probe the function of specific amino acids.")
    
    print("\n--- Final Recommended Mutation Plan ---")
    
    original_sequence = []
    mutated_sequence = []
    
    # Use a loop to print the mutation for each position
    sorted_positions = sorted(original_patch.keys())
    for position in sorted_positions:
        original_aa = original_patch[position]
        print(f"Position {position}: Mutate {original_aa['name']} ({original_aa['one_letter_code']}) to {replacement_amino_acid['name']} ({replacement_amino_acid['one_letter_code']})")
        original_sequence.append(original_aa['one_letter_code'])
        mutated_sequence.append(replacement_amino_acid['one_letter_code'])

    print(f"\nOriginal Sequence Patch: {'-'.join(original_sequence)}")
    print(f"Proposed Mutant Patch:   {'-'.join(mutated_sequence)}")

# Execute the function to print the plan
suggest_mutagenesis()
<<<A-A-A-A>>>