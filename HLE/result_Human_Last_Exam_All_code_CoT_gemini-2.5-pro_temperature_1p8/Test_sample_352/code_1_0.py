def recommend_mutagenesis():
    """
    Provides a recommendation for a site-directed mutagenesis experiment
    to neutralize a negatively charged patch in protein x.
    """
    
    # Define the original state
    original_patch = {
        47: 'Serine (S)',
        48: 'Glutamate (E)',
        49: 'Glutamate (E)',
        50: 'Aspartate (D)'
    }
    
    # Define the proposed replacement
    replacement_aa = 'Alanine (A)'
    
    # Print the plan and rationale
    print("Plan for Site-Directed Mutagenesis of Protein X")
    print("=" * 50)
    print("Objective: To relieve the autoinhibitory effect of the negatively charged patch at amino acids 47-50.")
    print("\nRationale for Replacement:")
    print(f"The best replacement amino acid is {replacement_aa}. Alanine is small, has a neutral charge,")
    print("and its side chain cannot be phosphorylated. This makes it ideal for neutralizing the")
    print("negative charge from both the acidic residues (E, D) and the phosphoserine (pS)")
    print("while minimizing structural disruption.\n")
    
    print("Proposed Mutations:")
    
    # Print each individual mutation
    for position, original in original_patch.items():
        print(f"  - Position {position}: Replace {original} with {replacement_aa}")
        
    print("\nSummary of Change:")
    original_seq = "".join([val.split(' ')[1][1] for val in original_patch.values()])
    proposed_seq = "A" * len(original_patch)
    print(f"The original sequence patch '{original_seq}' at positions 47-50 should be mutated to '{proposed_seq}'.")
    print("\nThis S47A/E48A/E49A/D50A quadruple mutant will effectively test the hypothesis.")

# Execute the function to print the recommendation
recommend_mutagenesis()