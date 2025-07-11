def design_mutagenesis_experiment():
    """
    Designs and explains the site-directed mutagenesis experiment to relieve
    a negatively charged patch in protein x.
    """

    # Define the original amino acid patch
    original_positions = [47, 48, 49, 50]
    original_amino_acids = ['S', 'E', 'E', 'D']
    
    # Define the proposed replacement amino acid
    replacement_amino_acid = 'A'
    proposed_amino_acids = [replacement_amino_acid] * 4

    # Explain the rationale
    print("--- Site-Directed Mutagenesis Experimental Design ---")
    print("\nObjective: To relieve the autoinhibitory, negatively charged patch (aa 47-50).\n")
    print("Strategy: Replace the native amino acids with Alanine (A), a small, neutral, and non-phosphorylatable residue.\n")
    
    # --- Display the changes ---
    print("--- Proposed Amino Acid Changes ---")
    print(f"{'Position':<10} | {'Original':<10} | {'Proposed':<10} | {'Justification'}")
    print("-" * 70)
    
    # Position 47
    print(f"{original_positions[0]:<10} | {original_amino_acids[0]} (Serine)   | {proposed_amino_acids[0]} (Alanine)  | Prevents phosphorylation, removing potential negative charge.")
    
    # Position 48
    print(f"{original_positions[1]:<10} | {original_amino_acids[1]} (Glutamate)| {proposed_amino_acids[1]} (Alanine)  | Replaces a negatively charged residue with a neutral one.")

    # Position 49
    print(f"{original_positions[2]:<10} | {original_amino_acids[2]} (Glutamate)| {proposed_amino_acids[2]} (Alanine)  | Replaces a negatively charged residue with a neutral one.")

    # Position 50
    print(f"{original_positions[3]:<10} | {original_amino_acids[3]} (Aspartate)| {proposed_amino_acids[3]} (Alanine)  | Replaces a negatively charged residue with a neutral one.")
    
    # --- Display the final sequence comparison ---
    print("\n--- Sequence Summary ---")
    
    original_sequence_str = '-'.join(original_amino_acids)
    proposed_sequence_str = '-'.join(proposed_amino_acids)

    # We still need to output the original numbers and characters for the equation part of the prompt
    print(f"The original sequence at positions 47-48-49-50 is {original_amino_acids[0]}-{original_amino_acids[1]}-{original_amino_acids[2]}-{original_amino_acids[3]}.")
    print(f"The best replacement sequence at positions 47-48-49-50 would be {proposed_amino_acids[0]}-{proposed_amino_acids[1]}-{proposed_amino_acids[2]}-{proposed_amino_acids[3]}.")

    # Final Answer in specified format
    final_answer = 'A-A-A-A'
    print(f"\nFinal recommended replacement sequence for positions 47-50 is {final_answer}")

# Run the function to display the plan
design_mutagenesis_experiment()
<<<A-A-A-A>>>