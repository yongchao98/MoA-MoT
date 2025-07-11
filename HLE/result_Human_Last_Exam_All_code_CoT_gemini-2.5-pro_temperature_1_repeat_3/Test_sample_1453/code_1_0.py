def analyze_coiled_coil(sequence):
    """
    Analyzes a coiled-coil sequence to predict its oligomeric state.
    """
    print(f"Analyzing sequence: {sequence}\n")

    # Step 1: Identify the Heptad Register
    # The key to coiled-coils is the (abcdefg)n heptad repeat.
    # We must find the register that places hydrophobic residues at the core 'a' and 'd' positions.
    # A visual inspection of the sequence G-E-I-A-Q-S-L-K-E-I-A-K-S-L-K-E-I-A-W-S-L-K-E-I-A-Q-S-L-K-G
    # reveals a repeating pattern of large hydrophobic residues (I, L, W).
    # A register starting with L at the first 'a' position (index 6) and I at the first 'd' position (index 2)
    # aligns the hydrophobic residues correctly.
    # Register:  b c d e f g a b c d e f g a b c d e f g a b c d e f g a
    # Sequence: G E I A Q S L K E I A K S L K E I A W S L K E I A Q S L K G

    a_indices = [6, 13, 20, 27]
    d_indices = [2, 9, 16, 23]

    a_residues = [sequence[i] for i in a_indices]
    d_residues = [sequence[i] for i in d_indices]

    print("Step 1 & 2: Identify Heptad Register and Extract Core Residues")
    print("----------------------------------------------------------")
    print(f"The most likely heptad register places L/I/W residues in the core.")
    print(f"Residues at 'a' positions: {a_residues}")
    print(f"Residues at 'd' positions: {d_residues}\n")

    # Step 3 & 4: Analyze Packing Efficiency and Impact of Specific Residues
    print("Step 3 & 4: Analyze Packing Efficiency")
    print("---------------------------------------")
    print("We will evaluate the likely oligomeric state based on the size of the core residues.")
    
    # Dimer Analysis
    print("\n*   Analysis for Dimer (2 helices):")
    print("    - Dimer cores are flexible and can accommodate various hydrophobic residues.")
    print("    - The pattern of Leucine (L) at 'a' and Isoleucine (I) at 'd' is a hallmark of dimeric coiled-coils.")
    print(f"    - Our sequence predominantly features L at 'a' ({a_residues}) and I at 'd' ({d_residues}). This strongly favors a dimer.")
    print("    - The single bulky Tryptophan (W) at a 'd' position packs against an 'a' position from the other helix. Packing W against L is sterically demanding but plausible in a flexible dimer.")
    
    # Trimer Analysis
    print("\n*   Analysis for Trimer (3 helices):")
    print("    - Trimer cores are tight. The 'd' positions are buried and require small side chains (e.g., Alanine).")
    print(f"    - Our sequence has large residues (I, W) at the 'd' positions: {d_residues}.")
    print("    - Packing three bulky Isoleucines, and especially three very bulky Tryptophans, into the core at the 'd' positions would cause severe steric clashes.")
    print("    - Therefore, a trimer is highly unlikely.")

    # Tetramer Analysis
    print("\n*   Analysis for Tetramer (4 helices):")
    print("    - Tetramer cores are also very tight and similarly require small residues at the 'd' positions.")
    print(f"    - The large I and W residues at our 'd' positions ({d_residues}) would lead to even greater steric hindrance than in a trimer.")
    print("    - Therefore, a tetramer is also highly unlikely.")

    # Step 5: Conclusion
    print("\nStep 5: Conclusion")
    print("--------------------")
    print("The core residue pattern (L at 'a', I at 'd') is characteristic of a 2-stranded coiled-coil (dimer).")
    print("Higher-order oligomers like trimers and tetramers are destabilized by the presence of large amino acids (I, W) at the sterically constrained 'd' positions.")
    print("The final predicted oligomeric state is 2.")
    
    # Final Answer
    final_answer_value = 2
    final_answer_choice = "A"
    print(f"\nThe oligomeric state is {final_answer_value}, which corresponds to choice {final_answer_choice}.")


protein_sequence = "GEIAQSLKEIAKSLKEIAWSLKEIAQSLKG"
analyze_coiled_coil(protein_sequence)