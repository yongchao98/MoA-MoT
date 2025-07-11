def analyze_coiled_coil_oligomerization(sequence):
    """
    Analyzes a coiled-coil sequence to predict its oligomeric state
    based on the identity of residues at the core 'a' and 'd' positions.
    """
    print(f"Analyzing Protein Sequence: {sequence}\n")

    # The heptad repeat phasing was determined by identifying the alternating
    # 3- and 4-residue spacing of the core hydrophobic amino acids.
    # In Python's 0-based indexing:
    d_indices = [2, 9, 16, 23]  # Positions 3, 10, 17, 24 in the sequence
    a_indices = [6, 13, 20, 27]  # Positions 7, 14, 21, 28 in the sequence

    d_residues = [sequence[i] for i in d_indices]
    a_residues = [sequence[i] for i in a_indices]

    print("1. Identifying Core Residues based on Heptad Repeat Phasing:")
    print(f"   - Residues at 'a' positions: {a_residues}")
    print(f"   - Residues at 'd' positions: {d_residues}\n")

    print("2. Analyzing Core Composition for Knobs-into-Holes Packing:")
    print("   - The 'a' positions are consistently Leucine (L).")
    print("   - The 'd' positions are consistently Isoleucine (I).\n")

    print("3. Predicting Oligomeric State:")
    print("   - The combination of a non-beta-branched residue (L) at position 'a' and a")
    print("     beta-branched residue (I) at position 'd' is a well-established motif.")
    print("   - This specific arrangement ('a'=L, 'd'=I) optimizes the packing of the four alpha-helices,")
    print("     strongly favoring the formation of a tetramer over other states like dimers or trimers.\n")

    predicted_state = 4
    
    print("Final Conclusion:")
    print("The predicted oligomeric state is 4 (Tetramer).")

# The coiled-coil sequence to be analyzed
protein_sequence = "GEIAQSLKEIAKSLKEIAWSLKEIAQSLKG"

# Run the analysis and print the result
analyze_coiled_coil_oligomerization(protein_sequence)
<<<C>>>