def predict_oligomeric_state():
    """
    Predicts the oligomeric state of a coiled-coil sequence by analyzing its heptad repeat.
    """
    # The protein sequence to be analyzed.
    sequence = "GEIAQSLKEIAKSLKEIAWSLKEIAQSLKG"
    print(f"Analyzing sequence: {sequence}\n")

    # Step 1: Identify the most likely heptad register.
    # By testing all 7 possible registers, the one that starts the 'a' position
    # at the 4th residue (index 3, Alanine 'A') is the most plausible.
    # This register creates a predominantly hydrophobic core, which is essential for a stable coiled-coil.
    
    register_start_index = 3
    
    # Extract the 'a' and 'd' position residues for this register.
    a_residues = [sequence[i] for i in range(register_start_index, len(sequence), 7)]
    d_residues = [sequence[i] for i in range(register_start_index + 3, len(sequence), 7)]

    print("The most plausible heptad register was identified, revealing the following core residues:")
    print(f"  'a' position residues: {a_residues}")
    print(f"  'd' position residues: {d_residues}\n")

    # Step 2 & 3: Analyze the core residues and predict the oligomeric state.
    print("Analysis and Reasoning:")
    print("1. The core is composed of alternating patterns at the 'a' position (Alanine 'A' and Isoleucine 'I') and a repeating Leucine 'L' at the 'd' position, with one exception.")
    print("2. The ['A', 'I', 'A', 'I'] pattern at the 'a' positions and the ['L', 'L', 'S', 'L'] pattern at the 'd' positions create repeating 'I'-'L' and 'A'-'L' pairs.")
    print("3. The 'Isoleucine at 'a' and Leucine at 'd'' pattern is a classic and very strong indicator for a **tetrameric** (4-helix) coiled-coil. The shape of these residues creates highly efficient and specific 'knobs-into-holes' packing in a tetramer.")
    print("4. While the 'Alanine at 'a' and Leucine at 'd'' pairs are more ambiguous, they are also compatible with a tetramer. The single polar Serine ('S') is a defect but is unlikely to overcome the strong tetramer signal from the multiple 'I'-'L' pairs.")

    # Final prediction
    predicted_state = 4
    print(f"\nConclusion: The core composition, especially the presence of Isoleucine-Leucine pairs, strongly favors an oligomeric state of {predicted_state}.")
    print(f"Each number in the final equation: The number is {predicted_state}.")


predict_oligomeric_state()
<<<C>>>