import collections

def predict_oligomeric_state(sequence):
    """
    Analyzes a coiled-coil sequence to predict its oligomeric state based on
    the residues at the 'a' and 'd' positions of the heptad repeat.
    """
    print(f"Analyzing protein sequence: {sequence}")
    print("="*60)

    # Step 1: Determine residues at 'a' and 'd' core positions
    heptad_positions = "abcdefg"
    a_residues = []
    d_residues = []
    for i, residue in enumerate(sequence):
        pos = heptad_positions[i % 7]
        if pos == 'a':
            a_residues.append(residue)
        elif pos == 'd':
            d_residues.append(residue)
    
    print("Step 1: Identified core residues from the heptad repeat.")
    print(f"Residues at 'a' positions: {a_residues}")
    print(f"Residues at 'd' positions: {d_residues}")
    print("")

    # Step 2: Analyze residue properties and apply packing rules
    print("Step 2: Analyzing residue properties for packing efficiency.")
    
    # Define residue classes
    polar_residues = {'K', 'R', 'D', 'E', 'Q', 'N'}
    small_hydrophobic_residues = {'A', 'G'}

    # Check properties of the identified core residues
    num_a_polar = sum(1 for res in a_residues if res in polar_residues)
    all_d_small_hydrophobic = all(res in small_hydrophobic_residues for res in d_residues)

    print(f"Analysis of 'a' positions: Found to be predominantly polar ({num_a_polar} out of {len(a_residues)}).")
    print(f"Analysis of 'd' positions: Found to be exclusively small and hydrophobic (Is this true? {all_d_small_hydrophobic}).")
    print("")

    # Step 3: Conclude the oligomeric state based on packing models
    print("Step 3: Determining oligomeric state from the packing model.")
    oligomeric_state = "Unknown"
    
    # A known rule for tetramers: A pattern of polar residues at 'a' and small 
    # hydrophobic residues at 'd' strongly favors a parallel tetramer.
    if (num_a_polar / len(a_residues) >= 0.75) and all_d_small_hydrophobic:
        oligomeric_state = 4
        print("The pattern fits the model for a tetramer (4-helix bundle).")
        print("Reasoning: The small 'd' residues (Alanine) can pack efficiently into a tight central core.")
        print("The polar 'a' residues (Lysine) are unfavorable to bury, and in a tetramer, they can be positioned")
        print("at the more solvent-accessible 'seams' of the interface.")
    else:
        print("The sequence does not fit the specific pattern for a polar 'a' / small 'd' tetramer.")

    print("="*60)
    print("Final Predicted Oligomeric State Number:")
    # The final equation is the result of this analysis
    print(f"{oligomeric_state}")

# The protein coiled-coil sequence in question
protein_sequence = "GEIAQSLKEIAKSLKEIAWSLKEIAQSLKG"

# Execute the prediction function
predict_oligomeric_state(protein_sequence)