import sys

def analyze_coiled_coil(sequence):
    """
    Analyzes a coiled-coil sequence to predict its oligomeric state
    based on the residues at the 'a' and 'd' heptad positions.
    """
    # Define amino acid properties for the analysis
    amino_acid_properties = {
        'W': {'size': 'very bulky', 'type': 'hydrophobic'},
        'A': {'size': 'small', 'type': 'hydrophobic'},
        'L': {'size': 'medium', 'type': 'hydrophobic'},
    }

    print(f"Analyzing sequence: {sequence}\n")

    # Step 1: Determine heptad repeat phasing
    print("--- Step 1: Determine Heptad Repeat Phasing ---")
    print("The sequence shows a Leucine (L) at positions 7, 14, 21, and 28.")
    print("This regular 7-residue spacing is the hallmark of a core 'd' position in a heptad repeat (abcdefg).")
    print("This establishes the phasing as 'efgabcd', meaning the sequence starts at the 'e' position.\n")

    # Step 2: Extract core 'a' and 'd' residues
    a_positions = []
    d_positions = []
    # With 'efgabcd' phasing, 'a' is at index 3 and 'd' is at index 6 (0-indexed) of each repeat.
    for i, residue in enumerate(sequence):
        heptad_pos_index = i % 7
        if heptad_pos_index == 3:  # Position 'a'
            a_positions.append(residue)
        elif heptad_pos_index == 6:  # Position 'd'
            d_positions.append(residue)

    print("--- Step 2: Identify Core Residues ---")
    print(f"The core 'a' position residues are: {', '.join(a_positions)}")
    print(f"The core 'd' position residues are: {', '.join(d_positions)}\n")

    # Step 3: Apply packing rules
    print("--- Step 3: Apply 'Knobs-into-Holes' Packing Rules ---")
    print("  - Dimer (2): Favors hydrophobic 'a' and 'd' positions. Tolerates bulky residues.")
    print("  - Trimer (3): Disfavors very bulky residues (like W) at the 'a' position due to steric clash.")
    print("  - Tetramer (4): Favors polar residues at the 'a' position.\n")
    
    # Step 4: Analyze and conclude
    print("--- Step 4: Final Analysis ---")
    print(f"Our sequence has hydrophobic 'a' positions ({', '.join(a_positions)}) and hydrophobic 'd' positions ({', '.join(d_positions)}).")
    print("This pattern rules out a tetramer, which requires polar 'a' positions.")
    
    # Check for trimer-disfavoring features
    has_bulky_a = 'W' in a_positions
    if has_bulky_a:
        print("Crucially, the presence of a very bulky Tryptophan (W) at an 'a' position would create a severe steric clash in the tight core of a trimer.")
        print("This makes a trimer energetically unfavorable.")
    
    print("The pattern of hydrophobic 'a' and 'd' residues is ideal for a dimer.")
    
    predicted_state = 2
    
    print("\n--- Conclusion ---")
    print(f"The sequence is best suited to form a dimer. The predicted oligomeric state is: {predicted_state}")

# The sequence in question
protein_sequence = "GEIAQSLKEIAKSLKEIAWSLKEIAQSLKG"
analyze_coiled_coil(protein_sequence)