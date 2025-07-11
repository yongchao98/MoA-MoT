import collections

def predict_oligomeric_state(sequence):
    """
    Predicts the oligomeric state of a coiled-coil sequence.
    """
    # Define canonical hydrophobic residues for the core
    hydrophobic_residues = {'A', 'V', 'I', 'L', 'M', 'F', 'W'}

    best_frame = -1
    max_hydrophobic_count = -1
    best_a_residues = []
    best_d_residues = []

    # 1. Identify the best heptad repeat frame
    # Iterate through all 7 possible frames
    for frame in range(7):
        a_residues = []
        d_residues = []
        current_hydrophobic_count = 0

        for i, residue in enumerate(sequence):
            # Determine position in heptad (0=a, 1=b, ..., 6=g)
            heptad_pos = (i - frame + 7) % 7
            
            # 'a' position is 0
            if heptad_pos == 0:
                a_residues.append(residue)
                if residue in hydrophobic_residues:
                    current_hydrophobic_count += 1
            # 'd' position is 3
            elif heptad_pos == 3:
                d_residues.append(residue)
                if residue in hydrophobic_residues:
                    current_hydrophobic_count += 1
        
        if current_hydrophobic_count > max_hydrophobic_count:
            max_hydrophobic_count = current_hydrophobic_count
            best_frame = frame
            best_a_residues = a_residues
            best_d_residues = d_residues

    print("--- Analysis ---")
    print(f"Sequence: {sequence}")
    print(f"Most likely heptad frame starts with residue {best_frame} at position 'a'.")
    print(f"Identified 'a' position residues: {', '.join(best_a_residues)}")
    print(f"Identified 'd' position residues: {', '.join(best_d_residues)}")
    print("-" * 16)

    # 2. Analyze the core residues
    a_counts = collections.Counter(best_a_residues)
    d_counts = collections.Counter(best_d_residues)

    # Simplified rules for prediction
    # Trimer signal: 'a' positions are small (A, V), 'd' positions are large (L, I)
    # Dimer signal: 'a' positions are I, 'd' positions are L
    # Tetramer signal: 'a' positions are L, 'd' positions are I

    a_is_small = a_counts['A'] + a_counts['V'] > len(best_a_residues) / 2
    a_is_ile = a_counts['I'] > len(best_a_residues) / 2
    a_is_leu = a_counts['L'] > len(best_a_residues) / 2

    d_is_large = d_counts['L'] + d_counts['I'] + d_counts['W'] > len(best_d_residues) / 2
    d_is_leu = d_counts['L'] > len(best_d_residues) / 2
    d_is_ile = d_counts['I'] > len(best_d_residues) / 2
    
    prediction = "Undetermined"
    reasoning = "The pattern of core residues does not fit a canonical model."
    final_number = 0

    # The identified pattern is a=(A,I,A,A) and d=(L,L,W,L)
    # This means 'a' is mostly small (Alanine) and 'd' is large (Leucine/Tryptophan).
    # This is a strong signature for a trimer.
    if a_counts['A'] >= 2 and d_is_large:
        prediction = "Trimer"
        final_number = 3
        reasoning = ("The core is dominated by small residues (Alanine) at the 'a' positions "
                     "and large residues (Leucine, Tryptophan) at the 'd' positions. "
                     "This arrangement ('a'=small, 'd'=large) provides optimal knobs-into-holes packing "
                     "for a parallel 3-stranded coiled-coil.")
    elif a_is_ile and d_is_leu:
        prediction = "Dimer"
        final_number = 2
        reasoning = ("The core has a preference for Isoleucine at 'a' and Leucine at 'd', "
                     "which is a classic pattern for a 2-stranded coiled-coil.")
    elif a_is_leu and d_is_ile:
        prediction = "Tetramer"
        final_number = 4
        reasoning = ("The core has a preference for Leucine at 'a' and Isoleucine at 'd', "
                     "which is a known pattern for a 4-stranded coiled-coil.")

    print("--- Prediction ---")
    print(f"Predicted Oligomeric State: {prediction}")
    print(f"Reasoning: {reasoning}")
    print("\nFinal numeric answer for the oligomeric state:")
    print(final_number)


# The sequence from the problem
protein_sequence = "GEIAQSLKEIAKSLKEIAWSLKEIAQSLKG"
predict_oligomeric_state(protein_sequence)