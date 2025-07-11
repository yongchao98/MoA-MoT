import collections

def predict_oligomeric_state(sequence):
    """
    Predicts the oligomeric state of a coiled-coil sequence based on the
    residues at the 'a' and 'd' positions of the heptad repeat.
    """
    
    print("Analyzing sequence:", sequence)
    print("-" * 30)

    # Simplified classification: 1 for hydrophobic, 0 for other
    hydrophobic_residues = "AGVILMFPW"
    hydrophobicity = {res: 1 if res in hydrophobic_residues else 0 for res in "ACDEFGHIKLMNPQRSTVWY"}

    best_frame = -1
    max_score = -1
    
    # Step 1: Find the best heptad repeat frame by maximizing hydrophobicity at a/d positions
    print("Step 1: Finding the optimal heptad repeat frame...")
    for frame_start in range(7):
        current_score = 0
        is_bad_frame = False
        # In a valid frame, 'g' position should ideally not be strongly hydrophobic
        # We also check that a/d are consistently hydrophobic
        
        a_and_d_positions = []
        for i in range(frame_start, len(sequence), 7):
            if i < len(sequence):
                a_and_d_positions.append(sequence[i])
            if i + 3 < len(sequence):
                a_and_d_positions.append(sequence[i+3])

        # A good frame has hydrophobic residues at a/d
        score = sum(hydrophobicity.get(res, 0) for res in a_and_d_positions)
        
        # A bad frame has a very hydrophobic residue at 'g'
        g_residues = [sequence[i+6] for i in range(frame_start, len(sequence) - 6, 7)]
        if any(res in 'LIVMW' for res in g_residues):
             score -= 5 # Penalize heavily

        if score > max_score:
            max_score = score
            best_frame = frame_start
    
    print(f"Optimal frame starts with residue {best_frame+1} ('{sequence[best_frame]}') at position 'a'.")
    print("-" * 30)

    # Step 2: Identify the dominant residues at the 'a' and 'd' positions for the best frame
    print("Step 2: Identifying core residues from the optimal frame...")
    a_residues = [sequence[i] for i in range(best_frame, len(sequence), 7)]
    d_residues = [sequence[i+3] for i in range(best_frame, len(sequence) - 3, 7)]

    dominant_a = collections.Counter(a_residues).most_common(1)[0][0]
    dominant_d = collections.Counter(d_residues).most_common(1)[0][0]

    print(f"Residues at core position 'a': {', '.join(a_residues)}")
    print(f"Residues at core position 'd': {', '.join(d_residues)}")
    print(f"This reveals a dominant pattern of '{dominant_a}' at the 'a' positions and '{dominant_d}' at the 'd' positions.")
    print("-" * 30)

    # Step 3: Apply prediction rules
    print("Step 3: Applying prediction rules for oligomeric state...")
    oligomer_state = 0
    reason = "No definitive rule matched."

    if dominant_a == 'I' and dominant_d == 'L':
        oligomer_state = 2
        reason = f"The Ia-Ld pattern ('{dominant_a}' at 'a', '{dominant_d}' at 'd') is characteristic of a dimer."
    elif dominant_a == 'L' and dominant_d == 'I':
        oligomer_state = 3
        reason = f"The La-Id pattern ('{dominant_a}' at 'a', '{dominant_d}' at 'd') strongly favors a trimer for efficient core packing."
    elif dominant_a == 'I' and dominant_d == 'I':
        oligomer_state = 3
        reason = f"The Ia-Id pattern ('{dominant_a}' at 'a', '{dominant_d}' at 'd') favors a trimeric state."
    elif dominant_a == 'L' and dominant_d == 'L':
        oligomer_state = 4
        reason = f"The La-Ld pattern ('{dominant_a}' at 'a', '{dominant_d}' at 'd') often forms a tetramer to relieve steric strain."

    print("Conclusion:")
    print(reason)
    # The final 'equation' as requested
    print(f"The rule '{dominant_a} at a' + '{dominant_d} at d' predicts an oligomeric state of: {oligomer_state}")

# The sequence provided by the user
protein_sequence = "GEIAQSLKEIAKSLKEIAWSLKEIAQSLKG"

# Run the prediction
predict_oligomeric_state(protein_sequence)
<<<B>>>