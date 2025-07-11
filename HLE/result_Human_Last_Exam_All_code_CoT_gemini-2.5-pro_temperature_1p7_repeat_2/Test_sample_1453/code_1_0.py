import collections

def predict_oligomeric_state(sequence):
    """
    Analyzes a peptide sequence to predict its coiled-coil oligomeric state.
    This is based on identifying the heptad repeat and analyzing the residues
    at the core 'a' and 'd' positions.
    """
    print(f"Analyzing sequence: {sequence}\n")
    
    hydrophobic_residues = {'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W'}
    register_template = "abcdefg"
    best_frame_info = {'score': -1}

    # 1. Find the optimal heptad repeat frame
    print("--- Step 1: Finding Optimal Heptad Register ---")
    print("Scoring frames by counting hydrophobic residues at core 'a' and 'd' positions...")
    
    # Iterate through all 7 possible start frames (phase shifts)
    for frame_offset in range(7):
        a_pos_residues = []
        d_pos_residues = []
        
        for i, amino_acid in enumerate(sequence):
            register_char = register_template[(i - frame_offset + 7) % 7]
            if register_char == 'a':
                a_pos_residues.append(amino_acid)
            elif register_char == 'd':
                d_pos_residues.append(amino_acid)

        # Calculate score for this frame
        score = sum(1 for r in a_pos_residues if r in hydrophobic_residues) + \
                sum(1 for r in d_pos_residues if r in hydrophobic_residues)
        
        # Keep track of the best frame found so far
        if score > best_frame_info['score']:
            best_frame_info = {
                'score': score,
                'offset': frame_offset,
                'a_residues': a_pos_residues,
                'd_residues': d_pos_residues
            }

    # 2. Display the optimal frame and core residues
    print("\n--- Step 2: Analyzing the Optimal Core ---")
    print(f"Optimal frame found with a score of {best_frame_info['score']} out of a possible {len(best_frame_info['a_residues']) + len(best_frame_info['d_residues'])}.")
    
    full_register = "".join([register_template[(i - best_frame_info['offset'] + 7) % 7] for i in range(len(sequence))])
    
    print("\nSequence aligned to the optimal heptad register:")
    print(f"Register: {full_register}")
    print(f"Sequence: {sequence}")
    
    a_res = best_frame_info['a_residues']
    d_res = best_frame_info['d_residues']
    print(f"\nIdentified core 'a' position residues: {a_res}")
    print(f"Identified core 'd' position residues: {d_res}")

    # 3. Predict the oligomeric state based on 'a' and 'd' residues
    print("\n--- Step 3: Predicting Oligomeric State ---")
    print("Applying 'knobs-into-holes' packing principles:")
    
    # Simplified rules:
    # - Beta-branched (I, V) at 'a' causes steric hindrance, favoring higher oligomers (trimers/tetramers).
    # - Smaller/non-beta-branched (A, L) at 'a' and 'd' fit well in dimers.
    # - The 'Ia-Ld' (Isoleucine at 'a', Leucine at 'd') pattern is a classic tetramer signal.
    # - The 'La-Id' pattern is a classic dimer signal.

    dominant_a = collections.Counter(a_res).most_common(1)[0][0]
    dominant_d = collections.Counter(d_res).most_common(1)[0][0]

    final_state = 0
    
    if dominant_a in ['A', 'L'] and dominant_d in ['A', 'L']:
        print(f"The core is dominated by Alanine ('{dominant_a}') at the 'a' position and Leucine ('{dominant_d}') at the 'd' position.")
        print("This combination of small and non-beta-branched residues packs most efficiently as a DIMER.")
        final_state = 2
    elif dominant_a == 'I' and dominant_d == 'L':
        print("The core has Isoleucine at 'a' and Leucine at 'd', which strongly predicts a TETRAMER.")
        final_state = 4
    elif dominant_a == 'L' and dominant_d == 'I':
        print("The core has Leucine at 'a' and Isoleucine at 'd', which strongly predicts a DIMER.")
        final_state = 2
    else:
        print("The core pattern is not a classic simple case, but the absence of strong trimer/tetramer signals (like Isoleucine at 'a') suggests a lower oligomeric state.")
        final_state = 2 # Defaulting to the most common state, Dimer.

    print("\n-------------------------------------------")
    print(f"Final Predicted Oligomeric State = {final_state}")
    print("-------------------------------------------")

# The protein sequence in question
protein_sequence = "GEIAQSLKEIAKSLKEIAWSLKEIAQSLKG"

# Run the prediction
predict_oligomeric_state(protein_sequence)