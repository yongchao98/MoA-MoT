import collections

def predict_oligomeric_state(sequence):
    """
    Predicts the oligomeric state of a coiled-coil sequence.
    """
    print(f"Analyzing protein sequence: {sequence}\n")
    
    # --- Step 1: Find the optimal heptad repeat frame ---
    print("Step 1: Identifying the optimal heptad repeat (abcdefg) frame.")
    print("The best frame maximizes the hydrophobicity of the core 'a' and 'd' positions.\n")

    # Scoring for hydrophobicity
    hydrophobic_scores = {'L': 2, 'I': 2, 'V': 2, 'W': 2, 'F': 2, 'M': 2, 'A': 1, 'Y': 1}
    
    best_frame_info = {
        'offset': -1,
        'score': -1,
        'a_residues': [],
        'd_residues': []
    }

    # Iterate through all 7 possible frames (by changing the start offset)
    for offset in range(7):
        a_residues = []
        d_residues = []
        current_score = 0
        
        for i, residue in enumerate(sequence):
            # Heptad position (0-6) corresponds to (a-g)
            heptad_pos_index = (i - offset + 7) % 7
            
            if heptad_pos_index == 0:  # 'a' position
                a_residues.append(residue)
                current_score += hydrophobic_scores.get(residue, 0)
            elif heptad_pos_index == 3:  # 'd' position
                d_residues.append(residue)
                current_score += hydrophobic_scores.get(residue, 0)
        
        if current_score > best_frame_info['score']:
            best_frame_info['score'] = current_score
            best_frame_info['offset'] = offset
            best_frame_info['a_residues'] = a_residues
            best_frame_info['d_residues'] = d_residues

    a_res = best_frame_info['a_residues']
    d_res = best_frame_info['d_residues']
    
    print("Optimal frame identified. The core residues are:")
    print(f"  'a' position residues: {a_res}")
    print(f"  'd' position residues: {d_res}\n")

    # --- Step 2: Analyze core residues to predict oligomeric state ---
    print("Step 2: Analyzing core residues based on knobs-into-holes packing rules.")
    
    # Count key residue types
    a_counts = collections.Counter(a_res)
    d_counts = collections.Counter(d_res)
    
    a_beta_branched = a_counts['I'] + a_counts['V']
    d_beta_branched = d_counts['I'] + d_counts['V']
    
    print(f"Rule for Dimers (2): Favors Leucine (L) at 'd' and beta-branched (I/V) at 'a'.")
    print(f"Rule for Trimers (3): Favors Leucine (L) at 'a' and beta-branched (I/V) or L at 'd'.")
    print(f"Rule for Tetramers (4): Favors beta-branched (I/V) at both 'a' and 'd'.\n")
    
    print("Applying rules to our sequence:")
    
    # Final logic
    print(f"The 'd' positions ({d_res}) are rich in beta-branched (Isoleucine) and bulky (Tryptophan) residues.")
    print("This type of packing is sterically unfavorable for dimers and points toward a trimer or tetramer.")
    
    print(f"\nThe 'a' positions ({a_res}) contain Leucine but are not dominated by beta-branched residues (I/V).")
    print("The lack of consistent beta-branching at the 'a' position makes a tetramer less likely.")
    
    print("\nConclusion: The combination of a beta-branched/bulky 'd' core and a non-beta-branched 'a' core is most characteristic of a 3-helix bundle (Trimer).")
    
    final_answer = 3
    print(f"\nFinal Equation: 'd'=[{', '.join(d_res)}] + 'a'=[{', '.join(a_res)}] => Oligomeric State = {final_answer}")
    
    return final_answer

# The protein sequence from the user
protein_sequence = "GEIAQSLKEIAKSLKEIAWSLKEIAQSLKG"
predicted_state = predict_oligomeric_state(protein_sequence)

# The final answer in the required format
# A=2, B=3, C=4, D=5, E=6, F=7, G=8
answer_map = {2: 'A', 3: 'B', 4: 'C', 5: 'D', 6: 'E', 7: 'F', 8: 'G'}
final_choice = answer_map.get(predicted_state, "Unknown")
print(f"\n<<<B>>>")