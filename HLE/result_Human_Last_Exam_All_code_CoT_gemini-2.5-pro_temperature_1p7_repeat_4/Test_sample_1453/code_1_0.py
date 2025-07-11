import collections

def predict_oligomeric_state(sequence):
    """
    Analyzes a coiled-coil sequence to predict its oligomeric state.

    Args:
        sequence (str): The amino acid sequence of the protein.
    """

    # Fauchere-Pliska hydrophobicity scale. More positive = more hydrophobic.
    hydrophobicity = {
        'A': 0.31, 'C': 0.77, 'D': -0.77, 'E': -0.67, 'F': 2.02, 'G': 0.07,
        'H': -0.40, 'I': 1.80, 'K': -1.10, 'L': 1.70, 'M': 1.28, 'N': -0.64,
        'P': 0.72, 'Q': -0.69, 'R': -1.81, 'S': -0.18, 'T': 0.00, 'V': 1.22,
        'W': 2.25, 'Y': 1.42
    }
    
    # Residues considered hydrophobic for penalty calculation
    hydrophobic_set = {'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W'}

    analysis_results = []

    print("Step 1: Analyzing all 7 possible heptad repeat frames.")
    for offset in range(7):
        core_a = [sequence[i] for i in range(len(sequence)) if (i + offset) % 7 == 0]
        core_d = [sequence[i] for i in range(len(sequence)) if (i + offset) % 7 == 3]
        
        exposed_e = [sequence[i] for i in range(len(sequence)) if (i + offset) % 7 == 4]
        exposed_g = [sequence[i] for i in range(len(sequence)) if (i + offset) % 7 == 6]
        
        # Calculate core hydrophobicity score
        core_score = sum(hydrophobicity.get(aa, 0) for aa in core_a) + \
                     sum(hydrophobicity.get(aa, 0) for aa in core_d)

        # Calculate exposure penalty (sum of hydrophobicity for hydrophobic residues at e/g)
        exposure_penalty = sum(hydrophobicity.get(aa, 0) for aa in exposed_e if aa in hydrophobic_set) + \
                           sum(hydrophobicity.get(aa, 0) for aa in exposed_g if aa in hydrophobic_set)
                           
        analysis_results.append({
            'offset': offset,
            'a': core_a,
            'd': core_d,
            'core_score': core_score,
            'exposure_penalty': exposure_penalty
        })

    # Sort results first by highest core_score, then by lowest exposure_penalty
    sorted_results = sorted(analysis_results, key=lambda x: (-x['core_score'], x['exposure_penalty']))
    
    best_frame = sorted_results[0]
    
    print("\nStep 2: Selecting the best frame based on core hydrophobicity and surface exposure.")
    print(f"The best frame found has offset {best_frame['offset']}.")
    print(f"  - Core Score (hydrophobicity of 'a'/'d' residues): {best_frame['core_score']:.2f} (Maximized)")
    print(f"  - Exposure Penalty (hydrophobicity of exposed 'e'/'g' residues): {best_frame['exposure_penalty']:.2f} (Minimized)")
    
    best_a = best_frame['a']
    best_d = best_frame['d']

    print("\nStep 3: Analyzing the core residue pattern of the best frame.")
    print(f"  - Residues at 'a' positions: {', '.join(best_a)}")
    print(f"  - Residues at 'd' positions: {', '.join(best_d)}")

    # Apply prediction rules
    prediction = "Unknown"
    reason = "The core residue pattern does not match a canonical motif."

    # Tetramer rule: a=L, d=I
    if all(res == 'L' for res in best_a) and all(res == 'I' for res in best_d):
        prediction = 4
        reason = "The pattern of non-branched Leucine (L) at 'a' and beta-branched Isoleucine (I) at 'd' is a canonical motif for forming a 4-helix bundle (tetramer)."
    # Dimer rule: d=L, a=small (like A)
    elif all(res == 'L' for res in best_d) and all(res == 'A' for res in best_a):
        prediction = 2
        reason = "The pattern of Leucine (L) at 'd' and the small Alanine (A) at 'a' is a common motif for forming a 2-helix coiled-coil (dimer)."
    # Trimer rule: a=I, d=L
    elif all(res == 'I' for res in best_a) and all(res == 'L' for res in best_d):
        prediction = 3
        reason = "The pattern of beta-branched Isoleucine (I) at 'a' and non-branched Leucine (L) at 'd' is a canonical motif for forming a 3-helix coiled-coil (trimer)."

    print("\nStep 4: Predicting the oligomeric state based on established packing rules.")
    print(f"Conclusion: {reason}")
    print("\nFinal predicted oligomeric state:")
    print(prediction)

if __name__ == '__main__':
    protein_sequence = "GEIAQSLKEIAKSLKEIAWSLKEIAQSLKG"
    predict_oligomeric_state(protein_sequence)
<<<C>>>