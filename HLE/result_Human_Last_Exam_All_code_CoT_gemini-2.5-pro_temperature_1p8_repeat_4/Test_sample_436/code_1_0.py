import pandas as pd

def rank_epitopes():
    """
    Ranks peptide epitopes based on their predicted binding affinity to H2-Kd.
    """
    # Define the epitopes
    epitopes = {
        'E1': 'TYQRTRALV', # High-affinity reference
        'E2': 'TFQRTRALV', # Sub-optimal P2 anchor
        'E3': 'TFQRTRALK', # Very poor P9 anchor
        'E4': 'TYQRMFALV', # Optimal anchors, sub-optimal core
        'E5': 'TYIPMRALV', # Optimal anchors, proline in core
    }

    # Simplified scoring matrix based on H2-Kd binding motifs
    # The scores represent relative binding contribution.
    p2_scores = {'Y': 10, 'F': 7, 'default': 0}
    p9_scores = {'V': 9, 'L': 10, 'K': 1, 'default': 0}

    # Penalties for core sequence deviations from the optimal E1
    # Reference core (P3-P7) is QRTRA
    reference_core = list(epitopes['E1'][2:7])

    results = []

    for name, seq in epitopes.items():
        # Get anchor residues
        p2 = seq[1]
        p9 = seq[8]

        # Calculate anchor scores
        p2_score = p2_scores.get(p2, p2_scores['default'])
        p9_score = p9_scores.get(p9, p9_scores['default'])
        anchor_score = p2_score + p9_score
        
        # Calculate core penalty
        core_penalty = 0
        current_core = list(seq[2:7])
        
        # Check for Proline in the core (P3-P7, index 2-6)
        if 'P' in current_core:
            core_penalty -= 5 # Significant penalty for proline's structural disruption
        
        # Penalize for other substitutions from the reference core
        subs = sum(1 for i in range(len(reference_core)) if current_core[i] != reference_core[i] and current_core[i] != 'P')
        core_penalty -= subs # Minor penalty for each substitution

        total_score = anchor_score + core_penalty
        
        results.append({
            'Epitope': name,
            'Sequence': seq,
            'P2_Score': p2_score,
            'P9_Score': p9_score,
            'Core_Penalty': core_penalty,
            'Total_Score': total_score,
            'Reasoning': []
        })

    # Sort by total score in descending order
    ranked_epitopes = sorted(results, key=lambda x: x['Total_Score'], reverse=True)

    # Build detailed reasoning for the printout
    print("Ranking epitopes based on predicted H2-Kd binding affinity:")
    print("=" * 60)
    print(f"{'Rank':<5} {'Epitope':<8} {'Sequence':<12} {'Score':<6} {'Reasoning'}")
    print("-" * 60)
    for i, res in enumerate(ranked_epitopes):
        reason = []
        if res['Epitope'] == 'E1': reason.append("Reference high-affinity epitope.")
        if res['Epitope'] == 'E3': reason.append("Very poor P9 anchor (K).")
        if res['P2_Score'] < 10 and res['Epitope'] != 'E3': reason.append("Sub-optimal P2 anchor (F).")
        if res['Core_Penalty'] != 0 and 'P' in res['Sequence']: reason.append("Proline in core sequence disrupts binding.")
        elif res['Core_Penalty'] != 0: reason.append("Substitutions in core sequence.")
        
        print(f"{i+1:<5} {res['Epitope']:<8} {res['Sequence']:<12} {res['Total_Score']:<6.1f} {' '.join(reason)}")
    print("=" * 60)

    # Print the final ranking equation
    final_ranking_str = " > ".join([f"{res['Epitope']} ({res['Sequence']})" for res in ranked_epitopes])
    print("\nFinal Predicted Binding Affinity Ranking:")
    print(final_ranking_str)


# Run the analysis and print the results
rank_epitopes()