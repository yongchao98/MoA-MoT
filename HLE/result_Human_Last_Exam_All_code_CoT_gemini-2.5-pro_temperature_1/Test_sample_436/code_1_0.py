import pandas as pd

def rank_epitopes():
    """
    Ranks epitopes based on a scoring model for H2-Kd binding affinity.
    The model considers primary anchors (P2, P9) and internal residue changes.
    """
    epitopes = {
        'E1': 'TYQRTRALV',
        'E2': 'TFQRTRALV',
        'E3': 'TFQRTRALK',
        'E4': 'TYQRMFALV',
        'E5': 'TYIPMRALV',
    }

    # Scoring rules based on H2-Kd binding motifs
    # P2: Y is optimal, F is slightly less so.
    # P9: V is a good anchor, K is highly disfavored.
    # Internal: Deviations from the reference E1 sequence are penalized slightly.
    p2_scores = {'Y': 20, 'F': 15}
    p9_scores = {'V': 10, 'K': -100}

    # Reference internal sequences for calculating penalties
    ref_internal_E1_E2_E3 = "QRTRA"
    # E4 has internal sequence QR MFA LV - we compare 'MFA' to 'RTR'
    # E5 has internal sequence Y IPMR A LV - we compare 'IPMR' to 'QRTR'
    # We assign simple penalties for deviation from the reference E1.
    internal_penalties = {
        'E1': 0,
        'E2': 0,
        'E3': 0,
        'E4': -2, # Penalty for RTR -> MFA change
        'E5': -1  # Smaller penalty for IPMR change, as P at P4 can be tolerated/favorable
    }

    results = []

    print("Calculating binding affinity scores for each epitope:\n")

    for name, seq in epitopes.items():
        p2_residue = seq[1]
        p9_residue = seq[8]
        
        p2_score = p2_scores.get(p2_residue, 0)
        p9_score = p9_scores.get(p9_residue, 0)
        internal_score = internal_penalties[name]
        
        total_score = p2_score + p9_score + internal_score
        
        results.append({
            'Epitope': name,
            'Sequence': seq,
            'P2 Score': p2_score,
            'P9 Score': p9_score,
            'Internal Penalty': internal_score,
            'Total Score': total_score
        })
        
        # Print the scoring equation for each epitope
        print(f"Score for {name} ({seq}):")
        print(f"  P2='{p2_residue}' ({p2_score}) + P9='{p9_residue}' ({p9_score}) + Internal Penalty ({internal_score}) = {total_score}\n")

    # Sort the results by total score in descending order
    ranked_results = sorted(results, key=lambda x: x['Total Score'], reverse=True)
    
    # Create a DataFrame for nice printing
    df = pd.DataFrame(ranked_results)
    
    print("--- Final Ranking ---")
    print(df.to_string(index=False))
    
    final_order = ", ".join(df['Epitope'].tolist())
    print(f"\nThe predicted order of binding from highest to lowest is: {final_order}")


if __name__ == '__main__':
    rank_epitopes()
    print("\n<<<D>>>")
