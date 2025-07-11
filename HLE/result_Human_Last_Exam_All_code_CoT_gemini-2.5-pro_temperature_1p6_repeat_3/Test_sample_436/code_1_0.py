import pandas as pd

def calculate_binding_score(sequence, reference_sequence="TYQRTRALV"):
    """
    Calculates a simplified binding score for a peptide to H2-Kd based on deviations
    from a reference high-affinity sequence.
    """
    score = 100
    p2_anchor = sequence[1]
    p9_anchor = sequence[8]
    
    # 1. Penalize for bad P9 anchor
    if p9_anchor not in ['V', 'L']:
        score -= 60  # Heavy penalty for bad P9 anchor
        
    # 2. Penalize for non-reference P2 anchor
    if p2_anchor != reference_sequence[1]:
        # F is a good anchor, but Y is reference. Small penalty.
        if p2_anchor == 'F':
            score -= 10 
        else:
            score -= 20 # Heavier penalty for other non-reference anchors
            
    # 3. Penalize for internal changes, with extra penalty for proline
    for i in range(2, 8): # Check internal residues P3-P8
        if sequence[i] != reference_sequence[i]:
            if sequence[i] == 'P':
                score -= 5 # Proline is structurally disruptive
            else:
                score -= 2 # Penalty for other internal changes
                
    return score

def rank_epitopes():
    """
    Ranks epitopes based on a calculated binding score and prints the result.
    """
    epitopes = {
        'E1': 'TYQRTRALV',
        'E2': 'TFQRTRALV',
        'E3': 'TFQRTRALK',
        'E4': 'TYQRMFALV',
        'E5': 'TYIPMRALV'
    }
    
    scores = {}
    for name, seq in epitopes.items():
        scores[name] = calculate_binding_score(seq)
        
    # Create a DataFrame for nice printing
    df = pd.DataFrame(list(scores.items()), columns=['Epitope', 'Score'])
    df_sorted = df.sort_values(by='Score', ascending=False).reset_index(drop=True)
    
    print("--- Epitope Binding Score Analysis ---")
    print("The ranking is based on a scoring system where deviations from the ideal reference epitope (E1) are penalized.")
    print("Major penalties are applied for non-optimal anchor residues (P2 and P9).\n")
    print(df_sorted.to_string())
    print("\n--- Final Ranking ---")
    
    ranked_list = df_sorted['Epitope'].tolist()
    print("The predicted rank from highest to lowest binding is:")
    print(', '.join(ranked_list))

# Run the ranking
rank_epitopes()

<<<A>>>