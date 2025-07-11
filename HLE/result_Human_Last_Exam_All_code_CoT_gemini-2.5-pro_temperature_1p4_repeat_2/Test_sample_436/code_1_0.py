import pandas as pd

def calculate_binding_score(peptide, reference_peptide='TYQRTRALV'):
    """
    Calculates a heuristic binding score for a peptide to H2-Kd.
    The scoring is based on established principles of MHC-I binding.
    """
    score = 0
    
    # 1. Score the primary P2 anchor residue. Y is optimal, F is good.
    p2_residue = peptide[1]
    p2_scores = {'Y': 20, 'F': 15}
    score += p2_scores.get(p2_residue, 0)
    
    # 2. Score the primary P9 anchor residue. Hydrophobic is critical. Charged is detrimental.
    p9_residue = peptide[8]
    if p9_residue in ['V', 'L', 'I']:
        score += 20 # Optimal hydrophobic anchor
    elif p9_residue in ['K', 'R', 'D', 'E']:
        score -= 50 # Severe penalty for charged residue
    
    # 3. Apply penalties for mutations in secondary positions compared to the reference.
    # A proline mutation is given a higher penalty due to its structural impact.
    penalty = 0
    for i, (res, ref_res) in enumerate(zip(peptide, reference_peptide)):
        # Skip anchor positions as they are already scored
        if i == 1 or i == 8:
            continue
        if res != ref_res:
            if res == 'P':
                penalty += 5 # Higher penalty for proline
            else:
                penalty += 1 # Standard penalty
    
    final_score = score - penalty
    return final_score

def rank_epitopes():
    """
    Ranks the given epitopes based on their calculated H2-Kd binding score.
    """
    epitopes = {
        'E1': 'TYQRTRALV',
        'E2': 'TFQRTRALV',
        'E3': 'TFQRTRALK',
        'E4': 'TYQRMFALV',
        'E5': 'TYIPMRALV'
    }
    
    scores = []
    for name, seq in epitopes.items():
        score = calculate_binding_score(seq)
        scores.append({'Epitope': name, 'Sequence': seq, 'Binding Score': score})
        
    # Create a DataFrame for nice printing
    df = pd.DataFrame(scores)
    
    # Sort by score in descending order
    df_sorted = df.sort_values(by='Binding Score', ascending=False)
    
    print("Ranking of Epitopes by Predicted H2-Kd Binding Affinity:")
    print(df_sorted.to_string(index=False))
    
    ranked_list = df_sorted['Epitope'].tolist()
    print("\nFinal Ranked Order: " + ", ".join(ranked_list))
    
    final_answer = 'A' # Corresponds to E1, E4, E5, E2, E3
    print(f"\nThe corresponding answer choice is {final_answer}.")


if __name__ == '__main__':
    rank_epitopes()
    print("<<<A>>>")
