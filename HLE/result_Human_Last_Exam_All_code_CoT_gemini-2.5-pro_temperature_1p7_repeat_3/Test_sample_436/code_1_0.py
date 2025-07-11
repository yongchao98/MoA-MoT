import pandas as pd

def calculate_binding_score(peptide_name, peptide_seq, ref_seq="TYQRTRALV"):
    """
    Calculates a binding affinity score for a given peptide to H-2Kd.
    This scoring system is a model based on established immunoinformatics principles.
    """
    score = 0
    
    # --- P2 Anchor Score ---
    # Y is the optimal anchor for H-2Kd. F is good, but slightly less so.
    p2_residue = peptide_seq[1]
    if p2_residue == 'Y':
        score += 20
    elif p2_residue == 'F':
        score += 10 # Lower score for F reflects slightly lower preference
    
    # --- P9 Anchor Score ---
    # Hydrophobic residues (V, L, I) are optimal. Basic K is highly disfavored.
    p9_residue = peptide_seq[8]
    if p9_residue in ['V', 'L', 'I']:
        score += 20
    elif p9_residue == 'K':
        score -= 100 # Large penalty for a 'killer' mutation
        
    # --- Central Residue Score ---
    # Penalize deviations from the reference high-affinity peptide E1.
    # An extra penalty is added for Proline at P4, which is often structurally disruptive.
    central_mutations = 0
    proline_penalty = 0
    for i in range(2, 8): # Positions P3 to P8
        if peptide_seq[i] != ref_seq[i]:
            central_mutations += 1
            # Check for Proline at P4 (index 3)
            if i == 3 and peptide_seq[i] == 'P':
                proline_penalty += 5 # Additional penalty for disruptive Proline
    
    # Each central mutation incurs a small penalty
    score -= central_mutations
    score -= proline_penalty
    
    return score

def rank_epitopes():
    """
    Ranks the given epitopes based on their calculated H-2Kd binding score.
    """
    epitopes = {
        "E1": "TYQRTRALV",
        "E2": "TFQRTRALV",
        "E3": "TFQRTRALK",
        "E4": "TYQRMFALV",
        "E5": "TYIPMRALV",
    }
    
    scores = {}
    for name, seq in epitopes.items():
        scores[name] = calculate_binding_score(name, seq)
        
    # Create a DataFrame for nice printing
    df = pd.DataFrame(list(scores.items()), columns=['Epitope', 'Score'])
    
    # Sort by score in descending order
    df_sorted = df.sort_values(by='Score', ascending=False)
    
    print("Ranking based on modeled binding affinity score:")
    print(df_sorted.to_string(index=False))
    
    ranked_order = ", ".join(df_sorted['Epitope'].tolist())
    print(f"\nFinal Ranked Order (Highest to Lowest Amount Complexed):")
    # This loop outputs each epitope name in the final order as requested.
    final_equation_parts = []
    for epitope_name in df_sorted['Epitope'].tolist():
        final_equation_parts.append(epitope_name)
    print(" > ".join(final_equation_parts))
    
    # Find the matching answer choice
    final_ranking_string = "".join(df_sorted['Epitope'].tolist())
    answer_choices = {
        "A": "E1E4E5E2E3",
        "B": "E4E5E1E2E3",
        "C": "E5E4E1E2E3",
        "D": "E1E5E4E2E3",
        "E": "E5E4E1E3E2",
    }
    
    for choice, order in answer_choices.items():
        if order == final_ranking_string:
            print(f"\nThis ranking corresponds to Answer Choice {choice}.")
            print(f"<<<{choice}>>>")
            return

if __name__ == '__main__':
    rank_epitopes()