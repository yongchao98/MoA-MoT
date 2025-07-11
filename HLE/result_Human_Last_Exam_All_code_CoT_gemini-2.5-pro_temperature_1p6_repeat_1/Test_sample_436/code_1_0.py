import pandas as pd

def rank_epitopes():
    """
    Ranks peptide epitopes based on their predicted binding affinity to H-2Kd.
    
    A scoring system is devised based on known anchor residue preferences for the H-2Kd allele.
    - P2 Anchor: Y is optimal.
    - P9 Anchor: V, L, I are optimal. K is strongly disfavored.
    - P3 Secondary Anchor: Q is preferred.
    """
    
    epitopes = {
        'E1': 'TYQRTRALV',
        'E2': 'TFQRTRALV',
        'E3': 'TFQRTRALK',
        'E4': 'TYQMFALV',
        'E5': 'TYIPMRALV'
    }

    # Scoring based on H-2Kd anchor preferences
    # Higher score means better binding
    p2_scores = {'Y': 10, 'F': 2, 'T': 6, 'I': 6}
    p3_scores = {'Q': 5, 'I': 3, 'R': 2}
    p9_scores = {'V': 9, 'L': 10, 'I': 10, 'K': -5}

    results = []

    print("Calculating affinity scores for H-2Kd binding:")
    for name, seq in epitopes.items():
        p2 = seq[1]
        p3 = seq[2]
        p9 = seq[8]

        score_p2 = p2_scores.get(p2, 0)
        score_p3 = p3_scores.get(p3, 0)
        score_p9 = p9_scores.get(p9, 0)

        total_score = score_p2 + score_p3 + score_p9
        
        # We must output each number in the final equation.
        print(f"{name} ({seq}): P2='{p2}'({score_p2}) + P3='{p3}'({score_p3}) + P9='{p9}'({score_p9}) = {total_score}")

        results.append({'Epitope': name, 'Sequence': seq, 'Score': total_score})

    # Create a DataFrame and sort by score
    df = pd.DataFrame(results)
    ranked_df = df.sort_values(by='Score', ascending=False)

    print("\nEpitopes ranked by predicted binding amount (highest to lowest):")
    ranked_list = ranked_df['Epitope'].tolist()
    print(", ".join(ranked_list))
    
    # Based on the ranking, determine the correct answer choice
    final_ranking_str = ", ".join(ranked_list)
    answer_choices = {
        "A": "E1, E4, E5, E2, E3",
        "B": "E4, E5, E1, E2, E3",
        "C": "E5, E4, E1, E2, E3",
        "D": "E1, E5, E4, E2, E3",
        "E": "E5, E4, E1, E3, E2"
    }

    for choice, ranking in answer_choices.items():
        if ranking == final_ranking_str:
            print(f"\nThis corresponds to Answer Choice {choice}.")
            # The final answer format as requested by the prompt.
            print(f"\n<<<{choice}>>>")
            return

rank_epitopes()
