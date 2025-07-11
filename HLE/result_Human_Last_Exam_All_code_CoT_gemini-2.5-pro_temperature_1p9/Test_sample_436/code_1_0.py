def rank_epitopes():
    """
    Ranks epitopes based on a scoring model for H2-Kd binding affinity.
    The model is based on known anchor residue preferences for the H2-Kd allele.
    """
    epitopes = {
        "E1": "TYQRTRALV",
        "E2": "TFQRTRALV",
        "E3": "TFQRTRALK",
        "E4": "TYQRMFALV",
        "E5": "TYIPMRALV",
    }

    # Simplified scoring rules based on H2-Kd binding motif
    # Primary anchors (P2, P9) are weighted heavily.
    # Secondary anchor (P5) has a smaller contribution.
    scoring_rules = {
        'P2': {'Y': 100, 'F': 30},
        'P5': {'R': 5, 'M': 2, 'P': -5},
        'P9': {'V': 100, 'K': -100}
    }

    epitope_scores = {}
    
    print("Calculating binding scores for each epitope:\n")

    for name, seq in epitopes.items():
        # Get amino acids at positions of interest (P2, P5, P9)
        p2_aa = seq[1]
        p5_aa = seq[4]
        p9_aa = seq[8]

        # Get scores from rules, defaulting to 0 if not found
        p2_score = scoring_rules['P2'].get(p2_aa, 0)
        p5_score = scoring_rules['P5'].get(p5_aa, 0)
        p9_score = scoring_rules['P9'].get(p9_aa, 0)

        total_score = p2_score + p5_score + p9_score
        epitope_scores[name] = total_score
        
        print(f"{name} ({seq}):")
        print(f"  Score = P2({p2_aa}) + P5({p5_aa}) + P9({p9_aa})")
        print(f"  Score = {p2_score} + {p5_score} + {p9_score} = {total_score}\n")


    # Sort the epitopes by score in descending order
    # The key for sorting is the second element (the score) of each item in the dictionary
    sorted_epitopes = sorted(epitope_scores.items(), key=lambda item: item[1], reverse=True)

    # Print the final ranked list
    print("="*40)
    print("Final Ranking (from highest to lowest binding affinity):")
    ranked_list = [name for name, score in sorted_epitopes]
    print(", ".join(ranked_list))
    print("="*40)

if __name__ == '__main__':
    rank_epitopes()
