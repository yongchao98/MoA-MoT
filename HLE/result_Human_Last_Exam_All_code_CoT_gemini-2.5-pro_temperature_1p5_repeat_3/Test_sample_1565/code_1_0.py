def calculate_fat_score():
    """
    Calculates the Fringe of Absolute Trust (FAT) score for user B based on the provided rules.
    """

    # --- Data for incoming connections to B ---
    # P1: Trust, 7 total, 4 trust, 3 distrust
    # P2: Trust, 6 total, 4 trust, 2 distrust
    # P3: Trust, 4 total, 2 trust, 2 distrust
    # N1: Distrust, 6 total, 3 trust, 3 distrust
    # N2: Distrust, 4 total, 1 trust, 3 distrust

    # --- Rule 1: Positive Edges ---
    p1_total = 7
    p1_score = 1 / (p1_total + 1)

    p2_total = 6
    p2_score = 1 / (p2_total + 1)

    p3_total = 4
    p3_score = 1 / (p3_total + 1)

    # --- Rules 2 & 3: Negative Edges ---

    # N1: Negative edge with mixed ratings
    n1_total = 6
    n1_trust = 3
    n1_distrust = 3
    # Based on testing, the 'total' in 'trust_ratings/total' likely means 'total_relationships + 1'
    n1_score = -1 / (n1_total + 1) * (n1_trust / (n1_total + 1))
    # Rule 3 check: distrust > trust? (3 > 3 is false) No 1.5x multiplier.

    # N2: Negative edge with mixed ratings
    n2_total = 4
    n2_trust = 1
    n2_distrust = 3
    n2_base_score = -1 / (n2_total + 1) * (n2_trust / (n2_total + 1))
    # Rule 3 check: distrust > trust? (3 > 1 is true) Apply 1.5x multiplier.
    n2_score = n2_base_score * 1.5

    # --- Total Score Calculation ---
    total_score = p1_score + p2_score + p3_score + n1_score + n2_score

    # --- Output ---
    print("Calculating B's FAT Score:\n")
    print(f"P1 (Trust): 1 / ({p1_total} + 1) = {p1_score:.4f}")
    print(f"P2 (Trust): 1 / ({p2_total} + 1) = {p2_score:.4f}")
    print(f"P3 (Trust): 1 / ({p3_total} + 1) = {p3_score:.4f}")
    print(f"N1 (Distrust): -1 / ({n1_total} + 1) * ({n1_trust} / ({n1_total} + 1)) = {n1_score:.4f}")
    print(f"N2 (Distrust): [-1 / ({n2_total} + 1) * ({n2_trust} / ({n2_total} + 1))] * 1.5 = {n2_score:.4f}")
    print("-" * 30)
    print(f"Final Equation: {p1_score:.4f} + {p2_score:.4f} + {p3_score:.4f} + ({n1_score:.4f}) + ({n2_score:.4f})")
    print(f"Total Score = {total_score:.4f}")

calculate_fat_score()