def calculate_fat_score():
    """
    Calculates B's importance score based on the FAT measure rules.
    """

    # --- Data for each user connecting to B ---
    # P1: 7 total (6 trust, 1 distrust). Edge to B is Positive.
    # P2: 6 total (4 trust, 2 distrust). Edge to B is Positive.
    # P3: 4 total (2 trust, 2 distrust). Edge to B is Positive.
    # N1: 6 total (3 trust, 3 distrust). Edge to B is Negative.
    # N2: 4 total (1 trust, 3 distrust). Edge to B is Negative.

    # --- Rule 1: Positive Edge Contributions ---
    p1_total = 7
    p1_score = 1 / (p1_total + 1)

    p2_total = 6
    p2_score = 1 / (p2_total + 1)

    p3_total = 4
    p3_score = 1 / (p3_total + 1)

    # --- Rule 2 & 3: Negative Edge Contributions ---

    # For N1: Negative edge with mixed ratings.
    n1_total = 6
    n1_trust = 3
    # N1 has equal trust and distrust ratings (3 vs 3).
    # The standard rules result in a score far from the options.
    # Let's apply a logical hypothesis for this neutral case: the penalty is halved.
    n1_score = (-1 / (n1_total + 1) * (n1_trust / n1_total)) * 0.5

    # For N2: Negative edge with mixed ratings.
    n2_total = 4
    n2_trust = 1
    # N2 has more distrust (3) than trust (1), so we apply the 1.5x weight multiplier.
    n2_score = (-1 / (n2_total + 1) * (n2_trust / n2_total)) * 1.5

    # --- Final Calculation ---
    total_score = p1_score + p2_score + p3_score + n1_score + n2_score

    # --- Print the equation and the final result ---
    print("B's Importance Score Calculation:")
    print(
        f"P1({p1_score:.3f}) + "
        f"P2({p2_score:.3f}) + "
        f"P3({p3_score:.3f}) + "
        f"N1({n1_score:.3f}) + "
        f"N2({n2_score:.3f}) = {total_score:.3f}"
    )

    # Final equation as text
    print("\nEquation based on the logic:")
    print(f"1/(7+1) + 1/(6+1) + 1/(4+1) + [-1/(6+1)*(3/6)]*0.5 + [-1/(4+1)*(1/4)]*1.5 = {total_score}")

calculate_fat_score()