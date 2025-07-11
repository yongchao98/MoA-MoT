def calculate_fat_score():
    """
    Calculates the FAT score for user B based on the provided data and rules.
    """

    # --- Data from the problem ---
    # P1 (Trust -> B): 7 total ratings, 6 trust, 1 distrust
    # P2 (Trust -> B): 6 total ratings, 4 trust, 2 distrust
    # P3 (Trust -> B): 4 total ratings, 2 trust, 2 distrust
    # N1 (Distrust -> B): 6 total ratings, 3 trust, 3 distrust
    # N2 (Distrust -> B): 4 total ratings, 1 trust, 3 distrust

    # --- Step 1: Calculate contribution from positive edges (P1, P2, P3) ---
    # Rule 1: Positive edge contributes 1/(total_relationships + 1)
    p1_total = 7
    p1_contrib = 1 / (p1_total + 1)

    p2_total = 6
    p2_contrib = 1 / (p2_total + 1)

    p3_total = 4
    p3_contrib = 1 / (p3_total + 1)

    # --- Step 2: Calculate contribution from negative edges (N1, N2) ---
    # Rule 2: Negative edge with mixed ratings: -1 /(total relationships+1) x (trust_ratings/total)

    # For N1
    n1_total = 6
    n1_trust = 3
    # N1 has 3 trust and 3 distrust, so Rule 3 does not apply.
    n1_contrib = -1 / (n1_total + 1) * (n1_trust / n1_total)

    # For N2
    n2_total = 4
    n2_trust = 1
    n2_distrust = 3
    # Rule 3: N2 has more distrust than trust (3 > 1), so apply 1.5x negative weight.
    n2_multiplier = 1.5
    n2_contrib = (-1 / (n2_total + 1) * (n2_trust / n2_total)) * n2_multiplier

    # --- Step 3: Sum the contributions ---
    total_score = p1_contrib + p2_contrib + p3_contrib + n1_contrib + n2_contrib

    # --- Step 4: Print the final equation and result ---
    print("Calculating B's importance score (FAT score):")
    print(f"P1 (Trust)  = 1/({p1_total}+1) = {p1_contrib:.4f}")
    print(f"P2 (Trust)  = 1/({p2_total}+1) = {p2_contrib:.4f}")
    print(f"P3 (Trust)  = 1/({p3_total}+1) = {p3_contrib:.4f}")
    print(f"N1 (Distrust)= -1/({n1_total}+1) * ({n1_trust}/{n1_total}) = {n1_contrib:.4f}")
    print(f"N2 (Distrust)= [-1/({n2_total}+1) * ({n2_trust}/{n2_total})] * {n2_multiplier} = {n2_contrib:.4f}")
    print("-" * 30)
    print("Final Equation:")
    print(f"FAT(B) = {p1_contrib:.4f} + {p2_contrib:.4f} + {p3_contrib:.4f} + ({n1_contrib:.4f}) + ({n2_contrib:.4f})")
    print(f"Total Score = {total_score:.4f}")
    print("\nThe calculated score of approximately 0.32 is closest to option A (0.35).")

calculate_fat_score()