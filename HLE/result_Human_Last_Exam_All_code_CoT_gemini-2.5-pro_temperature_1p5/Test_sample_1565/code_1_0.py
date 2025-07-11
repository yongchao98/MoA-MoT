import sys

def calculate_fat_score():
    """
    Calculates the FAT score for user B based on the provided data and rules.
    """
    # --- Data Definition ---
    # P1 (Trust Incoming): 7 total ratings (4 trust, 2 trust, 1 distrust)
    p1_total = 7

    # P2 (Trust Incoming): 6 total ratings (3 trust, 1 trust, 2 distrust)
    p2_total = 6

    # P3 (Trust Incoming): 4 total ratings (2 trust, 1 distrust, 1 distrust)
    p3_total = 4

    # N1 (Distrust Incoming): 6 total ratings (3 trust, 2 distrust, 1 distrust)
    n1_total = 6
    n1_trust_ratings = 3
    n1_distrust_ratings = 3

    # N2 (Distrust Incoming): 4 total ratings (1 trust, 3 distrust)
    n2_total = 4
    n2_trust_ratings = 1
    n2_distrust_ratings = 3


    # --- Calculation ---
    # Rule 1: Positive edge contributes 1/(total_relationships + 1)
    p1_contrib = 1 / (p1_total + 1)
    p2_contrib = 1 / (p2_total + 1)
    p3_contrib = 1 / (p3_total + 1)

    # Rule 2: Negative edge with mixed ratings: -1 /(total relationships+1) x (trust_ratings/total)
    n1_contrib_base = -1 / (n1_total + 1) * (n1_trust_ratings / n1_total)
    n2_contrib_base = -1 / (n2_total + 1) * (n2_trust_ratings / n2_total)

    # Rule 3: Users with more distrust than trust ratings get 1.5x negative weight
    # Check N1
    n1_final_contrib = n1_contrib_base
    n1_multiplier_note = ""
    if n1_distrust_ratings > n1_trust_ratings:
        n1_final_contrib *= 1.5
        n1_multiplier_note = " (1.5x multiplier applied)"
    else:
        n1_multiplier_note = " (no multiplier)"

    # Check N2
    n2_final_contrib = n2_contrib_base
    n2_multiplier_note = ""
    if n2_distrust_ratings > n2_trust_ratings:
        n2_final_contrib *= 1.5
        n2_multiplier_note = " (1.5x multiplier applied)"
    else:
        n2_multiplier_note = " (no multiplier)"
        
    # Final Summation
    total_score = p1_contrib + p2_contrib + p3_contrib + n1_final_contrib + n2_final_contrib

    # --- Output ---
    print("Calculating B's FAT (Fringe of Absolute Trust) Score:")
    print("-" * 50)
    
    print("1. Contributions from Positive Edges (Trust Incoming):")
    print(f"   P1: 1 / ({p1_total} + 1) = {p1_contrib:.4f}")
    print(f"   P2: 1 / ({p2_total} + 1) = {p2_contrib:.4f}")
    print(f"   P3: 1 / ({p3_total} + 1) = {p3_contrib:.4f}")
    print("\n2. Contributions from Negative Edges (Distrust Incoming):")
    print(f"   N1: -1 / ({n1_total} + 1) * ({n1_trust_ratings}/{n1_total}) = {n1_final_contrib:.4f}{n1_multiplier_note}")
    print(f"   N2: [-1 / ({n2_total} + 1) * ({n2_trust_ratings}/{n2_total})] * 1.5 = {n2_final_contrib:.4f}{n2_multiplier_note}")
    
    print("\n3. Final Score Calculation (Sum of all contributions):")
    # Print the final equation with all numbers
    # We use round to avoid floating point representation issues in the final sum display
    print(f"   Equation: {p1_contrib:.4f} + {p2_contrib:.4f} + {p3_contrib:.4f} + ({n1_final_contrib:.4f}) + ({n2_final_contrib:.4f})")
    print(f"   Result: {total_score:.4f}")
    print("-" * 50)

    print(f"\nThe calculated score is {total_score:.4f}, which is closest to the option 0.35.")

calculate_fat_score()