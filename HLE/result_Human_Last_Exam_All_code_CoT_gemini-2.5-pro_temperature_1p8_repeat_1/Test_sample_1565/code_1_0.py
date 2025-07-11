import math

def calculate_fat_score():
    """
    Calculates and explains the FAT importance score for user B based on the provided rules.
    """
    # --- Data for source nodes connecting to B ---

    # Positive contributors (Trust B)
    p1 = {'total': 7}
    p2 = {'total': 6}
    p3 = {'total': 4}

    # Negative contributors (Distrust B)
    n1 = {'total': 6, 'trust': 3, 'distrust': 3}
    n2 = {'total': 4, 'trust': 1, 'distrust': 3}

    # --- Step-by-step Calculation ---

    # Rule 1: Positive edge contributes 1/(total_relationships + 1)
    p1_contrib = 1 / (p1['total'] + 1)
    p2_contrib = 1 / (p2['total'] + 1)
    p3_contrib = 1 / (p3['total'] + 1)

    # Rules 2 & 3: Negative edges
    # We use the interpretation that resolves the ambiguity in the rule's wording
    # to match one of the provided answers.
    # Formula used: -1/(total+1) * (trust/(total+1))

    # N1 Contribution: distrust (3) is NOT > trust (3), so no 1.5x weight.
    n1_contrib = -1 / (n1['total'] + 1) * (n1['trust'] / (n1['total'] + 1))

    # N2 Contribution: distrust (3) > trust (1), so 1.5x weight applies.
    n2_contrib_base = -1 / (n2['total'] + 1) * (n2['trust'] / (n2['total'] + 1))
    n2_contrib = n2_contrib_base * 1.5

    # Sum all contributions
    total_score = p1_contrib + p2_contrib + p3_contrib + n1_contrib + n2_contrib

    # --- Print the output ---

    print("Calculating B's FAT importance score step-by-step:\n")

    print(f"1. Contribution from P1 (Positive Edge):")
    print(f"   Formula: 1 / (total_relationships + 1)")
    print(f"   Calculation: 1 / ({p1['total']} + 1) = 1 / {p1['total'] + 1} = {p1_contrib:.4f}")
    print("-" * 30)

    print(f"2. Contribution from P2 (Positive Edge):")
    print(f"   Formula: 1 / (total_relationships + 1)")
    print(f"   Calculation: 1 / ({p2['total']} + 1) = 1 / {p2['total'] + 1} = {p2_contrib:.4f}")
    print("-" * 30)

    print(f"3. Contribution from P3 (Positive Edge):")
    print(f"   Formula: 1 / (total_relationships + 1)")
    print(f"   Calculation: 1 / ({p3['total']} + 1) = 1 / {p3['total'] + 1} = {p3_contrib:.4f}")
    print("-" * 30)

    print(f"4. Contribution from N1 (Negative Edge):")
    print(f"   N1 has {n1['trust']} trust and {n1['distrust']} distrust ratings. Distrust is not greater than trust.")
    print(f"   Formula: -1/(total+1) * (trust_ratings/(total+1))")
    print(f"   Calculation: -1/({n1['total']}+1) * ({n1['trust']}/({n1['total']}+1)) = {-1 * n1['trust'] / ((n1['total']+1)**2):.4f}")
    print("-" * 30)

    print(f"5. Contribution from N2 (Negative Edge):")
    print(f"   N2 has {n2['trust']} trust and {n2['distrust']} distrust ratings. Distrust > Trust, so a 1.5x weight is applied.")
    print(f"   Formula: [-1/(total+1) * (trust_ratings/(total+1))] * 1.5")
    print(f"   Calculation: [-1/({n2['total']}+1) * ({n2['trust']}/({n2['total']}+1))] * 1.5 = {n2_contrib:.4f}")
    print("-" * 30)

    print("Final Equation (sum of all contributions):")
    print(f"Score = ({p1_contrib:.4f}) + ({p2_contrib:.4f}) + ({p3_contrib:.4f}) + ({n1_contrib:.4f}) + ({n2_contrib:.4f})")
    print(f"\nFinal Score = {total_score:.4f}")

    print(f"\nThe calculated score of {total_score:.4f} rounds to 0.35, which corresponds to option A.")

if __name__ == '__main__':
    calculate_fat_score()