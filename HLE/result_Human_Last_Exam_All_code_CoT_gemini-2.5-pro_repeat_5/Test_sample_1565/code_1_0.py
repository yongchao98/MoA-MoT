import sys

def calculate_fat_score():
    """
    Calculates the FAT score for user B based on the provided rules and data.
    """
    # Step 1: Define the data for each user connecting to B.
    # Data is parsed from the problem description:
    # P1 (7 total): 4 trust-senior, 2 trust-peer (6 trust), 1 distrust-peer (1 distrust)
    # P2 (6 total): 3 trust-junior, 1 trust-senior (4 trust), 2 distrust-peer (2 distrust)
    # P3 (4 total): 2 trust-peer (2 trust), 1 distrust-senior, 1 distrust-junior (2 distrust)
    # N1 (6 total): 3 trust-senior (3 trust), 2 distrust-peer, 1 distrust-junior (3 distrust)
    # N2 (4 total): 1 trust-peer (1 trust), 3 distrust-junior (3 distrust)
    nodes = {
        'P1': {'total': 7, 'trust': 6, 'distrust': 1},
        'P2': {'total': 6, 'trust': 4, 'distrust': 2},
        'P3': {'total': 4, 'trust': 2, 'distrust': 2},
        'N1': {'total': 6, 'trust': 3, 'distrust': 3},
        'N2': {'total': 4, 'trust': 1, 'distrust': 3}
    }

    total_score = 0.0
    
    print("Calculating B's FAT score step-by-step:\n")

    # Step 2: Calculate contribution from positive edges (P1, P2, P3)
    # P1
    p1_total = nodes['P1']['total']
    p1_contribution = 1 / (p1_total + 1)
    total_score += p1_contribution
    print(f"P1 (Trust Edge):")
    print(f"Formula: 1 / (total_relationships + 1)")
    print(f"Calculation: 1 / ({p1_total} + 1) = {p1_contribution:.4f}\n")

    # P2
    p2_total = nodes['P2']['total']
    p2_contribution = 1 / (p2_total + 1)
    total_score += p2_contribution
    print(f"P2 (Trust Edge):")
    print(f"Formula: 1 / (total_relationships + 1)")
    print(f"Calculation: 1 / ({p2_total} + 1) = {p2_contribution:.4f}\n")

    # P3
    p3_total = nodes['P3']['total']
    p3_contribution = 1 / (p3_total + 1)
    total_score += p3_contribution
    print(f"P3 (Trust Edge):")
    print(f"Formula: 1 / (total_relationships + 1)")
    print(f"Calculation: 1 / ({p3_total} + 1) = {p3_contribution:.4f}\n")

    # Step 3, 4, 5: Calculate contribution from negative edges (N1, N2)
    # N1
    n1_total = nodes['N1']['total']
    n1_trust = nodes['N1']['trust']
    n1_distrust = nodes['N1']['distrust']
    n1_contribution = -1 / (n1_total + 1) * (n1_trust / n1_total)
    # Rule 3 check: n1_distrust (3) is not > n1_trust (3), so no penalty.
    total_score += n1_contribution
    print(f"N1 (Distrust Edge):")
    print(f"Formula: -1 / (total_relationships + 1) * (trust_ratings / total)")
    print(f"Calculation: -1 / ({n1_total} + 1) * ({n1_trust} / {n1_total}) = {n1_contribution:.4f}\n")

    # N2
    n2_total = nodes['N2']['total']
    n2_trust = nodes['N2']['trust']
    n2_distrust = nodes['N2']['distrust']
    n2_contribution = -1 / (n2_total + 1) * (n2_trust / n2_total)
    # Rule 3 check: n2_distrust (3) is > n2_trust (1), so 1.5x penalty applies.
    n2_contribution *= 1.5
    total_score += n2_contribution
    print(f"N2 (Distrust Edge):")
    print(f"Note: User N2 has more distrust than trust ratings (3 > 1), so a 1.5x weight is applied.")
    print(f"Formula: [-1 / (total_relationships + 1) * (trust_ratings / total)] * 1.5")
    print(f"Calculation: [-1 / ({n2_total} + 1) * ({n2_trust} / {n2_total})] * 1.5 = {n2_contribution:.4f}\n")

    # Step 6: Final score
    print("-----------------------------------------")
    print("Final Equation:")
    print(f"Score = (1/{p1_total+1}) + (1/{p2_total+1}) + (1/{p3_total+1}) + (-1/{n1_total+1} * {n1_trust}/{n1_total}) + (-1/{n2_total+1} * {n2_trust}/{n2_total} * 1.5)")
    print("\nTotal Score:")
    print(f"Score = {p1_contribution:.4f} + {p2_contribution:.4f} + {p3_contribution:.4f} + ({n1_contribution:.4f}) + ({n2_contribution:.4f})")
    print(f"Score = {total_score:.4f}")
    
    options = {'A': 0.35, 'B': 0.45, 'C': 0.55, 'D': 0.64}
    closest_option = min(options, key=lambda k: abs(options[k] - total_score))
    print(f"\nThe calculated score {total_score:.4f} is closest to option {closest_option}) {options[closest_option]}.")
    # This is a special instruction for the AI model to output the final answer in a specific format.
    # It is not intended to be user-facing.
    sys.stdout.flush()
    print(f'<<<{closest_option}>>>')


calculate_fat_score()