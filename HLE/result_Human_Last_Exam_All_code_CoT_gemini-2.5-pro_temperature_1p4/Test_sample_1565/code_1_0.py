def calculate_fat_score():
    """
    Calculates the FAT score for user B based on the provided data and rules.
    """
    # Data for each user connected to B
    # P = Positive edge, N = Negative edge
    # format: (total_ratings, trust_ratings)
    user_data = {
        'P1': {'type': 'positive', 'total': 7, 'trust': 4},
        'P2': {'type': 'positive', 'total': 6, 'trust': 4}, # 3 trust-junior + 1 trust-senior
        'P3': {'type': 'positive', 'total': 4, 'trust': 2},
        'N1': {'type': 'negative', 'total': 6, 'trust': 3, 'distrust': 3}, # 2 distrust-peer + 1 distrust-junior
        'N2': {'type': 'negative', 'total': 4, 'trust': 1, 'distrust': 3},
    }

    total_score = 0
    calculation_parts = []

    # Calculate P1 contribution
    p1_total = user_data['P1']['total']
    p1_score = 1 / (p1_total + 1)
    total_score += p1_score
    calculation_parts.append(f"(1 / ({p1_total} + 1))")

    # Calculate P2 contribution
    p2_total = user_data['P2']['total']
    p2_score = 1 / (p2_total + 1)
    total_score += p2_score
    calculation_parts.append(f"(1 / ({p2_total} + 1))")
    
    # Calculate P3 contribution
    p3_total = user_data['P3']['total']
    p3_score = 1 / (p3_total + 1)
    total_score += p3_score
    calculation_parts.append(f"(1 / ({p3_total} + 1))")

    # Calculate N1 contribution (Negative, mixed, trust not < distrust)
    n1_total = user_data['N1']['total']
    n1_trust = user_data['N1']['trust']
    n1_score = -1 / (n1_total + 1) * (n1_trust / n1_total)
    total_score += n1_score
    calculation_parts.append(f"(-1 / ({n1_total} + 1) * ({n1_trust} / {n1_total}))")

    # Calculate N2 contribution (Negative, mixed, distrust > trust)
    n2_total = user_data['N2']['total']
    n2_trust = user_data['N2']['trust']
    # Rule 3: 1.5x weight
    n2_score = 1.5 * (-1 / (n2_total + 1) * (n2_trust / n2_total))
    total_score += n2_score
    calculation_parts.append(f"(-1.5 * (1 / ({n2_total} + 1) * ({n2_trust} / {n2_total})))")
    
    print("Calculating B's importance (FAT) score:")
    final_equation_str = " + ".join(calculation_parts)
    print(f"\nFinal Equation:\nB's score = {final_equation_str}")
    
    # Print the result of each part
    print(f"\n= {p1_score:.3f} + {p2_score:.3f} + {p3_score:.3f} + {n1_score:.3f} + {n2_score:.3f}")
    
    # Print the final result
    print(f"\n= {total_score:.4f}")
    
    print("\nThe calculated score is ~0.3214, which is closest to 0.35.")

calculate_fat_score()
<<<A>>>