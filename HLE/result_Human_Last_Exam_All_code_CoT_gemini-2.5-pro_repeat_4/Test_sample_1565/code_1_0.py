def calculate_fat_score():
    """
    Calculates the FAT score for user B based on the provided data and rules.
    """
    # Data for each user connecting to B
    # Format: [type, total_ratings, trust_ratings, distrust_ratings]
    users = {
        'P1': ['positive', 7, 6, 1],
        'P2': ['positive', 6, 4, 2],
        'P3': ['positive', 4, 2, 2],
        'N1': ['negative', 6, 3, 3],
        'N2': ['negative', 4, 1, 3]
    }

    total_score = 0
    equation_parts = []
    
    print("Calculating FAT score for user B...")

    # P1
    p1_total = users['P1'][1]
    p1_score = 1 / (p1_total + 1)
    total_score += p1_score
    equation_parts.append(str(p1_score))
    print(f"P1 (Positive): Contribution = 1/({p1_total} + 1) = {p1_score}")

    # P2
    p2_total = users['P2'][1]
    p2_score = 1 / (p2_total + 1)
    total_score += p2_score
    equation_parts.append(str(p2_score))
    print(f"P2 (Positive): Contribution = 1/({p2_total} + 1) = {p2_score}")

    # P3
    p3_total = users['P3'][1]
    p3_score = 1 / (p3_total + 1)
    total_score += p3_score
    equation_parts.append(str(p3_score))
    print(f"P3 (Positive): Contribution = 1/({p3_total} + 1) = {p3_score}")

    # N1
    n1_total = users['N1'][1]
    n1_trust = users['N1'][2]
    n1_distrust = users['N1'][3]
    n1_score = -1 / (n1_total + 1) * (n1_trust / n1_total)
    if n1_distrust > n1_trust:
        n1_score *= 1.5
        print(f"N1 (Negative): Contribution = (-1/({n1_total} + 1) * ({n1_trust}/{n1_total})) * 1.5 = {n1_score}")
    else:
        print(f"N1 (Negative): Contribution = -1/({n1_total} + 1) * ({n1_trust}/{n1_total}) = {n1_score}")
    total_score += n1_score
    equation_parts.append(str(n1_score))

    # N2
    n2_total = users['N2'][1]
    n2_trust = users['N2'][2]
    n2_distrust = users['N2'][3]
    n2_score = -1 / (n2_total + 1) * (n2_trust / n2_total)
    if n2_distrust > n2_trust:
        n2_score *= 1.5
        print(f"N2 (Negative): Contribution = (-1/({n2_total} + 1) * ({n2_trust}/{n2_total})) * 1.5 = {n2_score}")
    else:
        print(f"N2 (Negative): Contribution = -1/({n2_total} + 1) * ({n2_trust}/{n2_total}) = {n2_score}")
    total_score += n2_score
    equation_parts.append(str(n2_score))

    # Final result
    final_equation = " + ".join(equation_parts)
    print(f"\nFinal Equation: {final_equation}")
    print(f"B's Importance Score = {total_score}")
    
    options = {'A': 0.35, 'B': 0.45, 'C': 0.55, 'D': 0.64}
    closest_option = min(options, key=lambda k: abs(options[k] - total_score))
    print(f"\nThe calculated score {total_score:.4f} is closest to option {closest_option}) {options[closest_option]}.")


calculate_fat_score()