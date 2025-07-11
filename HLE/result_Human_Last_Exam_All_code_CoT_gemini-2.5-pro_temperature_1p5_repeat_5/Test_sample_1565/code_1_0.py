import sys

def calculate_fat_score():
    """
    Calculates the Fringe of Absolute Trust (FAT) score for user B based on provided data and rules.
    """

    # Data for users with an incoming edge to B
    # P1, P2, P3 have a positive (trust) edge to B.
    # N1, N2 have a negative (distrust) edge to B.
    users_data = {
        'P1': {'type': 'trust', 'total': 7, 'trust': 6, 'distrust': 1},
        'P2': {'type': 'trust', 'total': 6, 'trust': 4, 'distrust': 2},
        'P3': {'type': 'trust', 'total': 4, 'trust': 2, 'distrust': 2},
        'N1': {'type': 'distrust', 'total': 6, 'trust': 3, 'distrust': 3},
        'N2': {'type': 'distrust', 'total': 4, 'trust': 1, 'distrust': 3}
    }

    total_score = 0
    equation_parts = []
    
    # Calculate contribution for P1 (Trust)
    # Rule 1: Positive edge contributes 1/(total_relationships + 1)
    p1_total = users_data['P1']['total']
    p1_contribution = 1 / (p1_total + 1)
    total_score += p1_contribution
    equation_parts.append(f"1/({p1_total}+1)")

    # Calculate contribution for P2 (Trust)
    # Rule 1
    p2_total = users_data['P2']['total']
    p2_contribution = 1 / (p2_total + 1)
    total_score += p2_contribution
    equation_parts.append(f"1/({p2_total}+1)")

    # Calculate contribution for P3 (Trust)
    # Rule 1
    p3_total = users_data['P3']['total']
    p3_contribution = 1 / (p3_total + 1)
    total_score += p3_contribution
    equation_parts.append(f"1/({p3_total}+1)")

    # Calculate contribution for N1 (Distrust)
    # Rule 2: Negative edge with mixed ratings: -1 /(total+1) x (trust_ratings/total)
    n1_total = users_data['N1']['total']
    n1_trust = users_data['N1']['trust']
    n1_distrust = users_data['N1']['distrust']
    n1_contribution = -1 / (n1_total + 1) * (n1_trust / n1_total)
    # Rule 3: Check for 1.5x multiplier. (3 distrust is not > 3 trust)
    # No multiplier for N1.
    total_score += n1_contribution
    equation_parts.append(f"(-1/({n1_total}+1) * {n1_trust}/{n1_total})")

    # Calculate contribution for N2 (Distrust)
    # Rule 2
    n2_total = users_data['N2']['total']
    n2_trust = users_data['N2']['trust']
    n2_distrust = users_data['N2']['distrust']
    n2_contribution = -1 / (n2_total + 1) * (n2_trust / n2_total)
    # Rule 3: Check for 1.5x multiplier. (3 distrust > 1 trust)
    # Multiplier applies for N2.
    n2_contribution *= 1.5
    total_score += n2_contribution
    equation_parts.append(f"(-1/({n2_total}+1) * {n2_trust}/{n2_total} * 1.5)")

    # Print the full equation
    # To prevent floating point inaccuracies in the display, we show the fractions
    print("B's importance score calculation:")
    final_equation = f"Score = {equation_parts[0]} + {equation_parts[1]} + {equation_parts[2]} + {equation_parts[3]} + {equation_parts[4]}"
    # Replace the formula for N2 to make the multiplication explicit
    final_equation = final_equation.replace(f"(-1/({n2_total}+1) * {n2_trust}/{n2_total} * 1.5)", f"(-1/({n2_total}+1) * {n2_trust}/{n2_total})*1.5")
    
    print(final_equation)

    # Print the numeric breakdown
    print(f"      = (1/8) + (1/7) + (1/5) + (-1/14) + (-1/20)*1.5")
    print(f"      = 0.125 + {1/7:.3f} + 0.200 + ({-1/14:.3f}) + ({-1/20*1.5:.3f})")
    
    # Print the final score
    print(f"\nFinal Score = {total_score:.4f}")

# Execute the calculation
calculate_fat_score()
<<<A>>>