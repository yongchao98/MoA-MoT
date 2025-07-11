def calculate_fat_score():
    """
    Calculates the Fringe of Absolute Trust (FAT) score for user B based on the provided data and rules.
    """
    # Step 1: Parse and store data for each user connected to B.
    users = {
        'P1': {'trust': 6, 'distrust': 1, 'total': 7},
        'P2': {'trust': 4, 'distrust': 2, 'total': 6},
        'P3': {'trust': 2, 'distrust': 2, 'total': 4},
        'N1': {'trust': 3, 'distrust': 3, 'total': 6},
        'N2': {'trust': 1, 'distrust': 3, 'total': 4},
    }

    positive_edges_from = ['P1', 'P2', 'P3']
    negative_edges_from = ['N1', 'N2']

    total_score = 0
    equation_parts = []

    # Step 2: Calculate contributions from positive edges.
    for user_id in positive_edges_from:
        user_data = users[user_id]
        contribution = 1 / (user_data['total'] + 1)
        total_score += contribution
        equation_parts.append(f"(1 / ({user_data['total']} + 1))")

    # Step 3: Calculate contributions from negative edges.
    for user_id in negative_edges_from:
        user_data = users[user_id]
        
        # Rule 2: Calculate base negative contribution
        contribution = -1 / (user_data['total'] + 1) * (user_data['trust'] / user_data['total'])
        equation_part = f"(-1 / ({user_data['total']} + 1) * ({user_data['trust']} / {user_data['total']}))"

        # Rule 3: Check for 1.5x multiplier
        if user_data['distrust'] > user_data['trust']:
            contribution *= 1.5
            equation_part = f"({equation_part} * 1.5)"

        total_score += contribution
        equation_parts.append(equation_part)

    # Step 4: Display the results
    print("B's importance score calculation:")
    # The problem asks to output each number, so we replace the inner calculations with the actual terms.
    p1_term = f"1 / 8"
    p2_term = f"1 / 7"
    p3_term = f"1 / 5"
    n1_term = f"-1 / 14" # from -1/7 * 3/6
    n2_term = f"-0.075"   # from (-1/5 * 1/4) * 1.5 = -0.05 * 1.5
    
    # We will build the final equation string with all numbers.
    # The structure of the equation is shown step-by-step for clarity.
    p1_val = 1/(7+1)
    p2_val = 1/(6+1)
    p3_val = 1/(4+1)
    n1_val = -1/(6+1)*(3/6)
    n2_val = (-1/(4+1)*(1/4))*1.5
    
    print(f"Contribution from P1 = 1 / ({users['P1']['total']} + 1) = {p1_val:.4f}")
    print(f"Contribution from P2 = 1 / ({users['P2']['total']} + 1) = {p2_val:.4f}")
    print(f"Contribution from P3 = 1 / ({users['P3']['total']} + 1) = {p3_val:.4f}")
    print(f"Contribution from N1 = -1 / ({users['N1']['total']} + 1) * ({users['N1']['trust']} / {users['N1']['total']}) = {n1_val:.4f}")
    print(f"Contribution from N2 = (-1 / ({users['N2']['total']} + 1) * ({users['N2']['trust']} / {users['N2']['total']})) * 1.5 = {n2_val:.4f}")

    print("\nFinal Equation:")
    print(f"Score = {p1_val:.4f} + {p2_val:.4f} + {p3_val:.4f} + ({n1_val:.4f}) + ({n2_val:.4f})")
    
    print(f"\nFinal Score = {total_score:.4f}")

    # Step 5: Find the closest option
    options = {'A': 0.35, 'B': 0.45, 'C': 0.55, 'D': 0.64}
    closest_option = min(options, key=lambda k: abs(options[k] - total_score))
    
    print(f"The calculated score {total_score:.4f} is closest to option {closest_option}: {options[closest_option]}")
    return closest_option


if __name__ == "__main__":
    final_answer = calculate_fat_score()
    print(f"\nFinal Answer Option: {final_answer}")
    print(f"<<<{final_answer}>>>")
