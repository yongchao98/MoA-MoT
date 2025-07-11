def calculate_fat_score():
    """
    Calculates the FAT score for user B based on the provided data and rules.
    """
    # Data for users connected to B
    # Format: { 'name': {'total': int, 'trust': int, 'distrust': int} }
    trust_incoming = {
        'P1': {'total': 7, 'trust': 6, 'distrust': 1},
        'P2': {'total': 6, 'trust': 4, 'distrust': 2},
        'P3': {'total': 4, 'trust': 2, 'distrust': 2}
    }

    distrust_incoming = {
        'N1': {'total': 6, 'trust': 3, 'distrust': 3},
        'N2': {'total': 4, 'trust': 1, 'distrust': 3}
    }

    total_score = 0
    equation_parts = []

    print("Calculating positive contributions (Rule 1: 1 / (total + 1)):\n")
    # Calculate scores from positive (trust) edges
    for name, data in trust_incoming.items():
        total = data['total']
        score = 1 / (total + 1)
        total_score += score
        part = f"1/({total}+1)"
        equation_parts.append(part)
        print(f"Contribution from {name}: {part} = {score:.4f}")

    print("\nCalculating negative contributions (Rule 2: -1 / (total + 1) * (trust/total)):\n")
    # Calculate scores from negative (distrust) edges
    for name, data in distrust_incoming.items():
        total = data['total']
        trust = data['trust']
        # Apply Rule 2 for negative edges with mixed ratings.
        # As discussed in the plan, we interpret this rule to be the definitive one for N1 and N2
        # as they both have mixed ratings.
        score = -1 / (total + 1) * (trust / total)
        total_score += score
        part = f"(-1/({total}+1) * ({trust}/{total}))"
        equation_parts.append(part)
        print(f"Contribution from {name}: {part} = {score:.4f}")

    print("\n--- Final Calculation ---")
    final_equation = " + ".join(equation_parts).replace("+ (-", "- ")
    print(f"Final Equation: {final_equation}")
    print(f"Total Score = {total_score:.4f}")

    # Determine the closest answer
    options = {'A': 0.35, 'B': 0.45, 'C': 0.55, 'D': 0.64}
    closest_option = min(options, key=lambda x: abs(options[x] - total_score))
    print(f"\nThe calculated score {total_score:.4f} is closest to option {closest_option}: {options[closest_option]}")

calculate_fat_score()
<<<A>>>