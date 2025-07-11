def solve_chi_square_problem():
    """
    Calculates the chi-square value for each answer choice to determine which one
    is most likely to lead to the rejection of the null hypothesis of independent assortment.
    """

    # The 8 phenotypes in order are:
    # 1. Tall, round, yellow
    # 2. Tall, round, green
    # 3. Tall, wrinkled, yellow
    # 4. Tall, wrinkled, green
    # 5. dwarf, round, yellow
    # 6. dwarf, round, green
    # 7. dwarf, wrinkled, yellow
    # 8. dwarf, wrinkled, green

    # Observed counts from the answer choices
    choices = {
        'A': [140, 10, 10, 10, 10, 10, 10, 100],
        'B': [180, 0, 0, 0, 60, 0, 0, 60],
        'C': [144, 45, 45, 16, 52, 16, 16, 16],
        'D': [150, 60, 50, 40, 30, 40, 30, 50],
        'E': [0, 180, 0, 0, 0, 180, 0, 0]
    }

    # Expected phenotypic ratio for a trihybrid cross with independent assortment
    ratios = [27, 9, 9, 3, 9, 3, 3, 1]
    total_ratio_parts = sum(ratios)  # This is 64

    # Variables to store the results for the best choice
    max_chi_square = -1
    best_choice_label = ''
    best_choice_details = {}

    print("Analyzing which set of offspring phenotypes most deviates from the expected 27:9:9:3:9:3:3:1 ratio.\n")

    # Loop through each choice to calculate its chi-square value
    for label, observed_counts in choices.items():
        total_n = sum(observed_counts)

        # Calculate expected counts based on the total offspring
        expected_counts = [(part / total_ratio_parts) * total_n for part in ratios]

        # Calculate the chi-square value
        chi_square_value = 0
        for obs, exp in zip(observed_counts, expected_counts):
            # This check prevents division by zero if an expected count was 0
            if exp > 0:
                chi_square_value += ((obs - exp)**2) / exp

        # Check if this is the largest chi-square value found so far
        if chi_square_value > max_chi_square:
            max_chi_square = chi_square_value
            best_choice_label = label
            best_choice_details = {
                'observed': observed_counts,
                'expected': expected_counts,
                'chi_square': chi_square_value,
                'total_n': total_n
            }
        
    # --- Final Output ---
    print(f"The analysis shows that Choice {best_choice_label} has the highest chi-square value ({best_choice_details['chi_square']:.2f}), indicating the largest deviation from the expected ratio.")
    print("This combination would be the most likely to lead to the rejection of the null hypothesis of independent assortment.\n")
    print(f"Here is the detailed calculation for Choice {best_choice_label}:\n")
    
    obs = best_choice_details['observed']
    exp = best_choice_details['expected']

    print(f"Observed Counts: {obs}")
    print(f"Total Offspring: {best_choice_details['total_n']}")
    print(f"Expected Counts: {[round(e, 2) for e in exp]}\n")

    print("Chi-Square Calculation: χ² = Σ [ (Observed - Expected)² / Expected ]\n")
    
    # Build and print the full equation string
    equation_parts = []
    for o, e in zip(obs, exp):
        equation_parts.append(f"({o} - {round(e, 2)})^2 / {round(e, 2)}")
    
    print("χ² = " + " + ".join(equation_parts))
    print(f"\nFinal Chi-Square Value for Choice {best_choice_label} = {best_choice_details['chi_square']:.2f}")


solve_chi_square_problem()
<<<E>>>