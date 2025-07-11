import scipy.stats

def solve_chi_square_problem():
    """
    Calculates the Chi-Square value for each answer choice to determine which
    is most likely to lead to the rejection of the null hypothesis of independent assortment.
    """
    # Phenotype labels in the order of the ratio
    phenotypes = [
        "Tall, round, yellow", "Tall, round, green", "Tall, wrinkled, yellow",
        "Tall, wrinkled, green", "Dwarf, round, yellow", "Dwarf, round, green",
        "Dwarf, wrinkled, yellow", "Dwarf, wrinkled, green"
    ]

    # Expected phenotypic ratio parts for a trihybrid cross (27:9:9:3:9:3:3:1)
    expected_ratio_parts = [27, 9, 9, 3, 9, 3, 3, 1]
    total_ratio_parts = sum(expected_ratio_parts)

    # Observed counts from the answer choices
    observed_data = {
        'A': [140, 10, 10, 10, 10, 10, 10, 100],
        'B': [180, 0, 0, 0, 60, 0, 0, 60],
        'C': [144, 45, 45, 16, 52, 16, 16, 16],
        'D': [150, 60, 50, 40, 30, 40, 30, 50],
        'E': [0, 180, 0, 0, 0, 180, 0, 0]
    }

    chi_square_results = {}
    print("Calculating Chi-Square (χ²) values for each choice...\n")

    # Calculate chi-square for each choice
    for choice, observed_counts in observed_data.items():
        total_observed = sum(observed_counts)
        chi_square_value = 0
        for i in range(len(observed_counts)):
            o = observed_counts[i]
            e = (total_observed * expected_ratio_parts[i]) / total_ratio_parts
            if e > 0:
                chi_square_value += ((o - e) ** 2) / e
        chi_square_results[choice] = chi_square_value
        print(f"Choice {choice}: Total Offspring = {total_observed}, χ² = {chi_square_value:.2f}")

    # Determine which choice is most likely to lead to rejection
    most_likely_choice = max(chi_square_results, key=chi_square_results.get)
    max_chi_square = chi_square_results[most_likely_choice]

    print("\n--- Conclusion ---")
    print(f"The highest Chi-Square value is {max_chi_square:.2f}, corresponding to Choice {most_likely_choice}.")
    print("This indicates the largest deviation from the expected ratio, making it the most likely to lead to a rejection of the hypothesis of independent assortment.\n")
    
    # Show the detailed calculation for the chosen answer
    print(f"--- Detailed Calculation for Choice {most_likely_choice} ---")
    observed = observed_data[most_likely_choice]
    total_observed = sum(observed)
    
    # To display the full equation, we format it as a string
    # χ² = Σ [ (Observed - Expected)² / Expected ]
    equation_str_list = []
    term_values = []

    for i in range(len(observed)):
        o = observed[i]
        e = (total_observed * expected_ratio_parts[i]) / total_ratio_parts
        term_value = ((o - e) ** 2) / e
        term_values.append(f"{term_value:.2f}")
        # Build the string for the equation with numbers
        equation_str_list.append(f"({o} - {e:.2f})² / {e:.2f}")

    print("The formula is: χ² = Σ [ (Observed - Expected)² / Expected ]\n")
    print("Plugging in the numbers for Choice E:")
    # Print the equation with all numbers filled in
    print("χ² = " + " + ".join(equation_str_list))
    # Print the value of each term and the final sum
    print("\nχ² = " + " + ".join(term_values) + f" = {max_chi_square:.2f}")

# Execute the function
solve_chi_square_problem()