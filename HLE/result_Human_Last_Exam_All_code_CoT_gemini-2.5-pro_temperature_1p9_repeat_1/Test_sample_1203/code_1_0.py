def solve_chi_square_problem():
    """
    Calculates the chi-square statistic for each answer choice to determine
    which is most likely to lead to the rejection of the hypothesis of independent assortment.
    """
    # Expected phenotypic ratio for a Tt Rr Yy x Tt Rr Yy cross
    expected_ratio = [27, 9, 9, 9, 3, 3, 3, 1]
    total_ratio_parts = sum(expected_ratio)  # 64

    # Observed counts from the five answer choices
    observed_data = {
        "A": [140, 10, 10, 10, 10, 10, 10, 100],
        "B": [180, 0, 0, 0, 60, 0, 0, 60],
        "C": [144, 45, 45, 16, 52, 16, 16, 16],
        "D": [150, 60, 50, 40, 30, 40, 30, 50],
        "E": [0, 180, 0, 0, 0, 180, 0, 0]
    }
    
    phenotype_names = [
        "Tall/Round/Yellow", "Tall/Round/Green", "Tall/Wrinkled/Yellow", "Tall/Wrinkled/Green",
        "Dwarf/Round/Yellow", "Dwarf/Round/Green", "Dwarf/Wrinkled/Yellow", "Dwarf/Wrinkled/Green"
    ]

    print("The null hypothesis predicts a 27:9:9:9:3:3:3:1 phenotypic ratio.\n")
    print("Calculating the chi-square value for each option...\n")
    
    results = {}

    for option, observed_counts in observed_data.items():
        total_observed = sum(observed_counts)
        
        # Calculate expected counts based on the total observed for this option
        expected_counts = [(count / total_ratio_parts) * total_observed for count in expected_ratio]
        
        chi_square_value = 0
        equation_terms = []

        for i in range(len(observed_counts)):
            O = observed_counts[i]
            E = expected_counts[i]
            
            # Prevent division by zero if an expected count is zero
            if E == 0:
                # This would mean an infinite contribution to chi-square if O > 0
                term_value = float('inf') if O > 0 else 0
            else:
                term_value = ((O - E)**2) / E
            
            chi_square_value += term_value
            # Build the string for each term in the equation, showing the numbers
            term_str = f"({O} - {E:.2f})^2 / {E:.2f}"
            equation_terms.append(term_str)

        # Store result for final comparison
        results[option] = chi_square_value
        
        # Print the detailed calculation for the current option
        print(f"--- Option {option} ---")
        print(f"Total Offspring: {total_observed}")
        print("Calculation: sum of (Observed - Expected)^2 / Expected for each phenotype")
        print(f"χ² = {' + '.join(equation_terms)}")
        print(f"Result for Option {option}: χ² = {chi_square_value:.2f}\n")

    # Find the option with the maximum chi-square value
    max_chi_option = max(results, key=results.get)

    print("--- Conclusion ---")
    print("A higher chi-square value indicates a greater deviation from the expected results.")
    print("The dataset most likely to lead to the rejection of the null hypothesis is the one with the largest chi-square value.")
    print(f"Comparing the final values: {results}")
    print(f"Option {max_chi_option} has the highest chi-square value and is the correct answer.")

# Run the analysis
solve_chi_square_problem()
<<<E>>>