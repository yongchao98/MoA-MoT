import sys

def solve_chi_square_problem():
    """
    Calculates the Chi-Square value for several sets of observed genetic cross data
    to determine which is most likely to lead to a rejection of the null hypothesis.
    """
    
    # The expected phenotypic ratio for a trihybrid cross assuming independent assortment.
    # The order corresponds to:
    # 1. Tall, Round, Yellow (27)
    # 2. Tall, Round, Green (9)
    # 3. Tall, Wrinkled, Yellow (9)
    # 4. Tall, Wrinkled, Green (3)
    # 5. Dwarf, Round, Yellow (9)
    # 6. Dwarf, Round, Green (3)
    # 7. Dwarf, Wrinkled, Yellow (3)
    # 8. Dwarf, Wrinkled, Green (1)
    expected_ratio = [27, 9, 9, 3, 9, 3, 3, 1]
    ratio_sum = sum(expected_ratio)  # This will be 64

    # Observed counts for each answer choice
    observed_counts = {
        'A': [140, 10, 10, 10, 10, 10, 10, 100],
        'B': [180, 0, 0, 0, 60, 0, 0, 60],
        'C': [144, 45, 45, 16, 52, 16, 16, 16],
        'D': [150, 60, 50, 40, 30, 40, 30, 50],
        'E': [0, 180, 0, 0, 0, 180, 0, 0]
    }
    
    chi_square_results = {}

    print("--- Chi-Square Calculation for Each Option ---\n")
    # Redirect stderr to stdout to ensure error messages are captured if any (e.g., division by zero)
    # Note: With these inputs, division by zero in the chi-square formula won't occur
    # because total observed counts > 0, making all expected counts > 0.
    original_stderr = sys.stderr
    sys.stderr = sys.stdout

    try:
        for option, observed in observed_counts.items():
            total_observed = sum(observed)
            
            # Calculate expected counts based on the total observed and the ratio
            expected = [(r / ratio_sum) * total_observed for r in expected_ratio]
            
            # Calculate the chi-square value and the individual terms for the equation
            chi_square_value = 0
            chi_square_terms = []
            
            for o, e in zip(observed, expected):
                # Handle cases where E might be zero to avoid division errors, although not expected here.
                if e == 0:
                    # If expected is 0 and observed is also 0, the contribution is 0.
                    # If expected is 0 and observed is > 0, chi-square is technically infinite.
                    term = 0 if o == 0 else float('inf')
                else:
                    term = (o - e)**2 / e
                chi_square_terms.append(term)
                chi_square_value += term

            chi_square_results[option] = chi_square_value
            
            # Print the detailed breakdown for the current option
            print(f"--- Option {option} (Total: {total_observed}) ---")
            
            # Create the string for the equation with each number
            equation_str = " + ".join([f"{term:.2f}" for term in chi_square_terms])
            
            print(f"Observed: {observed}")
            print(f"Expected: {[f'{e:.2f}' for e in expected]}")
            print(f"χ² = {equation_str}")
            print(f"Total χ² Value: {chi_square_value:.2f}\n")
    
    finally:
        # Restore stderr
        sys.stderr = original_stderr


    # Find the option with the highest Chi-Square value
    max_chi_option = max(chi_square_results, key=chi_square_results.get)
    max_chi_value = chi_square_results[max_chi_option]

    print("--- Conclusion ---")
    print("A higher Chi-Square value indicates a larger deviation from the expected ratio.")
    print("The data set most likely to lead to the rejection of the null hypothesis is the one with the highest Chi-Square value.")
    print(f"\nOption '{max_chi_option}' has the highest Chi-Square value ({max_chi_value:.2f}).")
    print("\nTherefore, it represents the combination of phenotypes that would most likely lead to rejection of the hypothesis of independent assortment.")


# Execute the function
solve_chi_square_problem()
# Final answer based on the calculation result.
print("<<<E>>>")