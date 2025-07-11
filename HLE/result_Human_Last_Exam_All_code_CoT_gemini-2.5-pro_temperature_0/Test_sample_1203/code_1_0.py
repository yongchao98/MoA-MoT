import sys
# It's better to use the math module for infinity if needed, but it's not required here.
# import math 

def solve_chi_square_problem():
    """
    Calculates the chi-square value for different sets of genetic cross data
    to determine which one most likely leads to the rejection of the null hypothesis
    of independent assortment.
    """
    
    # Phenotype order is consistent with the 27:9:9:3:9:3:3:1 ratio:
    # 1. Tall, round, yellow
    # 2. Tall, round, green
    # 3. Tall, wrinkled, yellow
    # 4. Tall, wrinkled, green
    # 5. dwarf, round, yellow
    # 6. dwarf, round, green
    # 7. dwarf, wrinkled, yellow
    # 8. dwarf, wrinkled, green

    # Observed counts for each answer choice from the problem
    observed_data = {
        "A": [140, 10, 10, 10, 10, 10, 10, 100],
        "B": [180, 0, 0, 0, 60, 0, 0, 60],
        "C": [144, 45, 45, 16, 52, 16, 16, 16],
        "D": [150, 60, 50, 40, 30, 40, 30, 50],
        "E": [0, 180, 0, 0, 0, 180, 0, 0]
    }

    # Expected Mendelian ratio for a trihybrid cross
    expected_ratio = [27, 9, 9, 3, 9, 3, 3, 1]
    total_ratio_parts = sum(expected_ratio)

    results = {}
    
    print("--- Chi-Square Test Calculations ---")
    # Loop through each answer choice to calculate its chi-square value
    for option, observed_counts in observed_data.items():
        total_observed = sum(observed_counts)
        
        # Calculate expected counts based on the total observed offspring
        expected_counts = []
        for part in expected_ratio:
            expected_counts.append((part / total_ratio_parts) * total_observed)
        
        # Calculate the chi-square value
        chi_square_value = 0
        for o, e in zip(observed_counts, expected_counts):
            # The chi-square test is not well-defined if an expected value is zero.
            # This is not an issue here, as all expected counts will be > 0.
            if e == 0:
                # This case should not be reached with the given ratios.
                # If it were, it would imply an infinite chi-square value if o > 0.
                # We can skip this for safety, but it won't affect this problem.
                continue
            chi_square_value += ((o - e)**2) / e
            
        results[option] = chi_square_value
        print(f"Option {option}: Total Offspring = {total_observed}, Chi-Square (χ²) = {chi_square_value:.2f}")

    # Identify the option with the highest chi-square value
    most_likely_rejection_option = max(results, key=results.get)
    
    print("\n--- Conclusion ---")
    print(f"The combination that would most likely lead to rejection of the hypothesis is Option {most_likely_rejection_option},")
    print("as it has the highest chi-square value, indicating the largest deviation from the expected 27:9:9:3:9:3:3:1 ratio.")

    # Display the full calculation for the winning option
    winning_observed = observed_data[most_likely_rejection_option]
    winning_total = sum(winning_observed)
    winning_expected = [(part / total_ratio_parts) * winning_total for part in expected_ratio]
    winning_value = results[most_likely_rejection_option]

    equation_parts = []
    for o, e in zip(winning_observed, winning_expected):
        # Format each term of the sum: (Observed - Expected)^2 / Expected
        equation_parts.append(f"({o} - {e:.2f})^2 / {e:.2f}")

    final_equation_str = " + ".join(equation_parts)

    print("\nThe detailed chi-square calculation for this option is:")
    print(f"χ² = {final_equation_str}")
    print(f"χ² = {winning_value:.2f}")

# Execute the function
solve_chi_square_problem()