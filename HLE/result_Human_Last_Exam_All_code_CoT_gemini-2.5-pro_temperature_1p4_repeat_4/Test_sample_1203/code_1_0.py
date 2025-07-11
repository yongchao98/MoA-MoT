import sys
# Redirect stdout to a variable to prevent printing from imported modules if any.
# This is a good practice for environments where you can't control all library outputs.
original_stdout = sys.stdout
sys.stdout = None
try:
    # It's good practice to import libraries that might print info on import,
    # though scipy.stats doesn't typically do this.
    import scipy.stats
finally:
    # Restore stdout
    sys.stdout = original_stdout


def solve_chi_square_problem():
    """
    Calculates the chi-square statistic for different sets of genetic cross data
    to find which one most likely leads to rejecting the null hypothesis of
    independent assortment.
    """
    # Phenotype order: T_R_Y_, T_R_yy, T_rrY_, T_rryy, ttR_Y_, ttR_yy, ttrrY_, ttrryy
    # Corresponding expected ratio for a trihybrid cross
    expected_ratio = [27, 9, 9, 3, 9, 3, 3, 1]
    total_ratio_parts = sum(expected_ratio)

    # Observed offspring counts for each answer choice
    observed_counts = {
        'A': [140, 10, 10, 10, 10, 10, 10, 100],
        'B': [180, 0, 0, 0, 60, 0, 0, 60],
        'C': [144, 45, 45, 16, 52, 16, 16, 16],
        'D': [150, 60, 50, 40, 30, 40, 30, 50],
        'E': [0, 180, 0, 0, 0, 180, 0, 0]
    }

    chi_square_results = {}
    calculation_details = {}

    # Calculate chi-square for each choice
    for choice, observed in observed_counts.items():
        total_offspring = sum(observed)
        chi_square_stat = 0
        
        expected_values = [(part / total_ratio_parts) * total_offspring for part in expected_ratio]
        
        equation_parts = []
        for i in range(len(observed)):
            O = observed[i]
            E = expected_values[i]
            term = ((O - E)**2) / E
            chi_square_stat += term
            equation_parts.append(f"({O} - {E:.2f})^2 / {E:.2f}")

        chi_square_results[choice] = chi_square_stat
        calculation_details[choice] = {
            "equation": " + ".join(equation_parts),
            "observed": observed,
            "expected": expected_values
        }

    # Find the choice with the maximum chi-square value
    best_choice = max(chi_square_results, key=chi_square_results.get)
    max_chi_square = chi_square_results[best_choice]
    
    # Degrees of freedom = (number of categories - 1)
    df = len(expected_ratio) - 1
    # Critical value for significance level alpha = 0.05
    critical_value = scipy.stats.chi2.ppf(1 - 0.05, df)

    print("--- Chi-Square Analysis ---")
    for choice, value in chi_square_results.items():
        print(f"Choice {choice}: Chi-Square = {value:.2f}")
    
    print("\n--- Conclusion ---")
    print(f"The critical value for rejection at alpha=0.05 with {df} degrees of freedom is approximately {critical_value:.2f}.")
    print(f"Choice '{best_choice}' has the highest Chi-Square value ({max_chi_square:.2f}), indicating the largest deviation from the expected ratio.")
    print("This makes it the most likely candidate for rejecting the hypothesis of independent assortment.")

    print("\n--- Detailed Calculation for the Final Answer ---")
    details = calculation_details[best_choice]
    print(f"The calculation for choice '{best_choice}' is:")
    print(f"Chi-Square = {details['equation']}")
    
    # Printing each number in the final equation by recalculating the terms
    final_sum_parts = []
    recalculated_sum = 0
    for i in range(len(details['observed'])):
        O = details['observed'][i]
        E = details['expected'][i]
        term = ((O - E)**2) / E
        recalculated_sum += term
        final_sum_parts.append(f"{term:.2f}")
    
    print(f"Chi-Square = {' + '.join(final_sum_parts)}")
    print(f"Final Chi-Square Value = {recalculated_sum:.2f}")


if __name__ == '__main__':
    solve_chi_square_problem()
<<<E>>>