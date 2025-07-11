import scipy.stats

def solve_chi_square_problem():
    """
    Calculates the chi-square statistic for different sets of observed genetic cross data
    to determine which is most likely to reject the null hypothesis of independent assortment.
    """

    # Phenotypes in order:
    # 1. Tall, round, yellow (T_ R_ Y_)
    # 2. Tall, round, green (T_ R_ yy)
    # 3. Tall, wrinkled, yellow (T_ rr Y_)
    # 4. Tall, wrinkled, green (T_ rr yy)
    # 5. Dwarf, round, yellow (tt R_ Y_)
    # 6. Dwarf, round, green (tt R_ yy)
    # 7. Dwarf, wrinkled, yellow (tt rr Y_)
    # 8. Dwarf, wrinkled, green (tt rr yy)
    
    observed_data = {
        'A': [140, 10, 10, 10, 10, 10, 10, 100],
        'B': [180, 0, 0, 0, 60, 0, 0, 60],
        'C': [144, 45, 45, 16, 52, 16, 16, 16],
        'D': [150, 60, 50, 40, 30, 40, 30, 50],
        'E': [0, 180, 0, 0, 0, 180, 0, 0]
    }

    # Expected phenotypic ratio for a trihybrid cross
    ratio_parts = [27, 9, 9, 3, 9, 3, 3, 1]
    ratio_total = sum(ratio_parts)
    
    results = {}

    print("--- Chi-Square Test Analysis ---\n")
    print(f"Null Hypothesis: Genes assort independently following a {':'.join(map(str, ratio_parts))} ratio.")

    for option, observed in observed_data.items():
        total_offspring = sum(observed)
        
        expected_counts = [(part / ratio_total) * total_offspring for part in ratio_parts]
        
        chi_square_value = 0
        calculation_terms = []
        for i in range(len(observed)):
            O = observed[i]
            E = expected_counts[i]
            # Avoid division by zero if an expected count is zero (not the case here)
            if E == 0:
                # If observed is also 0, term is 0. If observed > 0, chi-square is infinite.
                term = 0 if O == 0 else float('inf')
            else:
                term = (O - E)**2 / E
            chi_square_value += term
            calculation_terms.append(term)

        results[option] = {
            'chi_square': chi_square_value,
            'observed': observed,
            'expected': expected_counts,
            'terms': calculation_terms
        }

        print(f"\n--- Analysis for Option {option} ---")
        print(f"Observed Counts: {observed}")
        print(f"Total Offspring (N): {total_offspring}")
        print(f"Expected Counts: {[round(e, 2) for e in expected_counts]}")
        print(f"Calculated Chi-Square (χ²): {chi_square_value:.2f}")

    # Find the option with the maximum chi-square value
    max_chi_option = max(results, key=lambda k: results[k]['chi_square'])
    max_chi_value = results[max_chi_option]['chi_square']
    
    # Determine critical value for context
    degrees_freedom = len(ratio_parts) - 1
    alpha = 0.05
    critical_value = scipy.stats.chi2.ppf(1 - alpha, degrees_freedom)

    print("\n--- Conclusion ---")
    print(f"The critical value for χ² with {degrees_freedom} degrees of freedom at a significance level of {alpha} is {critical_value:.2f}.")
    print("Rejection of the null hypothesis is most likely for the data set with the largest chi-square value.")
    print(f"\nOption {max_chi_option} has the largest chi-square value of {max_chi_value:.2f}.")
    print("Therefore, this combination of phenotypes would most likely lead to rejection of the hypothesis of independent assortment.")

    # Show the full calculation for the winning option
    print(f"\nFinal Equation for Option {max_chi_option}:")
    calculation_str = []
    obs = results[max_chi_option]['observed']
    exp = results[max_chi_option]['expected']
    
    for i in range(len(obs)):
        calculation_str.append(f"({obs[i]} - {exp[i]:.2f})^2 / {exp[i]:.2f}")
    
    final_equation = " + ".join(calculation_str)
    print(f"χ² = {final_equation}")
    
    term_values = " + ".join([f"{term:.2f}" for term in results[max_chi_option]['terms']])
    print(f"   = {term_values}")
    print(f"   = {max_chi_value:.2f}")


solve_chi_square_problem()
<<<E>>>