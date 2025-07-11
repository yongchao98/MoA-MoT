import scipy.stats

def solve_chi_square():
    """
    Calculates the chi-square value for each set of offspring phenotypes
    and determines which is most likely to lead to rejection of the null hypothesis.
    """
    options = {
        'A': [140, 10, 10, 10, 10, 10, 10, 100],
        'B': [180, 0, 0, 0, 60, 0, 0, 60],
        'C': [144, 45, 45, 16, 52, 16, 16, 16],
        'D': [150, 60, 50, 40, 30, 40, 30, 50],
        'E': [0, 180, 0, 0, 0, 180, 0, 0]
    }

    # The expected ratio based on the phenotype order in the problem
    expected_ratio_parts = [27, 9, 9, 3, 9, 3, 3, 1]
    
    phenotypes = [
        "tall, round, yellow", "tall, round, green", "tall, wrinkled, yellow",
        "tall, wrinkled, green", "dwarf, round, yellow", "dwarf, round, green",
        "dwarf, wrinkled, yellow", "dwarf, wrinkled, green"
    ]
    
    results = {}
    
    print("Calculating Chi-Square for each option...\n")
    
    for key, observed_counts in options.items():
        total_offspring = sum(observed_counts)
        
        expected_counts = [(p / 64) * total_offspring for p in expected_ratio_parts]
        
        # Using scipy.stats.chisquare to get the value
        # The lambda correction is not needed here as degrees of freedom > 1
        chi2_value, p_value = scipy.stats.chisquare(f_obs=observed_counts, f_exp=expected_counts)
        
        results[key] = chi2_value
        
        print(f"--- Option {key} ---")
        print(f"Observed: {observed_counts}")
        print(f"Total Offspring: {total_offspring}")
        print(f"Expected: {[round(e, 2) for e in expected_counts]}")
        print(f"Chi-Square (χ²) value: {chi2_value:.2f}\n")

    # Find the option with the maximum chi-square value
    max_chi2_option = max(results, key=results.get)
    
    # Degrees of freedom = number of categories - 1
    df = len(expected_ratio_parts) - 1
    # Significance level
    alpha = 0.05
    # Critical value from chi-square distribution table
    critical_value = scipy.stats.chi2.ppf(1 - alpha, df)
    
    print("--- Conclusion ---")
    print(f"The critical value for a significance level of {alpha} with {df} degrees of freedom is approximately {critical_value:.2f}.")
    print("To reject the null hypothesis of independent assortment, the calculated χ² value must be greater than this critical value.")
    print(f"\nOption {max_chi2_option} has the highest chi-square value ({results[max_chi2_option]:.2f}), indicating the largest deviation from the expected 27:9:9:3:9:3:3:1 ratio.")
    print(f"Therefore, the data from Option {max_chi2_option} would most likely lead to the rejection of the hypothesis.\n")

    # Display the final calculation for the winning option
    print(f"Detailed Calculation for Option {max_chi2_option}:")
    observed = options[max_chi2_option]
    total = sum(observed)
    expected = [(p / 64) * total for p in expected_ratio_parts]
    
    equation_terms = []
    for i in range(len(observed)):
        term = f"({observed[i]} - {expected[i]:.2f})^2 / {expected[i]:.2f}"
        equation_terms.append(term)
        
    final_equation = " + ".join(equation_terms)
    print(f"χ² = {final_equation}")
    print(f"χ² = {results[max_chi2_option]:.2f}")


solve_chi_square()
<<<E>>>