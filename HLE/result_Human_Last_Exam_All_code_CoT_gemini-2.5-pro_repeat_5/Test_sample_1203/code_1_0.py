import scipy.stats

def solve_chi_square_problem():
    """
    Calculates the chi-square statistic for each answer choice to determine
    which is most likely to lead to the rejection of the hypothesis of
    independent assortment.
    """
    # Phenotypic classes in order:
    # 1. Tall, round, yellow
    # 2. Tall, round, green
    # 3. Tall, wrinkled, yellow
    # 4. Tall, wrinkled, green
    # 5. dwarf, round, yellow
    # 6. dwarf, round, green
    # 7. dwarf, wrinkled, yellow
    # 8. dwarf, wrinkled, green
    expected_ratio = [27, 9, 9, 3, 9, 3, 3, 1]
    total_ratio = sum(expected_ratio)
    
    # Observed counts from the answer choices
    options = {
        "A": [140, 10, 10, 10, 10, 10, 10, 100],
        "B": [180, 0, 0, 0, 60, 0, 0, 60],
        "C": [144, 45, 45, 16, 52, 16, 16, 16],
        "D": [150, 60, 50, 40, 30, 40, 30, 50],
        "E": [0, 180, 0, 0, 0, 180, 0, 0]
    }
    
    # Degrees of freedom = number of categories - 1
    df = len(expected_ratio) - 1
    # Significance level
    alpha = 0.05
    # Critical value from chi-square distribution
    critical_value = scipy.stats.chi2.ppf(1 - alpha, df)

    print(f"The null hypothesis is that the genes assort independently, following a {':'.join(map(str, expected_ratio))} ratio.")
    print(f"We will reject this hypothesis if the calculated Chi-Square value is greater than the critical value.")
    print(f"Degrees of Freedom (df) = {df}")
    print(f"Significance Level (alpha) = {alpha}")
    print(f"Critical Chi-Square Value = {critical_value:.2f}\n")

    results = {}

    for option, observed_counts in options.items():
        total_observed = sum(observed_counts)
        chi_square_sum = 0
        
        # Build the equation string for printing
        equation_str = f"Option {option}: Chi-Square = "
        terms = []
        
        for i in range(len(observed_counts)):
            observed = observed_counts[i]
            expected = (expected_ratio[i] / total_ratio) * total_observed
            
            # Avoid division by zero if expected is 0 (though not possible with this ratio)
            if expected == 0:
                component = 0
            else:
                component = ((observed - expected)**2) / expected
            
            chi_square_sum += component
            terms.append(f"(({observed} - {expected:.2f})^2 / {expected:.2f})")
            
        equation_str += " + ".join(terms)
        equation_str += f" = {chi_square_sum:.2f}"
        
        print(equation_str)
        results[option] = chi_square_sum
        
    # Find the option with the highest chi-square value
    most_likely_rejection = max(results, key=results.get)
    
    print("\n--- Analysis ---")
    print("The chi-square statistic measures the deviation of observed counts from expected counts.")
    print("A higher chi-square value indicates a larger deviation from the expected ratio.")
    print(f"Comparing the calculated values, Option {most_likely_rejection} has the highest Chi-Square value ({results[most_likely_rejection]:.2f}).")
    print(f"This value is significantly greater than the critical value of {critical_value:.2f}, making it the most likely to lead to the rejection of the hypothesis of independent assortment.")
    
    return most_likely_rejection

final_answer = solve_chi_square_problem()
print(f"\n<<<E>>>")