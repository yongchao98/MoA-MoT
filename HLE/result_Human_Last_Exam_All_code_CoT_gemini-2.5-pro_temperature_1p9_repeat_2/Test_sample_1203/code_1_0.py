import scipy.stats

def solve_chi_square_problem():
    """
    Calculates the chi-square statistic for each answer choice to determine which is most likely
    to lead to the rejection of the hypothesis of independent assortment.
    """
    # Define the observed counts for each answer choice
    observed_data = {
        'A': [140, 10, 10, 10, 10, 10, 10, 100],
        'B': [180, 0, 0, 0, 60, 0, 0, 60],
        'C': [144, 45, 45, 16, 52, 16, 16, 16],
        'D': [150, 60, 50, 40, 30, 40, 30, 50],
        'E': [0, 180, 0, 0, 0, 180, 0, 0]
    }

    # Expected phenotypic ratio for a trihybrid cross (Tt Rr Yy x Tt Rr Yy)
    expected_ratio_parts = [27, 9, 9, 3, 9, 3, 3, 1]
    total_ratio = sum(expected_ratio_parts)
    expected_proportions = [part / total_ratio for part in expected_ratio_parts]

    results = {}
    print("Calculating Chi-Square values for each option:\n")

    for option, observed_counts in observed_data.items():
        total_observed = sum(observed_counts)
        expected_counts = [prop * total_observed for prop in expected_proportions]

        # Calculate chi-square statistic
        chi_square = 0
        for i in range(len(observed_counts)):
            # The chi-square test is not reliable if an expected count is less than 5.
            # However, for this problem, we calculate it as is.
            if expected_counts[i] == 0:
                # Avoid division by zero, though not expected in this specific problem
                term = 0
            else:
                term = (observed_counts[i] - expected_counts[i])**2 / expected_counts[i]
            chi_square += term

        results[option] = chi_square
        print(f"Option {option}: Total Offspring = {total_observed}, Chi-Square (χ²) = {chi_square:.2f}")

    # Find the option with the maximum chi-square value
    max_chi_option = max(results, key=results.get)
    max_chi_value = results[max_chi_option]
    
    # Degrees of freedom = number of categories - 1
    df = len(expected_ratio_parts) - 1
    # Significance level
    alpha = 0.05
    # Critical value from chi-square distribution
    critical_value = scipy.stats.chi2.ppf(1 - alpha, df)

    print("\n" + "="*50)
    print(f"\nThe critical value for a significance level of {alpha} and {df} degrees of freedom is ~{critical_value:.2f}.")
    print(f"Any Chi-Square value greater than {critical_value:.2f} would lead to a rejection of the null hypothesis.")
    print(f"\nOption '{max_chi_option}' has the highest Chi-Square value ({max_chi_value:.2f}).")
    print("This represents the greatest deviation from the expected 27:9:9:3:9:3:3:1 ratio and is the most likely to lead to rejection of the hypothesis of independent assortment.")

    # Show the full calculation for the winning option
    print("\nFinal Answer Calculation:")
    print(f"For option {max_chi_option}, the Chi-Square value is calculated as follows:")
    
    # Recalculate values for the winning option to print the full equation
    observed_counts = observed_data[max_chi_option]
    total_observed = sum(observed_counts)
    expected_counts = [prop * total_observed for prop in expected_proportions]
    
    equation_parts = []
    for obs, exp in zip(observed_counts, expected_counts):
        equation_parts.append(f"(({obs} - {exp:.2f})^2 / {exp:.2f})")
    
    final_equation = " + ".join(equation_parts)
    print(f"χ² = {final_equation} = {max_chi_value:.2f}")


solve_chi_square_problem()
<<<E>>>