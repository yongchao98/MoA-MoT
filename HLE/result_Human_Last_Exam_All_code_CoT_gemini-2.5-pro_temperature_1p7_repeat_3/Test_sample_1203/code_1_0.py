import scipy.stats

def solve_chi_square():
    """
    Calculates the chi-square value for each answer choice to determine which one
    is most likely to lead to the rejection of the hypothesis of independent assortment.
    """
    # Phenotype labels in the standard 27:9:9:3:9:3:3:1 order
    phenotypes = [
        "Tall, round, yellow", "Tall, round, green", "Tall, wrinkled, yellow",
        "Tall, wrinkled, green", "Dwarf, round, yellow", "Dwarf, round, green",
        "Dwarf, wrinkled, yellow", "Dwarf, wrinkled, green"
    ]
    # Corresponding ratio parts for the 9:3:3:1 pattern extended to three genes
    expected_ratio_parts = [27, 9, 9, 3, 9, 3, 3, 1]

    # Observed data from the answer choices
    # Note: The order in the problem statement needs to be matched to the standard ratio order
    # Problem order: T-R-Y, T-R-y, T-rr-Y, T-rr-y, t-R-Y, t-R-y, t-rr-Y, t-rr-y
    # Standard order: T-R-Y, T-R-y, T-rr-Y, T-rr-y, t-R-Y, t-R-y, t-rr-Y, t-rr-y
    # In this case, the problem order is consistent with a common listing format, let's use that order.
    # To be clear, let's reorder the ratio parts to match the problem's text exactly:
    # 1. T-R-Y: 27
    # 2. T-R-y: 9
    # 3. T-rr-Y: 9
    # 4. T-rr-y: 3
    # 5. t-R-Y: 9
    # 6. t-R-y: 3
    # 7. t-rr-Y: 3
    # 8. t-rr-y: 1
    expected_ratio_parts_ordered = [27, 9, 9, 3, 9, 3, 3, 1]

    observed_data = {
        'A': [140, 10, 10, 10, 10, 10, 10, 100],
        'B': [180, 0, 0, 0, 60, 0, 0, 60],
        'C': [144, 45, 45, 16, 52, 16, 16, 16],
        'D': [150, 60, 50, 40, 30, 40, 30, 50],
        'E': [0, 180, 0, 0, 0, 180, 0, 0]
    }

    results = {}
    total_ratio = sum(expected_ratio_parts_ordered)
    df = len(phenotypes) - 1
    p_value = 0.05
    critical_value = scipy.stats.chi2.ppf(1 - p_value, df)

    print(f"Null Hypothesis: The genes assort independently (expected ratio 27:9:9:3:9:3:3:1).")
    print(f"Significance Level (alpha) = {p_value}")
    print(f"Degrees of Freedom (df) = {df}")
    print(f"Critical Chi-Square Value = {critical_value:.2f}\n")


    for option, observed_counts in observed_data.items():
        total_offspring = sum(observed_counts)
        chi_square_value = 0
        
        print(f"--- Calculating for Option {option} ---")
        print(f"Observed Counts: {observed_counts}")
        print(f"Total Offspring: {total_offspring}")
        
        expected_counts = [(part / total_ratio) * total_offspring for part in expected_ratio_parts_ordered]
        
        # We will now print the full chi-square equation for this option.
        equation_parts = []
        for i in range(len(observed_counts)):
            obs = observed_counts[i]
            exp = expected_counts[i]
            # Avoid division by zero if an expected count is zero (though not possible with this ratio)
            if exp == 0:
                # If observed is also 0, contribution is 0. If observed > 0, chi2 is infinite.
                term_val = 0 if obs == 0 else float('inf')
            else:
                term_val = ((obs - exp)**2) / exp
            chi_square_value += term_val
            equation_parts.append(f"({obs} - {exp:.2f})^2 / {exp:.2f}")

        equation_str = " + ".join(equation_parts)
        print(f"Chi-Square Calculation:\nχ² = {equation_str}")
        print(f"Calculated Chi-Square Value for Option {option}: {chi_square_value:.2f}\n")
        results[option] = chi_square_value

    most_likely_rejection = max(results, key=results.get)
    
    print("--- Conclusion ---")
    print("Comparison of Chi-Square Values:")
    for option, value in results.items():
        print(f"Option {option}: χ² = {value:.2f}")

    print(f"\nOption {most_likely_rejection} has the highest chi-square value ({results[most_likely_rejection]:.2f}), "
          "indicating the largest deviation from the expected ratio.")
    print("This option would most likely lead to the rejection of the null hypothesis of independent assortment.")

solve_chi_square()
<<<E>>>