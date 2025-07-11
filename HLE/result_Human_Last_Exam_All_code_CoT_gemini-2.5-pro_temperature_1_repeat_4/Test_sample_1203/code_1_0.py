import scipy.stats

def solve_chi_square():
    """
    Calculates the chi-square value for each answer choice to determine which one
    most strongly contradicts the hypothesis of independent assortment.
    """
    # Phenotype order:
    # 1. Tall, round, yellow (T_ R_ Y_)
    # 2. Tall, round, green (T_ R_ yy)
    # 3. Tall, wrinkled, yellow (T_ rr Y_)
    # 4. Tall, wrinkled, green (T_ rr yy)
    # 5. dwarf, round, yellow (tt R_ Y_)
    # 6. dwarf, round, green (tt R_ yy)
    # 7. dwarf, wrinkled, yellow (tt rr Y_)
    # 8. dwarf, wrinkled, green (tt rr yy)
    
    expected_ratio = [27, 9, 9, 3, 9, 3, 3, 1]
    total_ratio_parts = sum(expected_ratio)

    phenotype_labels = [
        "Tall, round, yellow", "Tall, round, green", "Tall, wrinkled, yellow",
        "Tall, wrinkled, green", "dwarf, round, yellow", "dwarf, round, green",
        "dwarf, wrinkled, yellow", "dwarf, wrinkled, green"
    ]

    observed_data = {
        'A': [140, 10, 10, 10, 10, 10, 10, 100],
        'B': [180, 0, 0, 0, 60, 0, 0, 60],
        'C': [144, 45, 45, 16, 52, 16, 16, 16],
        'D': [150, 60, 50, 40, 30, 40, 30, 50],
        'E': [0, 180, 0, 0, 0, 180, 0, 0]
    }

    results = {}
    
    print("--- Chi-Square Test for Independent Assortment ---")
    print(f"Expected Phenotypic Ratio: {' : '.join(map(str, expected_ratio))}\n")

    for option, observed_counts in observed_data.items():
        total_offspring = sum(observed_counts)
        
        expected_counts = [(part / total_ratio_parts) * total_offspring for part in expected_ratio]
        
        chi_square_value = sum(
            ((obs - exp)**2) / exp for obs, exp in zip(observed_counts, expected_counts) if exp > 0
        )
        
        results[option] = {
            "chi_square": chi_square_value,
            "observed": observed_counts,
            "expected": expected_counts
        }

        print(f"--- Option {option} ---")
        print(f"Observed Counts: {observed_counts}")
        print(f"Total Offspring (N): {total_offspring}")
        print(f"Expected Counts: {[round(c, 2) for c in expected_counts]}")
        print(f"Calculated Chi-Square (Χ²) Value: {chi_square_value:.2f}\n")

    # Determine critical value
    df = len(expected_ratio) - 1
    alpha = 0.05
    critical_value = scipy.stats.chi2.ppf(1 - alpha, df)
    
    print("--- Conclusion ---")
    print(f"Degrees of Freedom (df): {df}")
    print(f"Significance Level (p): {alpha}")
    print(f"Critical Chi-Square Value: {critical_value:.2f}\n")

    # Find the option with the highest chi-square value
    best_option = max(results, key=lambda k: results[k]['chi_square'])
    best_result = results[best_option]

    print(f"The highest chi-square value is found in Option {best_option} (Χ² = {best_result['chi_square']:.2f}).")
    print("This indicates the largest deviation from the expected ratio, making it the most likely to be rejected.")
    
    print("\nThe full chi-square calculation for the winning option is:")
    
    equation_parts = []
    for obs, exp in zip(best_result['observed'], best_result['expected']):
        equation_parts.append(f"({obs} - {exp:.2f})^2 / {exp:.2f}")

    final_equation = " + ".join(equation_parts)
    print(f"Χ² = {final_equation}")
    print(f"Χ² = {best_result['chi_square']:.2f}")

solve_chi_square()