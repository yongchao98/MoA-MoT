import scipy.stats

def solve_chi_square():
    """
    Calculates the chi-square value for each answer choice to determine which one
    would most likely lead to the rejection of the hypothesis of independent assortment.
    """
    # Phenotype order:
    # 1. Tall, round, yellow
    # 2. Tall, round, green
    # 3. Tall, wrinkled, yellow
    # 4. Tall, wrinkled, green
    # 5. dwarf, round, yellow
    # 6. dwarf, round, green
    # 7. dwarf, wrinkled, yellow
    # 8. dwarf, wrinkled, green
    
    # Expected phenotypic ratio for a trihybrid cross (Tt Rr Yy x Tt Rr Yy)
    expected_ratio_parts = [27, 9, 9, 3, 9, 3, 3, 1]
    total_ratio_parts = sum(expected_ratio_parts) # Should be 64

    # Observed counts from the answer choices
    choices = {
        'A': [140, 10, 10, 10, 10, 10, 10, 100],
        'B': [180, 0, 0, 0, 60, 0, 0, 60],
        'C': [144, 45, 45, 16, 52, 16, 16, 16],
        'D': [150, 60, 50, 40, 30, 40, 30, 50],
        'E': [0, 180, 0, 0, 0, 180, 0, 0]
    }
    
    results = {}

    print("--- Chi-Square Test for Independent Assortment ---\n")
    
    for choice, observed_counts in choices.items():
        total_offspring = sum(observed_counts)
        expected_counts = [(total_offspring * part) / total_ratio_parts for part in expected_ratio_parts]
        
        chi_square_value = 0
        for i in range(len(observed_counts)):
            observed = observed_counts[i]
            expected = expected_counts[i]
            # To avoid division by zero if an expected count were zero (not the case here)
            if expected == 0:
                continue
            chi_square_value += ((observed - expected)**2) / expected
        
        results[choice] = chi_square_value
        print(f"Choice {choice}:")
        print(f"  Observed counts: {observed_counts}")
        print(f"  Total offspring: {total_offspring}")
        print(f"  Chi-square (χ²) value = {chi_square_value:.2f}\n")

    # Determine degrees of freedom and the critical value
    df = len(expected_ratio_parts) - 1
    alpha = 0.05
    critical_value = scipy.stats.chi2.ppf(1 - alpha, df)

    # Find the choice with the highest chi-square value
    most_likely_rejected = max(results, key=results.get)

    print("--- Conclusion ---")
    print(f"The degrees of freedom (df) is {df}.")
    print(f"The critical χ² value for a significance level of {alpha} is {critical_value:.2f}.")
    print("\nA calculated χ² value greater than this critical value leads to rejection of the null hypothesis.")
    print(f"\nThe combination of phenotypes from Choice '{most_likely_rejected}' gives the highest χ² value of {results[most_likely_rejected]:.2f}.")
    print(f"This value is significantly larger than the critical value, making it the most likely to be rejected.\n")
    
    # Display the final equation for the winning choice
    print(f"--- Calculation for the winning Choice '{most_likely_rejected}' ---")
    total_offspring_winner = sum(choices[most_likely_rejected])
    expected_counts_winner = [(total_offspring_winner * part) / total_ratio_parts for part in expected_ratio_parts]
    observed_counts_winner = choices[most_likely_rejected]
    
    equation_parts = []
    for obs, exp in zip(observed_counts_winner, expected_counts_winner):
         equation_parts.append(f"({obs} - {exp:.2f})²/ {exp:.2f}")

    print("χ² = " + " + ".join(equation_parts))
    print(f"χ² = {results[most_likely_rejected]:.2f}")


solve_chi_square()
<<<E>>>