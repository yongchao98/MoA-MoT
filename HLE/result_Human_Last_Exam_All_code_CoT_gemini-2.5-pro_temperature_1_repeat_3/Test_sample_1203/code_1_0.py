import scipy.stats

def solve_chi_square():
    """
    Calculates the chi-square value for each set of observed offspring
    phenotypes to determine which is most likely to lead to the rejection
    of the null hypothesis of independent assortment.
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
    phenotype_labels = [
        "Tall, round, yellow", "Tall, round, green", "Tall, wrinkled, yellow",
        "Tall, wrinkled, green", "dwarf, round, yellow", "dwarf, round, green",
        "dwarf, wrinkled, yellow", "dwarf, wrinkled, green"
    ]

    # Expected phenotypic ratio for a trihybrid cross (Tt Rr Yy x Tt Rr Yy)
    expected_ratio = [27, 9, 9, 3, 9, 3, 3, 1]
    total_ratio_parts = sum(expected_ratio)

    # Observed counts from answer choices
    observed_data = {
        'A': [140, 10, 10, 10, 10, 10, 10, 100],
        'B': [180, 0, 0, 0, 60, 0, 0, 60],
        'C': [144, 45, 45, 16, 52, 16, 16, 16],
        'D': [150, 60, 50, 40, 30, 40, 30, 50],
        'E': [0, 180, 0, 0, 0, 180, 0, 0]
    }

    results = {}

    print("--- Chi-Square Calculations ---")
    print("Hypothesis: Genes assort independently (Expected Ratio 27:9:9:3:9:3:3:1)\n")

    for option, observed_counts in observed_data.items():
        total_observed = sum(observed_counts)
        expected_counts = [(count / total_ratio_parts) * total_observed for count in expected_ratio]
        
        # Chi-square formula: sum of ( (O-E)^2 / E )
        chi_square_value = sum([(o - e)**2 / e for o, e in zip(observed_counts, expected_counts)])
        
        results[option] = chi_square_value
        print(f"Option {option}: Total Offspring = {total_observed}, Chi-Square (χ²) = {chi_square_value:.2f}")

    # Find the option with the highest chi-square value
    most_likely_rejection = max(results, key=results.get)
    
    print("\n--- Conclusion ---")
    print(f"The highest chi-square value is found in Option {most_likely_rejection}.")
    print("A higher chi-square value indicates a greater deviation from the expected ratio.")
    print("Therefore, this option is the most likely to lead to the rejection of the hypothesis of independent assortment.\n")
    
    # Show the detailed calculation for the winning option
    print(f"--- Detailed Calculation for Winning Option ({most_likely_rejection}) ---")
    
    winning_observed = observed_data[most_likely_rejection]
    total_observed_win = sum(winning_observed)
    winning_expected = [(count / total_ratio_parts) * total_observed_win for count in expected_ratio]

    print(f"χ² = Σ [ (Observed - Expected)² / Expected ]\n")
    total_chi_sq = 0
    equation_terms = []
    for i in range(len(winning_observed)):
        O = winning_observed[i]
        E = winning_expected[i]
        term = (O - E)**2 / E
        total_chi_sq += term
        print(f"Term for '{phenotype_labels[i]}':")
        print(f"  ({O:.2f} - {E:.2f})² / {E:.2f} = {term:.2f}\n")
        equation_terms.append(f"{term:.2f}")
    
    print("Final Equation:")
    print(f"χ² = {' + '.join(equation_terms)} = {total_chi_sq:.2f}")


solve_chi_square()
<<<E>>>