import scipy.stats

def solve_chi_square():
    """
    Calculates the chi-square statistic for each of the five answer choices
    to determine which is most likely to lead to rejection of the null hypothesis
    of independent assortment.
    """

    # Phenotype labels in the order given by the 27:9:9:3:9:3:3:1 expected ratio
    phenotypes = [
        "Tall, round, yellow", "Tall, round, green", "Tall, wrinkled, yellow",
        "Tall, wrinkled, green", "Dwarf, round, yellow", "Dwarf, round, green",
        "Dwarf, wrinkled, yellow", "Dwarf, wrinkled, green"
    ]
    
    # Expected ratio parts for a Tt Rr Yy x Tt Rr Yy cross
    expected_ratio_parts = [27, 9, 9, 3, 9, 3, 3, 1]
    total_ratio_parts = sum(expected_ratio_parts)

    # Observed data from the answer choices
    observed_data = {
        'A': [140, 10, 10, 10, 10, 10, 10, 100],
        'B': [180, 0, 0, 0, 60, 0, 0, 60],
        'C': [144, 45, 45, 16, 52, 16, 16, 16],
        'D': [150, 60, 50, 40, 30, 40, 30, 50],
        'E': [0, 180, 0, 0, 0, 180, 0, 0]
    }

    results = {}
    
    print("--- Chi-Square Calculation for Each Option ---")
    # Calculate chi-square for each option
    for option, observed in observed_data.items():
        total_offspring = sum(observed)
        
        # Calculate expected values
        expected_values = [(total_offspring * part) / total_ratio_parts for part in expected_ratio_parts]
        
        # Calculate chi-square statistic
        chi_square_value = sum([(o - e)**2 / e for o, e in zip(observed, expected_values) if e != 0])
        
        results[option] = chi_square_value
        print(f"Option {option}: Total Offspring = {total_offspring}, Chi-Square = {chi_square_value:.2f}")

    # Find the option with the highest chi-square value
    most_likely_rejection = max(results, key=results.get)
    max_chi_square = results[most_likely_rejection]
    
    # Degrees of freedom = number of categories - 1
    df = len(phenotypes) - 1
    # Critical value for alpha = 0.05 and df = 7
    critical_value = scipy.stats.chi2.ppf(0.95, df)
    
    print("\n--- Conclusion ---")
    print(f"The critical value for rejection at a significance level of 0.05 with {df} degrees of freedom is {critical_value:.2f}.")
    print(f"Option '{most_likely_rejection}' has the highest chi-square value ({max_chi_square:.2f}).")
    print("A higher chi-square value indicates a greater deviation from the expected outcome, making it the most likely to lead to the rejection of the null hypothesis of independent assortment.")

    print(f"\n--- Detailed Calculation for Winning Option '{most_likely_rejection}' ---")
    winning_observed = observed_data[most_likely_rejection]
    winning_total_offspring = sum(winning_observed)
    winning_expected = [(winning_total_offspring * part) / total_ratio_parts for part in expected_ratio_parts]

    equation_parts = []
    print("Chi-Square = Σ [ (Observed - Expected)² / Expected ]\n")
    for i in range(len(phenotypes)):
        obs = winning_observed[i]
        exp = winning_expected[i]
        component = (obs - exp)**2 / exp
        print(f"Phenotype: {phenotypes[i]}")
        print(f"   Observed (O) = {obs}, Expected (E) = {exp:.2f}")
        print(f"   Component = ({obs} - {exp:.2f})² / {exp:.2f} = {component:.2f}")
        equation_parts.append(f"({obs} - {exp:.2f})^2 / {exp:.2f}")
    
    print("\nFinal Equation:")
    print(f"χ² = {' + '.join(equation_parts)}")
    print(f"χ² = {max_chi_square:.2f}")

solve_chi_square()
<<<E>>>