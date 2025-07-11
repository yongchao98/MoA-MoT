import scipy.stats

def calculate_and_print_chi_square(option_name, observed_counts, pheno_labels):
    """
    Calculates the Chi-Square value for a given set of observed counts against the
    expected 27:9:9:3:9:3:3:1 ratio.
    """
    expected_ratio_parts = [27, 9, 9, 3, 9, 3, 3, 1]
    total_ratio_parts = sum(expected_ratio_parts)  # This is 64

    # Calculate the total number of observed offspring
    total_offspring = sum(observed_counts)

    # Calculate expected counts based on the total offspring and the ratio
    expected_counts = [(total_offspring * part) / total_ratio_parts for part in expected_ratio_parts]

    chi_square_sum = 0
    equation_terms_str = []
    
    print(f"--- Calculation for Option {option_name} ---")
    print(f"Total Observed Offspring (N): {total_offspring}")
    print("Equation: Chi-Square = sum[ (Observed - Expected)^2 / Expected ]")
    print("\nTerms of the sum:")
    
    # Calculate each term of the Chi-Square sum
    for i in range(len(observed_counts)):
        O = observed_counts[i]
        E = expected_counts[i]
        
        # Avoid division by zero if an expected count is zero
        term = (O - E)**2 / E if E > 0 else 0
        chi_square_sum += term
        
        # Print each number and calculation that goes into the final equation
        print(f"  {pheno_labels[i]:<25}: ({O:.2f} - {E:.2f})^2 / {E:.2f} = {term:.2f}")
        equation_terms_str.append(f"{term:.2f}")

    print(f"\nFinal Equation for Option {option_name}:")
    print(f"Chi-Square = {' + '.join(equation_terms_str)}")
    print(f"Total Chi-Square Value for Option {option_name}: {chi_square_sum:.2f}")
    print("-" * 40)
    
    return chi_square_sum

def main():
    """
    Main function to run the analysis for all answer choices.
    """
    phenotype_labels = [
        "Tall, round, yellow", "Tall, round, green", "Tall, wrinkled, yellow",
        "Tall, wrinkled, green", "Dwarf, round, yellow", "Dwarf, round, green",
        "Dwarf, wrinkled, yellow", "Dwarf, wrinkled, green"
    ]

    # Observed data from the answer choices
    observed_data = {
        'A': [140, 10, 10, 10, 10, 10, 10, 100],
        'B': [180, 0, 0, 0, 60, 0, 0, 60],
        'C': [144, 45, 45, 16, 52, 16, 16, 16],
        'D': [150, 60, 50, 40, 30, 40, 30, 50],
        'E': [0, 180, 0, 0, 0, 180, 0, 0]
    }

    # Degrees of freedom = number of phenotypic categories - 1
    df = len(phenotype_labels) - 1
    alpha = 0.05
    # Find the critical value for the given significance level
    critical_value = scipy.stats.chi2.ppf(1 - alpha, df)

    print("The null hypothesis of independent assortment predicts a 27:9:9:3:9:3:3:1 phenotypic ratio.")
    print(f"For this test (df={df}, alpha={alpha}), the Chi-Square critical value is {critical_value:.3f}.")
    print("The null hypothesis is rejected if the calculated Chi-Square > critical value.")
    print("The option MOST LIKELY to be rejected will have the HIGHEST Chi-Square value.\n")

    results = {}
    for option, counts in observed_data.items():
        results[option] = calculate_and_print_chi_square(option, counts, phenotype_labels)

    # Find the option with the maximum Chi-Square value
    max_chi_square = -1
    best_option = None
    for option, chi_val in results.items():
        if chi_val > max_chi_square:
            max_chi_square = chi_val
            best_option = option

    print("\n--- SUMMARY ---")
    for option, chi_val in sorted(results.items()):
        print(f"Option {option}: Chi-Square = {chi_val:.2f}")
    
    print(f"\nConclusion: Option {best_option} has the highest Chi-Square value ({max_chi_square:.2f}).")
    print("This indicates the greatest deviation from the expected ratio, making it the most likely")
    print("to lead to a rejection of the hypothesis of independent assortment.")


if __name__ == '__main__':
    main()
<<<E>>>