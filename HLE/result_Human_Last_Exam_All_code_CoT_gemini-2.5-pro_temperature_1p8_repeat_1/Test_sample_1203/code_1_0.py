import scipy.stats

def calculate_chi_square(option_label, observed_counts, phenotypes):
    """
    Calculates and prints the chi-square test results for a given set of observed counts.
    """
    # Expected phenotypic ratio for a trihybrid cross (TtRrYy x TtRrYy)
    expected_ratio = [27, 9, 9, 3, 9, 3, 3, 1]
    total_ratio_parts = sum(expected_ratio)

    # Calculate total observed offspring
    total_observed = sum(observed_counts)

    # Calculate expected counts
    expected_counts = [(part / total_ratio_parts) * total_observed for part in expected_ratio]

    # Calculate chi-square value
    chi_square_value = 0
    chi_square_components = []

    print(f"--- Analysis for Option {option_label} ---")
    print(f"Total Observed Offspring: {total_observed}\n")
    print("Chi-Square Calculation (Σ [(Observed - Expected)² / Expected]):\n")

    for i in range(len(observed_counts)):
        obs = observed_counts[i]
        exp = expected_counts[i]
        
        # Avoid division by zero if expected is 0, although not the case here.
        if exp == 0:
            component = float('inf') if obs > 0 else 0
        else:
            component = (obs - exp)**2 / exp
        
        chi_square_components.append(component)
        print(f"Phenotype: {phenotypes[i]}")
        print(f"Equation: ({obs} - {exp:.2f})**2 / {exp:.2f} = {component:.2f}")

    chi_square_value = sum(chi_square_components)

    # Degrees of freedom = number of categories - 1
    df = len(phenotypes) - 1
    
    # Significance level
    alpha = 0.05
    critical_value = scipy.stats.chi2.ppf(1 - alpha, df)

    print("\n----------------------------------")
    print(f"Total Chi-Square Value (χ²) for Option {option_label}: {chi_square_value:.2f}")
    print(f"Degrees of Freedom: {df}")
    print(f"Critical Value at α=0.05: {critical_value:.2f}")
    
    if chi_square_value > critical_value:
        print("Result: Reject the null hypothesis of independent assortment.")
    else:
        print("Result: Fail to reject the null hypothesis of independent assortment.")
    print("==================================\n")
    
    return chi_square_value, option_label

def main():
    phenotypes = [
        "Tall, round, yellow", "Tall, round, green", "Tall, wrinkled, yellow",
        "Tall, wrinkled, green", "Dwarf, round, yellow", "Dwarf, round, green",
        "Dwarf, wrinkled, yellow", "Dwarf, wrinkled, green"
    ]

    # Observed counts from the answer choices
    observed_data = {
        "A": [140, 10, 10, 10, 10, 10, 10, 100],
        "B": [180, 0, 0, 0, 60, 0, 0, 60],
        "C": [144, 45, 45, 16, 52, 16, 16, 16],
        "D": [150, 60, 50, 40, 30, 40, 30, 50],
        "E": [0, 180, 0, 0, 0, 180, 0, 0]
    }

    results = []
    for option_label, counts in observed_data.items():
        chi_square, label = calculate_chi_square(option_label, counts, phenotypes)
        results.append((chi_square, label))

    # Find the option with the highest chi-square value
    most_likely_rejection = max(results, key=lambda item: item[0])
    
    print("\n--- Final Conclusion ---")
    print(f"Option {most_likely_rejection[1]} has the highest chi-square value ({most_likely_rejection[0]:.2f}).")
    print("A higher chi-square value indicates a greater deviation from the expected outcome under independent assortment.")
    print(f"Therefore, the results in Option {most_likely_rejection[1]} would most likely lead to the rejection of the hypothesis.")

if __name__ == "__main__":
    main()
