import scipy.stats

def solve_chi_square_problem():
    """
    Calculates the chi-square value for each option to determine which is most likely
    to lead to the rejection of the hypothesis of independent assortment.
    """
    # Phenotype order corresponds to the 27:9:9:3:9:3:3:1 ratio
    phenotypes = [
        "Tall, round, yellow", "Tall, round, green", "Tall, wrinkled, yellow",
        "Tall, wrinkled, green", "Dwarf, round, yellow", "Dwarf, round, green",
        "Dwarf, wrinkled, yellow", "Dwarf, wrinkled, green"
    ]

    # Observed counts for each answer choice
    observed_counts = {
        'A': [140, 10, 10, 10, 10, 10, 10, 100],
        'B': [180, 0, 0, 0, 60, 0, 0, 60],
        'C': [144, 45, 45, 16, 52, 16, 16, 16],
        'D': [150, 60, 50, 40, 30, 40, 30, 50],
        'E': [0, 180, 0, 0, 0, 180, 0, 0]
    }

    # Expected ratio parts for a trihybrid cross
    expected_ratio = [27, 9, 9, 3, 9, 3, 3, 1]
    total_ratio_parts = sum(expected_ratio)

    chi_square_results = {}

    # Calculate chi-square for each option
    for option, observed in observed_counts.items():
        total_observed = sum(observed)
        expected = [(part / total_ratio_parts) * total_observed for part in expected_ratio]
        
        chi_square = sum(((o - e)**2) / e for o, e in zip(observed, expected) if e > 0)
        chi_square_results[option] = chi_square

    # Find the option with the highest chi-square value
    max_chi_square_option = max(chi_square_results, key=chi_square_results.get)

    print("The null hypothesis of independent assortment predicts a phenotypic ratio of 27:9:9:3:9:3:3:1.")
    print("A higher chi-square (χ²) value indicates a greater deviation from the expected ratio.\n")

    for option, value in chi_square_results.items():
        print(f"Chi-square for Option {option}: {value:.2f}")

    print(f"\nOption {max_chi_square_option} has the highest chi-square value ({chi_square_results[max_chi_square_option]:.2f}).")
    print("This indicates the greatest deviation and is therefore the most likely to lead to rejection of the hypothesis.\n")

    # Show the detailed calculation for the winning option
    print(f"Detailed calculation for Option {max_chi_square_option}:")
    winning_observed = observed_counts[max_chi_square_option]
    total_observed = sum(winning_observed)
    winning_expected = [(part / total_ratio_parts) * total_observed for part in expected_ratio]

    equation_parts = []
    for i in range(len(winning_observed)):
        o = winning_observed[i]
        e = winning_expected[i]
        part_str = f"({o} - {e:.2f})^2 / {e:.2f}"
        equation_parts.append(part_str)

    full_equation = " + ".join(equation_parts)
    print(f"χ² = {full_equation}")
    print(f"χ² = {chi_square_results[max_chi_square_option]:.2f}")
    
    # Compare to critical value
    df = len(expected_ratio) - 1
    p_value = 0.05
    critical_value = scipy.stats.chi2.ppf(1 - p_value, df)
    print(f"\nFor a significance level of {p_value} with {df} degrees of freedom, the critical χ² value is {critical_value:.2f}.")
    print(f"Since {chi_square_results[max_chi_square_option]:.2f} is much larger than {critical_value:.2f}, the null hypothesis would be rejected.")

solve_chi_square_problem()
<<<E>>>