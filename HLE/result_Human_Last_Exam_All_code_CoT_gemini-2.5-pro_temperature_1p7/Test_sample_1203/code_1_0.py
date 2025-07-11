import scipy.stats

def solve_chi_square_problem():
    """
    Calculates the chi-square statistic for five different sets of observed genetic cross data
    to determine which is most likely to lead to the rejection of the null hypothesis of
    independent assortment.
    """
    # Null hypothesis: Genes assort independently.
    # The expected phenotypic ratio for a trihybrid cross (Tt Rr Yy x Tt Rr Yy) is 27:9:9:3:9:3:3:1.
    expected_ratio = [27, 9, 9, 3, 9, 3, 3, 1]
    phenotypes = [
        "Tall, round, yellow", "Tall, round, green", "Tall, wrinkled, yellow",
        "Tall, wrinkled, green", "Dwarf, round, yellow", "Dwarf, round, green",
        "Dwarf, wrinkled, yellow", "Dwarf, wrinkled, green"
    ]

    # Observed counts from the answer choices
    observed_data = {
        'A': [140, 10, 10, 10, 10, 10, 10, 100],
        'B': [180, 0, 0, 0, 60, 0, 0, 60],
        'C': [144, 45, 45, 16, 52, 16, 16, 16],
        'D': [150, 60, 50, 40, 30, 40, 30, 50],
        'E': [0, 180, 0, 0, 0, 180, 0, 0]
    }

    # Chi-square test parameters
    alpha = 0.05
    df = len(phenotypes) - 1  # Degrees of freedom = (number of categories - 1)
    critical_value = scipy.stats.chi2.ppf(1 - alpha, df)

    print("--- Chi-Square Test for Independent Assortment ---")
    print(f"The null hypothesis is that the genes assort independently, giving an expected phenotypic ratio of {':'.join(map(str, expected_ratio))}.")
    print(f"With {df} degrees of freedom and a significance level of {alpha}, the critical Chi-Square value is {critical_value:.2f}.")
    print("A calculated Chi-Square value greater than this critical value would lead to rejection of the null hypothesis.")
    print("The option most likely to lead to rejection will have the HIGHEST calculated Chi-Square value.\n")

    results = {}
    for option, observed in observed_data.items():
        total_observed = sum(observed)
        total_ratio = sum(expected_ratio)
        expected_counts = [(r / total_ratio) * total_observed for r in expected_ratio]
        
        chi_square_stat = 0
        for o, e in zip(observed, expected_counts):
            if e > 0:
                chi_square_stat += ((o - e)**2) / e
        
        results[option] = chi_square_stat
        print(f"--- Option {option} ---")
        print(f"Calculated Chi-Square (χ²): {chi_square_stat:.2f}")

    # Find the option with the maximum Chi-Square value
    most_likely_rejection_option = max(results, key=results.get)
    highest_chi_square = results[most_likely_rejection_option]

    print("\n--- Conclusion ---")
    print(f"Option '{most_likely_rejection_option}' has the highest Chi-Square value ({highest_chi_square:.2f}), indicating the largest deviation from the expected ratio.")
    print("Therefore, this combination would most likely lead to rejection of the hypothesis of independent assortment.\n")
    
    print(f"Detailed Calculation for Final Answer (Option {most_likely_rejection_option}):")
    # Recalculate and print the step-by-step equation for the final answer
    obs = observed_data[most_likely_rejection_option]
    total_observed = sum(obs)
    total_ratio = sum(expected_ratio)
    exp = [(r / total_ratio) * total_observed for r in expected_ratio]
    
    chi_components_str = []
    for o, e in zip(obs, exp):
        component_str = f"({o} - {e:.2f})² / {e:.2f}"
        chi_components_str.append(component_str)
        
    print(f"χ² = {' + '.join(chi_components_str)}")

    chi_components_val = []
    for o, e in zip(obs, exp):
        val = ((o - e)**2) / e if e > 0 else 0
        chi_components_val.append(f"{val:.2f}")

    print(f"χ² = {' + '.join(chi_components_val)} = {highest_chi_square:.2f}")


solve_chi_square_problem()