import scipy.stats

def solve_chi_square_problem():
    """
    This script calculates the Chi-Square (χ²) value for each answer choice
    to determine which set of offspring phenotypes would most likely lead to
    the rejection of the hypothesis of independent assortment.
    """
    
    # --- Data ---
    phenotypes = [
        "Tall, round, yellow", "Tall, round, green", "Tall, wrinkled, yellow",
        "Tall, wrinkled, green", "dwarf, round, yellow", "dwarf, round, green",
        "dwarf, wrinkled, yellow", "dwarf, wrinkled, green",
    ]

    # Observed counts for each choice
    observed_data = {
        'A': [140, 10, 10, 10, 10, 10, 10, 100],
        'B': [180, 0, 0, 0, 60, 0, 0, 60],
        'C': [144, 45, 45, 16, 52, 16, 16, 16],
        'D': [150, 60, 50, 40, 30, 40, 30, 50],
        'E': [0, 180, 0, 0, 0, 180, 0, 0]
    }

    # Expected ratio parts for a trihybrid cross (27:9:9:3:9:3:3:1)
    expected_ratio_parts = [27, 9, 9, 3, 9, 3, 3, 1]
    total_ratio_parts = sum(expected_ratio_parts)

    # --- Calculations ---
    results = {}

    print("--- Chi-Square (χ²) Calculation for Each Option ---")
    print("The null hypothesis of independent assortment predicts a 27:9:9:3:9:3:3:1 phenotypic ratio.\n")
    
    for choice, observed in observed_data.items():
        total_n = sum(observed)
        
        # Calculate expected counts based on the total N
        expected = [(total_n * part) / total_ratio_parts for part in expected_ratio_parts]
        
        chi_square_terms = []
        equation_parts = []
        term_values_str_parts = []
        
        for o, e in zip(observed, expected):
            # The contribution to the Chi-Square sum is (Observed - Expected)^2 / Expected
            term = ((o - e) ** 2) / e
            chi_square_terms.append(term)
            
            # For printing the final equation with all its numbers
            equation_parts.append(f"({o} - {e:.2f})^2/{e:.2f}")
            term_values_str_parts.append(f"{term:.2f}")

        chi_square_val = sum(chi_square_terms)
        results[choice] = chi_square_val
        
        print("-" * 50)
        print(f"Analysis for Choice {choice}:")
        print(f"Observed Counts: {observed}")
        print(f"Total Offspring (N): {total_n}")
        print(f"Expected Counts: {[round(e, 2) for e in expected]}")
        
        # Print the Chi-Square equation with numbers
        print(f"\nχ² calculation:")
        print(" + ".join(equation_parts))
        print("= " + " + ".join(term_values_str_parts))
        print(f"= {chi_square_val:.2f}\n")

    # --- Conclusion ---
    # Find the choice with the maximum chi_square value
    max_choice = max(results, key=results.get)
    max_value = results[max_choice]

    # Determine the critical value for context
    df = len(phenotypes) - 1
    alpha = 0.05
    critical_value = scipy.stats.chi2.ppf(1 - alpha, df)

    print("-" * 50)
    print("--- Final Analysis ---")
    print(f"The degrees of freedom (df) is {len(phenotypes)} categories - 1 = {df}.")
    print(f"At a significance level (α) of {alpha}, the critical χ² value is {critical_value:.2f}.")
    print("A calculated χ² value greater than this critical value leads to rejection of the null hypothesis.\n")

    print("Summary of Chi-Square values:")
    for choice, value in sorted(results.items(), key=lambda item: item[1], reverse=True):
        print(f"  Choice {choice}: χ² = {value:.2f}")

    print(f"\nChoice {max_choice} has the highest Chi-Square value ({max_value:.2f}). This indicates the largest deviation from the expected Mendelian ratio.")
    print(f"\nTherefore, the data in Choice {max_choice} would most likely lead to the rejection of the hypothesis of independent assortment.")

if __name__ == "__main__":
    solve_chi_square_problem()