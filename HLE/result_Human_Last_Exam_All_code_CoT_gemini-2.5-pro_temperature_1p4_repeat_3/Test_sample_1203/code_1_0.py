import numpy as np

def solve_chi_square():
    """
    Calculates the Chi-Square value for each set of offspring phenotypes to determine
    which is most likely to lead to the rejection of the hypothesis of independent assortment.
    """
    # Phenotype labels in standard order
    phenotype_labels = [
        "tall, round, yellow", "tall, round, green", "tall, wrinkled, yellow",
        "tall, wrinkled, green", "dwarf, round, yellow", "dwarf, round, green",
        "dwarf, wrinkled, yellow", "dwarf, wrinkled, green"
    ]
    
    # Mendelian ratio for a trihybrid cross
    ratio_parts = np.array([27, 9, 9, 3, 9, 3, 3, 1])
    
    # Observed data from the answer choices
    observed_data = {
        'A': np.array([140, 10, 10, 10, 10, 10, 10, 100]),
        'B': np.array([180, 0, 0, 0, 60, 0, 0, 60]),
        'C': np.array([144, 45, 45, 16, 52, 16, 16, 16]),
        'D': np.array([150, 60, 50, 40, 30, 40, 30, 50]),
        'E': np.array([0, 180, 0, 0, 0, 180, 0, 0])
    }

    results = {}

    for option, observed in observed_data.items():
        total_observed = np.sum(observed)
        
        # Calculate expected counts
        expected = (ratio_parts / 64.0) * total_observed
        
        # Calculate Chi-Square value
        # To avoid division by zero, we replace expected values of 0 with a very small number.
        # However, in this Mendelian ratio, no expected count will be zero unless the total is zero.
        chi_square_terms = (observed - expected)**2 / expected
        chi_square_value = np.sum(chi_square_terms)
        results[option] = chi_square_value

        print(f"--- Option {option} ---")
        print(f"Observed counts: {list(observed)}")
        print(f"Total Offspring: {total_observed}")
        print(f"Expected counts: {[round(e, 2) for e in expected]}")
        
        # Print the detailed Chi-Square calculation
        calculation_str = "Chi-Square = "
        terms_str = []
        for i in range(len(observed)):
            terms_str.append(f"({observed[i]} - {expected[i]:.2f})^2 / {expected[i]:.2f}")
        calculation_str += " + ".join(terms_str)
        calculation_str += f" = {chi_square_value:.2f}"
        print(calculation_str)
        print("-" * 20 + "\n")

    # Find the option with the maximum Chi-Square value
    most_likely_rejection = max(results, key=results.get)
    
    print("\n--- Conclusion ---")
    print("The Chi-Square test measures the discrepancy between observed and expected frequencies.")
    print("A higher Chi-Square value indicates a larger discrepancy and a higher likelihood of rejecting the null hypothesis (of independent assortment).")
    print(f"\nComparing the results:")
    for option, value in results.items():
        print(f"Option {option}: Chi-Square = {value:.2f}")

    print(f"\nOption '{most_likely_rejection}' has the highest Chi-Square value, making it the combination of offspring phenotypes that would most likely lead to the rejection of the hypothesis of independent assortment.")
    
    print(f"\n<<<E>>>")

solve_chi_square()