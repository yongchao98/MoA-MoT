import numpy as np

def solve_chi_square():
    """
    Calculates the chi-square value for each answer choice to determine which
    is most likely to lead to the rejection of the hypothesis of independent assortment.
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

    observed_counts = {
        'A': [140, 10, 10, 10, 10, 10, 10, 100],
        'B': [180, 0, 0, 0, 60, 0, 0, 60],
        'C': [144, 45, 45, 16, 52, 16, 16, 16],
        'D': [150, 60, 50, 40, 30, 40, 30, 50],
        'E': [0, 180, 0, 0, 0, 180, 0, 0]
    }

    # Expected ratio for a trihybrid cross (Tt Rr Yy x Tt Rr Yy)
    expected_ratio = np.array([27, 9, 9, 3, 9, 3, 3, 1])
    total_ratio_parts = np.sum(expected_ratio)

    # Degrees of freedom for this test is (number of categories - 1) = 8 - 1 = 7.
    # The critical value for alpha = 0.05 and df = 7 is approximately 14.07.
    # A calculated chi-square value greater than 14.07 would lead to rejection of the null hypothesis.

    print("Calculating the Chi-Square value for each option...\n")

    chi_square_results = {}

    for option, observed in observed_counts.items():
        observed = np.array(observed)
        total_observed = np.sum(observed)

        # Calculate expected counts
        expected = (expected_ratio / total_ratio_parts) * total_observed

        # Calculate chi-square components
        # Handle cases where expected count is zero to avoid division by zero
        # Although not an issue with this specific ratio
        with np.errstate(divide='ignore', invalid='ignore'):
            chi_square_components = (observed - expected)**2 / expected
            chi_square_components[np.isnan(chi_square_components)] = 0 # Replace NaN with 0 if E=0

        chi_square_value = np.sum(chi_square_components)
        chi_square_results[option] = chi_square_value

        print(f"--- Option {option} ---")
        print(f"Observed counts: {observed.tolist()}")
        print(f"Total Observed: {total_observed}")
        print(f"Expected counts: {[f'{val:.2f}' for val in expected]}")
        
        # Build and print the equation string
        equation_parts = []
        for o, e in zip(observed, expected):
            equation_parts.append(f"(({o} - {e:.2f})^2 / {e:.2f})")
        
        component_values_str = " + ".join([f"{comp:.2f}" for comp in chi_square_components])
        
        print(f"Chi-Square Calculation:\n  = {component_values_str}")
        print(f"Total Chi-Square for Option {option}: {chi_square_value:.2f}\n")

    # Find the option with the highest chi-square value
    most_likely_rejection = max(chi_square_results, key=chi_square_results.get)
    max_chi_square_value = chi_square_results[most_likely_rejection]
    
    print("--- Conclusion ---")
    print("The null hypothesis of independent assortment is rejected when the calculated Chi-Square value is high, as this indicates a large deviation from the expected results.")
    print(f"Comparing the results, Option {most_likely_rejection} has the highest Chi-Square value of {max_chi_square_value:.2f}.")
    print("Therefore, this combination of offspring phenotypes would most likely lead to the rejection of the hypothesis.")
    print(f"<<<{most_likely_rejection}>>>")

solve_chi_square()