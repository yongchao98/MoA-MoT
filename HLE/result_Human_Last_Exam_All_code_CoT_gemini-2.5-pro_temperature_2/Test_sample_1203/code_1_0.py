import numpy as np

def calculate_chi_square():
    """
    Calculates the chi-square value for each set of observed offspring counts from a trihybrid cross.
    """
    # Phenotype order:
    # 1. tall, round, yellow
    # 2. tall, round, green
    # 3. tall, wrinkled, yellow
    # 4. tall, wrinkled, green
    # 5. dwarf, round, yellow
    # 6. dwarf, round, green
    # 7. dwarf, wrinkled, yellow
    # 8. dwarf, wrinkled, green

    observed_data = {
        'A': [140, 10, 10, 10, 10, 10, 10, 100],
        'B': [180, 0, 0, 0, 60, 0, 0, 60],
        'C': [144, 45, 45, 16, 52, 16, 16, 16],
        'D': [150, 60, 50, 40, 30, 40, 30, 50],
        'E': [0, 180, 0, 0, 0, 180, 0, 0]
    }

    # Expected ratio for a trihybrid cross (TtRrYy x TtRrYy)
    expected_ratio = np.array([27, 9, 9, 3, 9, 3, 3, 1])
    ratio_sum = np.sum(expected_ratio) # Should be 64

    results = {}

    for option, observed_counts in observed_data.items():
        observed_counts = np.array(observed_counts)
        total_offspring = np.sum(observed_counts)

        # Calculate expected counts
        expected_counts = (expected_ratio / ratio_sum) * total_offspring

        # Calculate chi-square value and build the equation string
        chi2_value = 0
        equation_terms = []
        for i in range(len(observed_counts)):
            obs = observed_counts[i]
            exp = expected_counts[i]
            # Avoid division by zero, though not expected here
            if exp == 0:
                # If exp is 0 and obs is also 0, the contribution is 0.
                # If exp is 0 and obs is not, chi-square is theoretically infinite.
                term_val = 0 if obs == 0 else float('inf')
            else:
                term_val = ((obs - exp)**2) / exp
            
            chi2_value += term_val
            # Add each number in the final equation
            equation_terms.append(f"({obs} - {exp:.2f})^2 / {exp:.2f}")

        # Store the final calculated value
        results[option] = chi2_value

        # Print the detailed calculation for the current option
        print(f"--- Choice {option} ---")
        print(f"Total Offspring: {total_offspring}")
        print(f"Observed Counts: {list(observed_counts)}")
        print(f"Expected Counts: {[float(f'{c:.2f}') for c in expected_counts]}")
        print(f"Chi-Square Calculation:")
        print(" + ".join(equation_terms))
        print(f"Chi-Square (χ²) Value for Option {option}: {chi2_value:.2f}\n")

    # Determine which option has the highest chi-square value
    most_likely_rejection = max(results, key=results.get)
    print("="*30)
    print(f"The highest chi-square value is from option {most_likely_rejection}.")
    print("A higher chi-square value indicates a greater deviation from the expected Mendelian ratio.")
    print("Therefore, this option is the most likely to lead to a rejection of the hypothesis of independent assortment.")
    
    # Degrees of freedom for this test is (number of categories - 1) = 8 - 1 = 7.
    # The critical value for α = 0.05 and df = 7 is approximately 14.07.
    # All options except C have chi-square values far exceeding this critical value.
    # Option E has the largest chi-square value by a significant margin.

if __name__ == '__main__':
    calculate_chi_square()
<<<E>>>