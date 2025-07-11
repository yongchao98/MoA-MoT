import numpy as np

def solve_chi_square():
    """
    Calculates the chi-square statistic for each answer choice to determine
    which is most likely to lead to the rejection of the hypothesis of
    independent assortment.
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
    
    observed_data = {
        'A': [140, 10, 10, 10, 10, 10, 10, 100],
        'B': [180, 0, 0, 0, 60, 0, 0, 60],
        'C': [144, 45, 45, 16, 52, 16, 16, 16],
        'D': [150, 60, 50, 40, 30, 40, 30, 50],
        'E': [0, 180, 0, 0, 0, 180, 0, 0]
    }
    
    expected_ratio = np.array([27, 9, 9, 3, 9, 3, 3, 1])
    ratio_sum = np.sum(expected_ratio)
    
    results = {}
    
    print("Calculating Chi-Square values for each option...\n")

    for option, observed_counts in observed_data.items():
        observed_counts = np.array(observed_counts)
        total_observed = np.sum(observed_counts)
        
        expected_counts = (expected_ratio / ratio_sum) * total_observed
        
        # To avoid division by zero if an expected count is zero, which is not the case here.
        # We replace it with a very small number if it occurs.
        expected_counts[expected_counts == 0] = 1e-9

        chi_square_components = (observed_counts - expected_counts)**2 / expected_counts
        chi_square_value = np.sum(chi_square_components)
        
        results[option] = {'value': chi_square_value, 'components': chi_square_components}

        print(f"--- Option {option} ---")
        print(f"Observed: {observed_counts.tolist()}")
        print(f"Expected: {[round(e, 2) for e in expected_counts]}")
        print(f"Chi-Square Value: {chi_square_value:.2f}\n")

    # Find the option with the highest chi-square value
    most_likely_option = max(results, key=lambda k: results[k]['value'])
    max_chi_square_info = results[most_likely_option]

    print("="*40)
    print(f"Conclusion:")
    print(f"Option {most_likely_option} has the highest Chi-Square value ({max_chi_square_info['value']:.2f}).")
    print("This indicates the largest deviation from the expected 27:9:9:3:9:3:3:1 ratio.")
    print("Therefore, it is the most likely to lead to the rejection of the null hypothesis.\n")

    print(f"Final equation for Option {most_likely_option}:")
    equation_str = " + ".join([f"{comp:.2f}" for comp in max_chi_square_info['components']])
    print(f"Chi-Square = {equation_str}")
    print(f"Total = {max_chi_square_info['value']:.2f}")

solve_chi_square()
<<<E>>>