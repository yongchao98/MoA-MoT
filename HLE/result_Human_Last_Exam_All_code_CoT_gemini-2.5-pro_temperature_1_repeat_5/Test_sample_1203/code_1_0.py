import numpy as np

def solve_chi_square():
    """
    Calculates the chi-square value for each answer choice to determine which is most likely
    to lead to the rejection of the null hypothesis of independent assortment.
    """
    # Phenotype order: Tall/Round/Yellow, Tall/Round/green, Tall/wrinkled/Yellow, Tall/wrinkled/green,
    #                  dwarf/Round/Yellow, dwarf/Round/green, dwarf/wrinkled/Yellow, dwarf/wrinkled/green
    observed_counts = {
        'A': np.array([140, 10, 10, 10, 10, 10, 10, 100]),
        'B': np.array([180, 0, 0, 0, 60, 0, 0, 60]),
        'C': np.array([144, 45, 45, 16, 52, 16, 16, 16]),
        'D': np.array([150, 60, 50, 40, 30, 40, 30, 50]),
        'E': np.array([0, 180, 0, 0, 0, 180, 0, 0])
    }

    # Expected phenotypic ratio for a trihybrid cross (Tt Rr Yy x Tt Rr Yy)
    expected_ratio = np.array([27, 9, 9, 3, 9, 3, 3, 1])
    ratio_total = expected_ratio.sum()

    results = {}
    
    print("Calculating Chi-Square values for each option:\n")

    for option, observed in observed_counts.items():
        total_offspring = observed.sum()
        expected_counts = (expected_ratio / ratio_total) * total_offspring
        
        # Chi-square calculation: sum(((O-E)^2)/E)
        # Add a small epsilon to expected_counts to avoid division by zero if an expected count is 0
        chi_square_components = (observed - expected_counts)**2 / (expected_counts + 1e-9)
        chi_square_value = chi_square_components.sum()
        
        results[option] = {
            'chi_square': chi_square_value,
            'observed': observed,
            'expected': expected_counts
        }
        print(f"Option {option}: χ² = {chi_square_value:.2f}")

    # Find the option with the maximum chi-square value
    best_option = max(results, key=lambda k: results[k]['chi_square'])
    
    print(f"\nOption {best_option} has the highest chi-square value, indicating the largest deviation from the expected ratio.")
    print("This option would most likely lead to the rejection of the null hypothesis of independent assortment.")

    # Display the full calculation for the winning option
    print(f"\nHere is the detailed calculation for Option {best_option}:")
    
    winner_data = results[best_option]
    obs = winner_data['observed']
    exp = winner_data['expected']
    
    equation_parts = []
    for i in range(len(obs)):
        equation_parts.append(f"({obs[i]} - {exp[i]:.2f})^2 / {exp[i]:.2f}")
    
    full_equation = " + ".join(equation_parts)
    print(f"\nχ² = {full_equation}")
    print(f"   = {winner_data['chi_square']:.2f}")
    
    # The degrees of freedom is number of categories - 1 = 8 - 1 = 7.
    # The critical value for p=0.05 and df=7 is 14.07.
    df = len(observed_counts[best_option]) - 1
    critical_value = 14.07
    print(f"\nWith {df} degrees of freedom, the critical value at a significance level of 0.05 is {critical_value}.")
    print(f"Since {winner_data['chi_square']:.2f} > {critical_value}, the null hypothesis would be rejected.")
    
    # Return the final letter answer
    print(f"\n<<<E>>>")

solve_chi_square()