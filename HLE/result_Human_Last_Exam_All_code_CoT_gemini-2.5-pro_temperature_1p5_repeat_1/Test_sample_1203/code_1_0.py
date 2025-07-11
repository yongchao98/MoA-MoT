import numpy as np

def solve_chi_square():
    """
    Calculates the Chi-Square value for each set of observed offspring counts
    and identifies the one most likely to reject the null hypothesis of
    independent assortment.
    """
    # The expected ratio for the 8 phenotypes in the order they appear in the problem
    expected_ratio = np.array([27, 9, 9, 3, 9, 3, 3, 1])
    
    # Observed counts for each answer choice
    observed_counts = {
        'A': np.array([140, 10, 10, 10, 10, 10, 10, 100]),
        'B': np.array([180, 0, 0, 0, 60, 0, 0, 60]),
        'C': np.array([144, 45, 45, 16, 52, 16, 16, 16]),
        'D': np.array([150, 60, 50, 40, 30, 40, 30, 50]),
        'E': np.array([0, 180, 0, 0, 0, 180, 0, 0]),
    }

    chi_square_results = {}

    print("Calculating Chi-Square values for each option:")
    print("-" * 50)
    for option, observed in observed_counts.items():
        total_offspring = np.sum(observed)
        
        # Calculate expected counts based on the total offspring and the expected ratio
        expected = (expected_ratio / 64) * total_offspring
        
        # Calculate the Chi-Square value
        # Handle cases where expected count is 0 to avoid division by zero, though not applicable here
        # Add a small epsilon where observed is 0 and expected is 0
        chi_square_components = (observed - expected)**2 / (expected + 1e-9)
        chi_square_value = np.sum(chi_square_components)
        
        chi_square_results[option] = chi_square_value
        print(f"Option {option}: Total Offspring = {total_offspring}, Chi-Square = {chi_square_value:.2f}")

    # Find the option with the maximum Chi-Square value
    max_chi_option = max(chi_square_results, key=chi_square_results.get)
    max_chi_value = chi_square_results[max_chi_option]

    print("-" * 50)
    print(f"\nOption {max_chi_option} has the highest Chi-Square value ({max_chi_value:.2f}), indicating the largest deviation from the expected ratio. This would most likely lead to the rejection of the null hypothesis.")

    # Display the full calculation for the winning option
    print("\nThe Chi-Square calculation for this option is:")
    
    obs_winner = observed_counts[max_chi_option]
    total_winner = np.sum(obs_winner)
    exp_winner = (expected_ratio / 64) * total_winner
    
    equation_parts = []
    for o, e in zip(obs_winner, exp_winner):
        equation_parts.append(f"({o} - {e:.2f})^2 / {e:.2f}")
    
    equation_string = " + ".join(equation_parts)
    print(f"\nChi-Square = {equation_string} = {max_chi_value:.2f}")


solve_chi_square()