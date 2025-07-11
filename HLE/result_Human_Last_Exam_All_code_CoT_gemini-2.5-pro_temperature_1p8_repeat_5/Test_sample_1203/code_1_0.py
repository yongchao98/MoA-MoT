import numpy as np

def solve_chi_square():
    """
    Calculates the chi-square value for each answer choice to find which
    one would most likely lead to the rejection of the hypothesis of independent assortment.
    """
    # Phenotype order: T_R_Y_, T_R_yy, T_rrY_, T_rryy, ttR_Y_, ttR_yy, ttrrY_, ttrryy
    # Observed data for each choice
    observed_data = {
        'A': [140, 10, 10, 10, 10, 10, 10, 100],
        'B': [180, 0, 0, 0, 60, 0, 0, 60],
        'C': [144, 45, 45, 16, 52, 16, 16, 16],
        'D': [150, 60, 50, 40, 30, 40, 30, 50],
        'E': [0, 180, 0, 0, 0, 180, 0, 0]
    }

    # Expected Mendelian ratio parts for a trihybrid cross
    expected_ratio = np.array([27, 9, 9, 3, 9, 3, 3, 1])

    # Dictionary to store chi-square results
    chi_square_results = {}
    calculation_details = {}

    print("--- Chi-Square (χ²) Calculation Summary ---")
    
    # Loop through each answer choice
    for choice, observed_counts in observed_data.items():
        observed = np.array(observed_counts)
        total_offspring = np.sum(observed)
        
        # Calculate expected counts based on the total offspring and the expected ratio
        total_ratio_parts = np.sum(expected_ratio)
        expected = total_offspring * (expected_ratio / total_ratio_parts)
        
        # Calculate individual chi-square terms and the total value
        # Avoid division by zero (not an issue here as all expected values are > 0)
        chi_square_terms = ((observed - expected)**2) / expected
        chi_square_value = np.sum(chi_square_terms)
        
        chi_square_results[choice] = chi_square_value
        calculation_details[choice] = {
            "terms": chi_square_terms,
            "total": chi_square_value
        }
        
        print(f"Choice {choice}: χ² = {chi_square_value:.2f}")

    # Find the choice with the maximum Chi-square value
    max_chi_choice = max(chi_square_results, key=chi_square_results.get)
    max_chi_details = calculation_details[max_chi_choice]

    print("\n--- Detailed Calculation for the Most Likely Rejection ---")
    print(f"Choice {max_chi_choice} has the highest Chi-square value, indicating the largest deviation from the expected ratio.")
    print("This choice would most likely lead to the rejection of the hypothesis of independent assortment.")
    
    # Format the equation string showing each term being added
    summation_string = " + ".join([f"{term:.2f}" for term in max_chi_details["terms"]])
    
    print("\nThe final equation for this choice is:")
    print(f"χ² = {summation_string} = {max_chi_details['total']:.2f}")

# Execute the function
solve_chi_square()