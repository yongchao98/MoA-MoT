import scipy.stats
import numpy as np
import math

def find_minimum_observations():
    """
    Calculates the minimum number of observations per group to achieve
    statistical significance with a Bonferroni-corrected p-value.
    """
    original_alpha = 0.05
    num_tests = 5
    corrected_alpha = original_alpha / num_tests

    # Start with the initial sample size and iterate upwards
    n = 5
    
    while True:
        # Generate data for Group 1 by repeating the pattern [1, 2, 3, 4, 5]
        base_g1 = [1, 2, 3, 4, 5]
        group1 = np.tile(base_g1, n // len(base_g1)).tolist() + base_g1[:n % len(base_g1)]

        # Generate data for Group 2 with ~60% 3s and ~40% 2s
        num_threes = round(0.6 * n)
        num_twos = n - num_threes
        group2 = [3] * num_threes + [2] * num_twos

        # Perform the two-sided Mann-Whitney U test
        u_statistic, p_value = scipy.stats.mannwhitneyu(group1, group2, alternative='two-sided')

        # Check if the p-value is below the corrected alpha
        if p_value < corrected_alpha:
            print("Significance found.")
            print(f"The Bonferroni-corrected p-value threshold is: {corrected_alpha}")
            print(f"The minimum number of observations per group required is: {n}")
            print("\n--- Final Statistical Test Details ---")
            print(f"Sample Size (n) per group: {n}")
            print(f"Mann-Whitney U-statistic: {u_statistic}")
            print(f"p-value: {p_value}")
            
            # This is the final numerical answer
            return n

        # Increment n and continue the search
        n += 1
        
        # Safety break to prevent an accidental infinite loop
        if n > 500:
            print("Could not find a solution for n up to 500.")
            return None

# Run the function and store the result
min_n = find_minimum_observations()

# The final answer is formatted as requested by the user.
print(f"\n<<<>>>")
print(f"<<<{min_n}>>>")