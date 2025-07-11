import numpy as np
from scipy.stats import mannwhitneyu
import math

def solve_min_sample_size():
    """
    Calculates the minimum sample size per group to achieve statistical
    significance with a Bonferroni-corrected p-value.
    """
    original_alpha = 0.05
    num_tests = 5
    corrected_alpha = original_alpha / num_tests

    # We start checking from n=6 since n=5 is given and is not significant.
    # We set an upper limit to prevent an infinite loop, though it should be found quickly.
    for n in range(6, 101):
        # Generate Group 1 data: repeating pattern of [1, 2, 3, 4, 5]
        base_g1 = [1, 2, 3, 4, 5]
        group1 = np.tile(base_g1, n // len(base_g1)).tolist() + base_g1[:n % len(base_g1)]

        # Generate Group 2 data: 60% are '3', 40% are '2'
        num_threes = round(n * 0.6)
        num_twos = n - num_threes
        group2 = [3] * int(num_threes) + [2] * int(num_twos)

        # Perform the two-sided Mann-Whitney U test
        # Use method='asymptotic' for continuity correction, which is default for n > 20
        # or when ties are present. 'auto' will handle this correctly.
        u_statistic, p_value = mannwhitneyu(group1, group2, alternative='two-sided')

        # Check if the p-value is below the corrected alpha
        if p_value < corrected_alpha:
            print(f"The minimum number of observations per group is: {n}")
            print("\nThis is the first sample size where the p-value is less than the Bonferroni-corrected alpha.")
            print("\nFinal Equation:")
            # Output each number in the final equation as requested
            print(f"p-value ({p_value:.4f}) < corrected alpha ({corrected_alpha})")
            return n

    print("No sample size found within the tested range (up to 100).")
    return None

solve_min_sample_size()