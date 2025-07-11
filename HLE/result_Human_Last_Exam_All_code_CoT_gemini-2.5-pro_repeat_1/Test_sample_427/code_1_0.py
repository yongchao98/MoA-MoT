import numpy as np
from scipy.stats import mannwhitneyu

def find_min_observations():
    """
    Calculates the minimum number of observations per group to achieve statistical
    significance after Bonferroni correction for a Mann-Whitney U test.
    """
    # --- Step 1: Define parameters and calculate the corrected alpha ---
    original_alpha = 0.05
    num_tests = 5
    bonferroni_alpha = original_alpha / num_tests

    # The base observations added in each session
    base_group1 = [1, 2, 3, 4, 5]
    base_group2 = [3, 3, 3, 2, 2]

    # --- Step 2 & 3: Iteratively increase sessions and perform the test ---
    num_sessions = 0
    while True:
        num_sessions += 1
        
        # Build the full datasets for the current number of sessions
        group1 = base_group1 * num_sessions
        group2 = base_group2 * num_sessions
        
        # Perform the two-sided Mann-Whitney U test
        # 'auto' method correctly handles ties and sample size
        res = mannwhitneyu(group1, group2, alternative='two-sided')
        p_value = res.pvalue

        # --- Step 4: Check if significance is met ---
        if p_value < bonferroni_alpha:
            n_per_group = len(group1)
            
            # --- Output the results ---
            print("To achieve the desired statistical significance, we must adjust the p-value using Bonferroni correction.")
            print("\nThe Bonferroni correction equation is:")
            print(f"Corrected Alpha = Original Alpha / Number of Tests")
            # The prompt asks to "output each number in the final equation"
            print(f"Corrected Alpha = {original_alpha} / {num_tests} = {bonferroni_alpha}\n")

            print(f"We need to find the minimum sample size where the p-value is < {bonferroni_alpha}.")
            print("-" * 50)
            print(f"The minimum number of observations per group required is: {n_per_group}")
            print("-" * 50)
            print(f"This is achieved with {num_sessions} repeated sessions.")
            print(f"At this sample size, the Mann-Whitney U test yields a p-value of {p_value:.5f}, which is less than the corrected alpha of {bonferroni_alpha}.")

            # Return the final answer as requested by the prompt format
            print(f"\n<<<{n_per_group}>>>")
            return

        # Safety break to prevent an infinite loop in case of an unexpected issue
        if num_sessions > 100:
            print("Could not find a significant p-value within 100 sessions (500 observations per group).")
            print(f"\n<<<-1>>>") # Indicate failure
            return

# Run the function
find_min_observations()