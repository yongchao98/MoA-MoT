import scipy.stats

def find_minimum_observations():
    """
    Calculates the minimum number of observations per group to achieve a
    statistically significant result in a Mann-Whitney U test after
    Bonferroni correction.
    """
    # --- 1. Define initial parameters ---
    # The fixed scores obtained in one session for each group
    group1_base_scores = [1, 2, 3, 4, 5]
    group2_base_scores = [3, 3, 3, 2, 2]

    # The original significance level and the number of tests for correction
    original_alpha = 0.05
    num_tests = 5

    # --- 2. Calculate the corrected significance level (alpha) ---
    corrected_alpha = original_alpha / num_tests

    print(f"The original desired p-value was < {original_alpha}")
    print(f"A Bonferroni correction is needed for {num_tests} tests.")
    print(f"The new significance threshold is: {original_alpha} / {num_tests} = {corrected_alpha}")
    print("-" * 50)

    # --- 3. Iteratively find the required sample size ---
    num_sessions = 0
    p_value = 1.0  # Initialize p-value to a value greater than alpha

    while p_value >= corrected_alpha:
        num_sessions += 1
        
        # Create the pooled data for the current number of sessions
        current_group1 = group1_base_scores * num_sessions
        current_group2 = group2_base_scores * num_sessions
        n_obs_per_group = len(current_group1)

        # Perform the two-sided Mann-Whitney U test
        u_statistic, p_value = scipy.stats.mannwhitneyu(
            current_group1,
            current_group2,
            alternative='two-sided'
        )

        print(f"Test with {num_sessions} session(s) ({n_obs_per_group} observations per group)...")
        print(f"Resulting p-value: {p_value:.5f}")

        if p_value < corrected_alpha:
            print(f"SUCCESS: The p-value ({p_value:.5f}) is less than the required threshold ({corrected_alpha}).")
        else:
            print(f"FAIL: The p-value ({p_value:.5f}) is not less than the required threshold ({corrected_alpha}).")
        print("-" * 50)

    # --- 4. Output the final result ---
    min_obs_per_group = num_sessions * len(group1_base_scores)
    print(f"\nThe minimum number of observations per group to achieve statistical significance is {min_obs_per_group}.")

if __name__ == '__main__':
    find_minimum_observations()
<<<30>>>