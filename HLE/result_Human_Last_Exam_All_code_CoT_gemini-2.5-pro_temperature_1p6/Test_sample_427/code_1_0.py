import scipy.stats

def find_minimum_observations():
    """
    Calculates the minimum number of observations per group to achieve
    statistical significance with Bonferroni correction.
    """
    # --- Step 1: Define parameters and calculate the corrected p-value ---
    original_alpha = 0.05
    num_tests = 5
    
    # This is the final equation for the corrected alpha
    corrected_alpha = original_alpha / num_tests
    
    print("--- Bonferroni Correction ---")
    print("The original desired p-value (alpha) is 0.05.")
    print("The number of statistical tests is 5.")
    print("The equation for the corrected alpha is: original_alpha / num_tests")
    print(f"The resulting corrected alpha is: {original_alpha} / {num_tests} = {corrected_alpha}")
    print("\nWe must find the minimum sample size where p-value < 0.01.\n")

    # --- Step 2: Define the data patterns for each group ---
    group1_pattern = [1, 2, 3, 4, 5]
    group2_pattern = [3, 3, 3, 2, 2]

    # --- Step 3: Iteratively increase sample size and test for significance ---
    print("--- Searching for Minimum Sample Size (n) ---")
    multiplier = 1
    while True:
        n_per_group = len(group1_pattern) * multiplier
        
        sample1 = group1_pattern * multiplier
        sample2 = group2_pattern * multiplier
        
        # Perform the two-sided Mann-Whitney U test
        u_statistic, p_value = scipy.stats.mannwhitneyu(sample1, sample2, alternative='two-sided')
        
        # Check if significance is met
        if p_value < corrected_alpha:
            print(f"Testing n = {n_per_group} per group... Significance FOUND.")
            print("\n--- Final Result ---")
            print(f"The minimum number of observations per group required is {n_per_group}.")
            print(f"With n = {n_per_group}, the Mann-Whitney U test yields a p-value of {p_value:.5f}.")
            print("\nThe final equation for the significance check is: p-value < corrected_alpha")
            print(f"Numerically: {p_value:.5f} < {corrected_alpha}")
            
            # The final answer to be extracted
            final_answer = n_per_group
            return final_answer
        else:
            print(f"Testing n = {n_per_group} per group... p-value={p_value:.5f}, significance not met.")

        multiplier += 1

if __name__ == '__main__':
    minimum_n = find_minimum_observations()
    # The final answer is enclosed in <<<>>> as requested.
    # The value is printed inside the function, but we add the tag here.
    print(f"\n<<<>>>\n{minimum_n}\n<<<>>>")
