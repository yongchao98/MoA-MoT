import scipy.stats

def find_minimum_sample_size():
    """
    Calculates the minimum number of observations per group to achieve
    statistical significance for a Mann-Whitney U test with Bonferroni correction.
    """
    # --- Part 1: Define problem parameters ---
    original_alpha = 0.05
    num_tests = 5
    group1_pattern = [1, 2, 3, 4, 5]
    group2_pattern = [3, 3, 3, 2, 2]

    # --- Part 2: Calculate Bonferroni corrected alpha ---
    corrected_alpha = original_alpha / num_tests
    
    print("--- Problem Setup ---")
    print(f"Test: Two-sided Mann-Whitney U test")
    print(f"Group 1 data pattern: {group1_pattern}")
    print(f"Group 2 data pattern: {group2_pattern}")
    print("\n--- Bonferroni Correction ---")
    print("The final equation for the corrected p-value is:")
    print(f"{original_alpha} (original alpha) / {num_tests} (number of tests) = {corrected_alpha}")
    print(f"We need to find the minimum sample size 'n' per group for the p-value to be < {corrected_alpha}.")
    print("\n--- Searching for Minimum Sample Size ---")

    # --- Part 3: Iterate to find the minimum sample size n ---
    n = 5  # Start from the initial number of observations
    while True:
        # Generate data for the current sample size 'n' by repeating the pattern
        repeats = n // len(group1_pattern)
        remainder = n % len(group1_pattern)
        group1_data = group1_pattern * repeats + group1_pattern[:remainder]
        group2_data = group2_pattern * repeats + group2_pattern[:remainder]

        # Perform the two-sided Mann-Whitney U test
        res = scipy.stats.mannwhitneyu(group1_data, group2_data, alternative='two-sided')
        p_value = res.pvalue
        
        # Check if the p-value is statistically significant
        if p_value < corrected_alpha:
            print(f"Found! With n = {n}, the p-value is {p_value:.6f}, which is less than {corrected_alpha}.")
            print("\n--- Conclusion ---")
            print("The minimum number of observations per group to achieve statistical significance is:")
            print(n)
            # Final answer format
            print(f"\n<<<{n}>>>")
            break
        
        # Add a safeguard to prevent an excessively long or infinite loop
        if n > 500:
             print("Could not find a solution for n up to 500. Please check the problem parameters.")
             break
        
        n += 1

# Execute the function to find the answer
find_minimum_sample_size()