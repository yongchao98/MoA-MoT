import scipy.stats as stats

def find_minimum_sample_size():
    """
    Calculates the minimum number of observations per group for a two-sided
    Mann-Whitney U test to be significant after Bonferroni correction.
    """
    # Base samples for each group
    g1_base = [1, 2, 3, 4, 5]
    g2_base = [3, 3, 3, 2, 2]

    # --- Step 1: Calculate Bonferroni corrected alpha ---
    original_alpha = 0.05
    num_tests = 5
    corrected_alpha = original_alpha / num_tests

    # --- Step 2 & 3: Iterate to find the minimum sample size ---
    n = 5  # Start with the initial number of observations
    while True:
        # Construct samples of size n by replicating the base pattern
        num_replications = n // 5
        remainder = n % 5
        
        group1 = g1_base * num_replications + g1_base[:remainder]
        group2 = g2_base * num_replications + g2_base[:remainder]
        
        # Perform the two-sided Mann-Whitney U test
        u_statistic, p_value = stats.mannwhitneyu(group1, group2, alternative='two-sided')

        # Check if the p-value is less than the corrected alpha
        if p_value < corrected_alpha:
            print("To solve this, we first calculate the Bonferroni-corrected significance level (alpha).")
            print(f"The original desired p-value was < {original_alpha}.")
            print("The final equation for the corrected alpha is:")
            print(f"{original_alpha} / {num_tests} = {corrected_alpha}")
            
            print("\nNext, we iterate, increasing the sample size per group (n) until the p-value from the Mann-Whitney U test is less than this corrected alpha.")
            print("\n-------------------------------------------")
            print(f"FOUND: Minimum sample size is n = {n}")
            print("-------------------------------------------")

            print(f"For n = {n}, the groups are:")
            print(f"Group 1 (n={n}): {group1}")
            print(f"Group 2 (n={n}): {group2}")
            
            print(f"\nThe Mann-Whitney U test for n={n} yields a p-value of {p_value:.6f}.")
            print(f"Since {p_value:.6f} < {corrected_alpha}, this is the minimum required number of observations per group to achieve statistical significance.")
            
            return n

        # Safety break to prevent an infinite loop in case of an issue
        if n > 1000:
            print("Could not find a solution for n up to 1000. Please check the parameters.")
            return None

        n += 1

# Run the function to find and print the result.
minimum_n = find_minimum_sample_size()

# The final answer must be a single value, as requested by the user format.
# We extract just the numerical value for the final output.
# print(f"\n<<< {minimum_n} >>>") # This is for internal check, will not be in the final output block.

if __name__ == '__main__':
    pass # The function call above will print the output.
