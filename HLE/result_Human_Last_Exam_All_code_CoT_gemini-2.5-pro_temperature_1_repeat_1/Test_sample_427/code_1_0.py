import scipy.stats

def find_minimum_sample_size():
    """
    Calculates the minimum number of observations per group to achieve
    statistical significance for a Mann-Whitney U test after Bonferroni correction.
    """
    # 1. Define problem parameters
    base_group1 = [1, 2, 3, 4, 5]
    base_group2 = [3, 3, 3, 2, 2]
    num_tests = 5
    original_alpha = 0.05

    # 2. Calculate the Bonferroni-corrected alpha
    corrected_alpha = original_alpha / num_tests

    # 3. Iterate, increasing the sample size multiplier `k`
    # The sample size `n` will be 5 * k.
    k = 1
    while True:
        n = len(base_group1) * k
        current_group1 = base_group1 * k
        current_group2 = base_group2 * k

        # Perform the two-sided Mann-Whitney U test
        u_statistic, p_value = scipy.stats.mannwhitneyu(
            current_group1,
            current_group2,
            alternative='two-sided'
        )

        # 4. Check if the p-value meets the corrected significance level
        if p_value < corrected_alpha:
            print(f"The original significance level is alpha = {original_alpha}.")
            print(f"After Bonferroni correction for {num_tests} tests, the new significance level is {original_alpha} / {num_tests} = {corrected_alpha}.")
            print(f"\nStarting with n=5 observations per group, we increase the sample size until the p-value is less than {corrected_alpha}.")
            print(f"\nAt n = {n} observations per group, the p-value becomes {p_value:.6f}.")
            print("This is the minimum sample size where the test is statistically significant.")
            
            print("\nFinal Equation:")
            # Output each number in the final equation
            print(f"p-value ({p_value:.6f}) < corrected alpha ({corrected_alpha})")
            
            print(f"\nThe minimum number of observations required per group is {n}.")
            
            # Final answer in the required format
            print(f"\n<<< {n} >>>")
            break
        
        # Increment the multiplier to increase the sample size for the next iteration
        k += 1

if __name__ == '__main__':
    find_minimum_sample_size()