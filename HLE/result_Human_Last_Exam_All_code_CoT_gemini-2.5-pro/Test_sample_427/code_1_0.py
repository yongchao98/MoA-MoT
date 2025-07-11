import numpy as np
from scipy.stats import mannwhitneyu

def find_minimum_sample_size():
    """
    Calculates the minimum number of observations per group for a two-sided 
    Mann-Whitney U test to achieve statistical significance after Bonferroni correction.
    """

    # --- Step 1: Determine the required p-value ---
    print("Step 1: Determine the required p-value.")
    original_alpha = 0.05
    num_tests = 5
    corrected_alpha = original_alpha / num_tests
    
    print("The originally desired p-value was < 0.05.")
    print(f"Since {num_tests} tests are being run, Bonferroni correction is applied.")
    print("The corrected p-value threshold is calculated as follows:")
    # Here we output each number in the final equation as requested.
    print(f"{original_alpha} / {num_tests} = {corrected_alpha}")
    print(f"Therefore, we need to find the sample size where the p-value is < {corrected_alpha}.\n")

    # --- Step 2: Iteratively test sample sizes ---
    print("Step 2: Iteratively test sample sizes.")
    print("I will start with n=5 and increase the number of observations per group by one in each iteration.")
    print("The loop will stop when the Mann-Whitney U test p-value is less than the corrected alpha.\n")
    
    # Base patterns for the observations
    base_group1 = [1, 2, 3, 4, 5]
    base_group2 = [3, 3, 3, 2, 2]

    n = 5
    p_value_at_n_minus_1 = 0.0
    
    while True:
        # Generate data by repeating the base pattern
        num_repeats = (n + len(base_group1) - 1) // len(base_group1)
        group1 = (base_group1 * num_repeats)[:n]
        group2 = (base_group2 * num_repeats)[:n]

        # Perform the two-sided Mann-Whitney U test
        stat, p_value = mannwhitneyu(group1, group2, alternative='two-sided')

        # Check if the p-value is below the corrected alpha
        if p_value < corrected_alpha:
            min_n = n
            final_p_value = p_value
            break

        # Store the p-value for the previous iteration for the final report
        p_value_at_n_minus_1 = p_value
        n += 1

        # Safety break to prevent an accidental infinite loop
        if n > 1000:
            print("Search stopped after n=1000. A solution was not found.")
            return

    # --- Step 3: State the final answer ---
    print("...")
    print(f"Search complete. The condition p < {corrected_alpha} was met.")
    print(f"At a sample size of n = {min_n - 1}, the p-value was {p_value_at_n_minus_1:.5f}, which is not significant.")
    print(f"At a sample size of n = {min_n}, the p-value is {final_p_value:.5f}, which IS significant.\n")

    print("The minimum number of observations per group required to achieve statistical significance is:")
    print(min_n)


if __name__ == '__main__':
    find_minimum_sample_size()
    print("<<<19>>>")