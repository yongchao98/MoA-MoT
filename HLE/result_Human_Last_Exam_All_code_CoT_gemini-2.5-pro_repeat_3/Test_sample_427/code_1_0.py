import scipy.stats

def find_minimum_sample_size():
    """
    Calculates the minimum number of observations per group to achieve
    statistical significance with a two-sided Mann-Whitney U test after
    Bonferroni correction.
    """
    # Step 1: Define parameters and calculate the Bonferroni-corrected alpha
    initial_alpha = 0.05
    num_tests = 5
    corrected_alpha = initial_alpha / num_tests

    # Step 2: Define the base observations for each group
    group1_base = [1, 2, 3, 4, 5]
    group2_base = [3, 3, 3, 2, 2]

    print("Searching for the minimum number of observations per group (n)...")
    print("-" * 60)
    print(f"The initial desired p-value was < {initial_alpha}.")
    print(f"With {num_tests} tests, the Bonferroni corrected p-value threshold is {initial_alpha} / {num_tests} = {corrected_alpha}.")
    print("-" * 60)

    # Step 3: Iterate through sample sizes, starting from the base size of 5
    # We set a reasonable upper limit (e.g., 500) to prevent an infinite loop.
    for n in range(len(group1_base), 501):
        # Step 4: Generate the samples for the current 'n' by repeating the base patterns
        group1 = group1_base * (n // len(group1_base)) + group1_base[:n % len(group1_base)]
        group2 = group2_base * (n // len(group2_base)) + group2_base[:n % len(group2_base)]

        # Step 5: Perform the two-sided Mann-Whitney U test
        u_statistic, p_value = scipy.stats.mannwhitneyu(group1, group2, alternative='two-sided')

        # Step 6: Check if the p-value is below the corrected significance level
        if p_value < corrected_alpha:
            # Step 7: If significance is achieved, print the results and stop.
            print(f"Success! Found the minimum number of observations.\n")
            print(f"The minimum number of observations per group required is: {n}")
            print(f"For n = {n}, the Mann-Whitney U test yields a p-value of {p_value:.6f}.")
            
            print("\nFinal Equation:")
            # The final equation showing the comparison
            print(f"{p_value:.6f} < {corrected_alpha}")
            
            # This is the final answer for automated extraction
            print(f"\n<<<{n}>>>")
            return

    # This part will be reached if no solution is found within the loop's range
    print("Could not find a sufficient sample size within the tested range (up to n=500).")

# Execute the function
find_minimum_sample_size()