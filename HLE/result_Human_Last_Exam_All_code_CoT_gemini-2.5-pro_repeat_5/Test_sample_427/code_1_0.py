import scipy.stats

def find_minimum_observations():
    """
    Calculates the minimum number of observations per group to achieve
    statistical significance with Bonferroni correction.
    """
    # --- Step 1: Define Test Parameters ---

    # The observations obtained in a single session for each group.
    base_group1 = [1, 2, 3, 4, 5]
    base_group2 = [3, 3, 3, 2, 2]

    # Original desired p-value (alpha) and number of tests for Bonferroni correction.
    original_alpha = 0.05
    num_tests = 5

    # --- Step 2: Calculate Bonferroni Corrected Alpha ---

    # The corrected p-value threshold is the original alpha divided by the number of tests.
    corrected_alpha = original_alpha / num_tests
    print("Calculating the Bonferroni-corrected significance level (alpha):")
    print(f"Original alpha / Number of tests = Corrected alpha")
    print(f"{original_alpha} / {num_tests} = {corrected_alpha}\n")


    # --- Step 3: Iteratively Find the Required Sample Size ---

    print("Searching for the minimum number of observations per group...")
    # Start with one session (k=1) and increase until significance is met.
    k = 1
    while True:
        # Construct the full data sets by repeating the base observations k times.
        current_group1 = base_group1 * k
        current_group2 = base_group2 * k
        n_per_group = len(current_group1)

        # Perform the two-sided Mann-Whitney U test.
        # 'auto' method correctly handles ties and switches to normal approximation for large samples.
        u_statistic, p_value = scipy.stats.mannwhitneyu(
            current_group1,
            current_group2,
            alternative='two-sided'
        )

        # Check if the p-value is below our corrected threshold.
        if p_value < corrected_alpha:
            print(f"Success! Significance found at k = {k} sessions.")
            print(f"With {n_per_group} observations per group, the p-value is {p_value:.5f}, which is less than {corrected_alpha}.")
            print("\n--- Final Answer ---")
            print(f"The minimum number of observations per group required is: {n_per_group}")
            return n_per_group

        # Increment the session multiplier and continue the search.
        k += 1

        # A safety break to prevent an accidental infinite loop.
        if k > 200:
            print("Search stopped: Could not find a solution within 200 iterations.")
            break

# Run the function to find the answer.
find_minimum_observations()