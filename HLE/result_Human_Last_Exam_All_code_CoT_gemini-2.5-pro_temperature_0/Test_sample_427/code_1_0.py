import scipy.stats

def find_min_observations():
    """
    This function calculates the minimum number of observations per group
    for a two-sided Mann-Whitney U test to be able to achieve a
    statistically significant result, considering a Bonferroni correction.

    The problem specifies:
    - Original alpha: 0.05
    - Number of tests: 5
    - This requires a corrected alpha for significance.
    """
    original_alpha = 0.05
    num_tests = 5
    corrected_alpha = original_alpha / num_tests

    print("This script determines the minimum number of observations per group ('n') required")
    print("to achieve statistical significance with a Bonferroni-corrected p-value.")
    print("-" * 60)

    # Step 1: Explain and show the Bonferroni correction
    print("Step 1: Calculate the corrected p-value threshold.")
    print(f"The original alpha was {original_alpha} and the number of tests is {num_tests}.")
    # Outputting each number in the equation as requested
    print(f"The required p-value is calculated as: {original_alpha} / {num_tests} = {corrected_alpha}")
    print("-" * 60)

    # Step 2: Iterate to find the minimum sample size
    print(f"Step 2: Find the smallest 'n' where a p-value < {corrected_alpha} is possible.")
    print("The smallest possible p-value for a given 'n' occurs with perfect separation of ranks.")
    print("We will now test sample sizes incrementally.")
    print("-" * 60)

    # We will test n from 2 up to a reasonable limit (e.g., 20)
    for n in range(2, 21):
        # Create two perfectly separated groups for the best-case scenario.
        # The actual values don't matter, only their ranks.
        # e.g., for n=4: group1=[1,2,3,4], group2=[5,6,7,8]
        group1 = list(range(1, n + 1))
        group2 = list(range(n + 1, 2 * n + 1))

        # Perform the two-sided Mann-Whitney U test
        # For small samples, scipy calculates the exact p-value.
        u_stat, p_value = scipy.stats.mannwhitneyu(group1, group2, alternative='two-sided')

        print(f"Testing with n = {n} observations per group...")
        print(f"   - The smallest possible p-value is {p_value:.4f}")

        # Check if this p-value meets the corrected significance level
        if p_value < corrected_alpha:
            print(f"\nAt n = {n}, the p-value ({p_value:.4f}) is less than the threshold ({corrected_alpha}).")
            print("-" * 60)
            print(f"Therefore, the minimum number of observations per group is {n}.")
            # We found the minimum, so we can stop.
            return n
        else:
            print(f"   - This is not less than {corrected_alpha}.")

    # This part would be reached if no solution is found in the loop range
    print("No solution found within the tested range of n.")
    return None

# Run the function to find the answer
min_n_result = find_min_observations()

if min_n_result is not None:
    print(f"\n<<<" + str(min_n_result) + ">>>")