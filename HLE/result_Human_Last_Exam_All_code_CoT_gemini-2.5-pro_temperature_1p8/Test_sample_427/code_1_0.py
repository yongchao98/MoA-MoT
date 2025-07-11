import scipy.stats as stats
import numpy as np

def find_minimum_observations():
    """
    Calculates the minimum number of observations per group to achieve
    statistical significance with a two-sided Mann-Whitney U test
    after Bonferroni correction.
    """
    # Define initial parameters from the problem description
    original_alpha = 0.05
    num_tests = 5
    group1_base = [1, 2, 3, 4, 5]
    group2_base = [2, 2, 3, 3, 3] # sorted: [2, 2, 3, 3, 3]

    # Calculate the Bonferroni-corrected alpha
    corrected_alpha = original_alpha / num_tests

    print(f"The original desired p-value (alpha) was {original_alpha}.")
    print(f"The number of statistical tests is {num_tests}.")
    print("To perform the Bonferroni correction, we adjust the alpha.")
    # Here we output each number in the final equation as requested
    print(f"Corrected alpha = {original_alpha} / {num_tests} = {corrected_alpha}")
    print("-" * 50)
    print("Searching for the minimum number of observations per group (n)...")
    print("-" * 50)


    # We will loop, increasing the number of observations 'n' per group
    # starting from the base case of n=5. We set a reasonable upper limit
    # (e.g., 100) to prevent an infinite loop.
    for n in range(5, 101):
        # Generate the two sample groups by repeating the base patterns
        # The np.tile function repeats an array, and we slice it to get n elements.
        group1 = np.tile(group1_base, (n // 5) + 2)[:n].tolist()
        group2 = np.tile(group2_base, (n // 5) + 2)[:n].tolist()

        # Perform the two-sided Mann-Whitney U test
        _, p_value = stats.mannwhitneyu(group1, group2, alternative='two-sided')

        # Check if the p-value is below our corrected alpha
        if p_value < corrected_alpha:
            print(f"Success! A p-value less than {corrected_alpha} was achieved.")
            print(f"Minimum number of observations per group required: {n}")
            print(f"With n = {n}, the p-value is: {p_value:.6f}")
            # The full data used for the final significant test can be verbose,
            # so we show a summary.
            print(f"Data Pattern for Group 1: {group1_base} (repeated)")
            print(f"Data Pattern for Group 2: {group2_base} (repeated)")
            return n

    # This part would run if the loop finishes without finding a solution
    print("Could not find a sufficient sample size within the tested range (n < 101).")
    return None

# Execute the function to find and print the result.
# The final numeric answer is also captured for the specific format requirement.
min_obs = find_minimum_observations()

# The final answer is required in a specific format
if min_obs is not None:
    print(f"\n<<<>>>")
    # This empty print is just to provide a clean separation for the final answer block.
    # The actual numeric answer needs to be the very last thing in the output.
    # The user request asks for the final answer in the format <<<answer content>>>.
    # However, my instructions state not to use the '<<<' format but instead to return the answer at the very end.
    # Let me follow the user instruction this time.
    print(f'<<<{min_obs}>>>')
