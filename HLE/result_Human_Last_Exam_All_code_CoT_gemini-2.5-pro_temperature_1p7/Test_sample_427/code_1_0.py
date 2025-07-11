import scipy.stats as stats

def find_min_observations():
    """
    Calculates the minimum number of observations per group needed to achieve
    statistical significance with a Bonferroni-corrected p-value in a
    two-sided Mann-Whitney U test.
    """
    original_alpha = 0.05
    num_tests = 5
    corrected_alpha = original_alpha / num_tests

    print(f"The original alpha was {original_alpha}.")
    print(f"With {num_tests} tests, the Bonferroni-corrected alpha is {original_alpha} / {num_tests} = {corrected_alpha}.")
    print("-" * 30)
    
    n = 5  # Start with the initial number of observations
    while True:
        # Generate Group 1 as a sequence from 1 to n
        group1 = list(range(1, n + 1))
        
        # Generate Group 2 maintaining a 40/60 split of 2s and 3s
        num_twos = round(0.4 * n)
        num_threes = n - num_twos
        group2 = [2] * num_twos + [3] * num_threes
        
        # Perform the two-sided Mann-Whitney U test
        # Use continuity correction, which is the default behavior
        stat, p_value = stats.mannwhitneyu(group1, group2, alternative='two-sided')
        
        # Check if the p-value meets the corrected significance level
        if p_value < corrected_alpha:
            print(f"Significance achieved at n = {n} observations per group.")
            print(f"\nFinal Groups (n={n}):")
            print(f"Group 1: {group1}")
            print(f"Group 2: {group2}")
            print(f"\nThe calculated p-value is {p_value:.6f}")
            print("\nFinal Equation:")
            print(f"{p_value:.6f} < {corrected_alpha}")
            return n
        
        n += 1

if __name__ == '__main__':
    min_n = find_min_observations()
    # The final result is requested in a specific format at the end of the response.
    # The printed output from the function provides the detailed answer.
    # The value returned here is for generating the final response tag.

