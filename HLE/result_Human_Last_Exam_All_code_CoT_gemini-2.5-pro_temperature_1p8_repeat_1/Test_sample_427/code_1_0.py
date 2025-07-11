import scipy.stats as stats

def find_minimum_observations():
    """
    Calculates the minimum number of observations per group to achieve a
    Bonferroni-corrected p-value < 0.05 for 5 tests using a two-sided
    Mann-Whitney U test.
    """
    original_alpha = 0.05
    num_tests = 5
    corrected_alpha = original_alpha / num_tests

    print(f"The original significance level (alpha) is {original_alpha}.")
    print(f"Number of tests to be performed is {num_tests}.")
    print(f"Applying Bonferroni correction, the new significance level is: {original_alpha} / {num_tests} = {corrected_alpha}")
    print("-" * 30)
    print("Searching for the minimum number of observations (n) per group...")

    # Start with n=2 since n=1 is not a meaningful sample for this test
    n = 2
    while True:
        # For a given n, the strongest effect (and lowest p-value) comes from
        # the most separated samples. We use Likert scores 1 and 5.
        group1 = [1] * n
        group2 = [5] * n

        # Perform a two-sided Mann-Whitney U test.
        # The 'exact' method is used for small sample sizes to get precise p-values.
        u_statistic, p_value = stats.mannwhitneyu(
            group1,
            group2,
            alternative='two-sided',
            method='exact'
        )

        # Check if the p-value meets our corrected significance level
        if p_value < corrected_alpha:
            print(f"Found minimum n = {n}")
            print(f"With n = {n}, we create two perfectly separated groups to test the best-case scenario:")
            print(f"Group 1: {group1}")
            print(f"Group 2: {group2}")
            print(f"The two-sided Mann-Whitney U test yields a p-value of {p_value:.5f}.")
            print("\nThis p-value is less than the corrected alpha.")
            print("The final comparison is:")
            # The user requested to output each number in the final equation.
            print(f"{p_value:.5f} < {original_alpha} / {num_tests}")
            return n

        # If not significant, increment n and continue the search
        n += 1

if __name__ == '__main__':
    min_n = find_minimum_observations()
    print(f"<<<{min_n}>>>")
