import numpy as np
from scipy import stats

def analyze_extinction_points(data_a, data_b, alpha=0.05):
    """
    Analyzes if the extinction points of two cell types are significantly different.

    This function follows the standard statistical procedure:
    1. Check if the data for each cell type is normally distributed using the Shapiro-Wilk test.
    2. If both are normal, perform an unpaired t-test.
    3. If at least one is not normal, perform a Wilcoxon rank-sum (Mann-Whitney U) test.
    4. Print the results and a conclusion.

    Args:
        data_a (list or np.array): Extinction point data for cell type A.
        data_b (list or np.array): Extinction point data for cell type B.
        alpha (float): The significance level for the tests.
    """
    print("--- Data ---")
    print(f"Cell Type A Extinction Times (hours): {np.round(data_a, 2)}")
    print(f"Cell Type B Extinction Times (hours): {np.round(data_b, 2)}\n")

    # Step 1: Check for normality using Shapiro-Wilk test
    print(f"--- Step 1: Normality Check (at alpha={alpha}) ---")
    shapiro_a_stat, shapiro_a_p = stats.shapiro(data_a)
    shapiro_b_stat, shapiro_b_p = stats.shapiro(data_b)

    print(f"Shapiro-Wilk test for Cell Type A: p-value = {shapiro_a_p:.4f}")
    print(f"Shapiro-Wilk test for Cell Type B: p-value = {shapiro_b_p:.4f}")

    # Step 2: Decide which test to use
    is_a_normal = shapiro_a_p > alpha
    is_b_normal = shapiro_b_p > alpha

    if is_a_normal and is_b_normal:
        print("\nConclusion: Both datasets appear to be normally distributed.")
        print("--- Step 2: Performing Unpaired T-test ---\n")
        # Perform unpaired t-test
        test_stat, p_value = stats.ttest_ind(data_a, data_b, equal_var=True) # Assuming equal variances for simplicity
        test_name = "Unpaired T-test"
    else:
        print("\nConclusion: At least one dataset does not appear to be normally distributed.")
        print("--- Step 2: Performing Wilcoxon Rank-Sum Test ---\n")
        # Perform Wilcoxon rank-sum test (Mann-Whitney U)
        test_stat, p_value = stats.mannwhitneyu(data_a, data_b, alternative='two-sided')
        test_name = "Wilcoxon Rank-Sum Test"

    # Step 3: Print results and interpret
    print(f"--- Step 3: Results ---")
    print(f"Test Used: {test_name}")
    print(f"Test Statistic: {test_stat:.4f}")
    print(f"P-value: {p_value:.4f}")

    if p_value < alpha:
        print(f"\nConclusion: Since the p-value ({p_value:.4f}) is less than {alpha}, we reject the null hypothesis.")
        print("The extinction points of the two cell types are significantly different.")
    else:
        print(f"\nConclusion: Since the p-value ({p_value:.4f}) is greater than {alpha}, we fail to reject the null hypothesis.")
        print("There is not enough evidence to say the extinction points are significantly different.")

# --- Example 1: Normally Distributed Data ---
print("==============================================")
print("Example 1: Data appears normally distributed")
print("==============================================")
# Simulate data for time to extinction (in hours) for 10 replicates of each cell type.
np.random.seed(0) # for reproducibility
# Cell type A: mean extinction time of 48 hours
cell_type_A_normal = np.random.normal(loc=48, scale=5, size=10)
# Cell type B: mean extinction time of 60 hours
cell_type_B_normal = np.random.normal(loc=60, scale=5, size=10)
analyze_extinction_points(cell_type_A_normal, cell_type_B_normal)


# --- Example 2: Non-Normally Distributed Data ---
print("\n==============================================")
print("Example 2: Data is not normally distributed")
print("==============================================")
np.random.seed(42)
# Cell type A: from an exponential distribution
cell_type_A_non_normal = np.random.exponential(scale=50, size=15)
# Cell type B: from an exponential distribution with a different scale
cell_type_B_non_normal = np.random.exponential(scale=80, size=15)
analyze_extinction_points(cell_type_A_non_normal, cell_type_B_non_normal)
