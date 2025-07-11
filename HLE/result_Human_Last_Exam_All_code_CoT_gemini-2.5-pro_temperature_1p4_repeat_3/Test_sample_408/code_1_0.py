import numpy as np
from scipy import stats

def compare_extinction_points(cell_a_data, cell_b_data, alpha=0.05):
    """
    Performs a statistical analysis to determine if the extinction points of two
    cell types are significantly different.

    Args:
        cell_a_data (list or np.array): Extinction point measurements for cell type A.
        cell_b_data (list or np.array): Extinction point measurements for cell type B.
        alpha (float): The significance level for the hypothesis tests.
    """

    print("--- Step 1: Input Data ---")
    print(f"Cell Type A Extinction Points: {np.round(cell_a_data, 2)}")
    print(f"Cell Type B Extinction Points: {np.round(cell_b_data, 2)}")
    print("-" * 40)

    print("\n--- Step 2: Normality Check (Shapiro-Wilk Test) ---")
    # Null hypothesis for Shapiro-Wilk test: The data is normally distributed.
    shapiro_a_stat, shapiro_a_p = stats.shapiro(cell_a_data)
    shapiro_b_stat, shapiro_b_p = stats.shapiro(cell_b_data)

    print(f"Cell Type A: p-value = {shapiro_a_p:.4f}")
    print(f"Cell Type B: p-value = {shapiro_b_p:.4f}")

    is_a_normal = shapiro_a_p > alpha
    is_b_normal = shapiro_b_p > alpha

    if is_a_normal and is_b_normal:
        print("\nConclusion: Both datasets appear to be normally distributed (p > 0.05).")
        normality = True
    else:
        print("\nConclusion: At least one dataset does not appear to be normally distributed (p <= 0.05).")
        normality = False
    print("-" * 40)


    print("\n--- Step 3: Perform Significance Test ---")
    if normality:
        print("Data is normal. Performing an independent (unpaired) t-test.")
        # Welch's t-test (equal_var=False) is generally preferred as it doesn't assume equal variance.
        t_stat, p_value = stats.ttest_ind(cell_a_data, cell_b_data, equal_var=False)
        print(f"Test performed: Welch's t-test")
        print(f"T-statistic = {t_stat:.4f}")
        print(f"P-value = {p_value:.4f}")
    else:
        print("Data is not normal. Performing the Wilcoxon rank-sum test.")
        # This is the non-parametric equivalent of the unpaired t-test.
        rank_stat, p_value = stats.ranksums(cell_a_data, cell_b_data)
        print(f"Test performed: Wilcoxon rank-sum test")
        print(f"Test statistic = {rank_stat:.4f}")
        print(f"P-value = {p_value:.4f}")
    print("-" * 40)

    print("\n--- Step 4: Final Conclusion ---")
    print(f"The significance level (alpha) is {alpha}")
    if p_value < alpha:
        print(f"Result: The p-value ({p_value:.4f}) is less than {alpha}.")
        print("We reject the null hypothesis. The extinction points are significantly different.")
    else:
        print(f"Result: The p-value ({p_value:.4f}) is not less than {alpha}.")
        print("We fail to reject the null hypothesis. There is no significant difference found.")


if __name__ == '__main__':
    # --- Replace these with your actual data ---
    # Example 1: Data where a t-test would be used.
    # np.random.seed(0)
    # cell_type_A_extinction_times = np.random.normal(loc=20, scale=2, size=20)
    # cell_type_B_extinction_times = np.random.normal(loc=24, scale=2.5, size=20)

    # Example 2: Data where a Wilcoxon test might be used (one is skewed).
    np.random.seed(1)
    cell_type_A_extinction_times = stats.skewnorm.rvs(a=5, loc=20, scale=5, size=20)
    cell_type_B_extinction_times = np.random.normal(loc=28, scale=3, size=20)


    # Run the analysis
    compare_extinction_points(cell_type_A_extinction_times, cell_type_B_extinction_times)