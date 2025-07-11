import numpy as np
from scipy import stats

def compare_extinction_points():
    """
    Analyzes and compares the extinction points of two microbial cell types
    to determine if they are significantly different.
    """
    # --- 1. Sample Data ---
    # These are the collected extinction times (e.g., in hours) from multiple
    # experimental replicates for two different cell types.
    # In a real scenario, you would replace this with your actual data.
    cell_type_A_extinction_times = [20.2, 22.5, 19.8, 24.1, 18.5, 21.7, 23.3, 19.1]
    cell_type_B_extinction_times = [25.5, 26.1, 28.3, 24.9, 23.8, 27.5, 26.8, 25.1, 27.9]

    print("--- Step 0: Input Data ---")
    print(f"Cell Type A Extinction Times (hours): {cell_type_A_extinction_times}")
    print(f"Cell Type B Extinction Times (hours): {cell_type_B_extinction_times}")
    print("-" * 50)

    # --- 2. Check for Normality using Shapiro-Wilk Test ---
    # The null hypothesis (H0) is that the data is normally distributed.
    # We choose a significance level (alpha) of 0.05 for this check.
    alpha_normality = 0.05
    shapiro_A_stat, shapiro_A_p = stats.shapiro(cell_type_A_extinction_times)
    shapiro_B_stat, shapiro_B_p = stats.shapiro(cell_type_B_extinction_times)

    is_A_normal = shapiro_A_p > alpha_normality
    is_B_normal = shapiro_B_p > alpha_normality

    print("--- Step 1: Check for Normality (Shapiro-Wilk Test) ---")
    print(f"Cell Type A: p-value = {shapiro_A_p:.4f}. Normal Distribution? {'Yes' if is_A_normal else 'No'}")
    print(f"Cell Type B: p-value = {shapiro_B_p:.4f}. Normal Distribution? {'Yes' if is_B_normal else 'No'}")
    print("-" * 50)

    # --- 3. Choose and Perform Significance Test ---
    # This logic directly follows answer choice D.
    alpha_significance = 0.05

    print("--- Step 2: Perform Significance Test ---")
    if is_A_normal and is_B_normal:
        print("Both datasets appear normal. Using Independent (Unpaired) T-test.")
        # Welch's T-test is used (equal_var=False) as it's more robust.
        test_statistic, p_value = stats.ttest_ind(
            cell_type_A_extinction_times,
            cell_type_B_extinction_times,
            equal_var=False
        )
        test_name = "Unpaired T-test"
        print(f"Comparing mean({np.mean(cell_type_A_extinction_times):.2f}) vs mean({np.mean(cell_type_B_extinction_times):.2f})")

    else:
        print("At least one dataset is not normal. Using Wilcoxon Rank-Sum Test.")
        # This is the non-parametric equivalent of the t-test.
        test_statistic, p_value = stats.mannwhitneyu(
            cell_type_A_extinction_times,
            cell_type_B_extinction_times,
            alternative='two-sided'
        )
        test_name = "Wilcoxon Rank-Sum Test"
        print(f"Comparing median({np.median(cell_type_A_extinction_times):.2f}) vs median({np.median(cell_type_B_extinction_times):.2f})")

    print(f"\nTest Used: {test_name}")
    print(f"Test Statistic: {test_statistic:.4f}")
    print(f"P-value: {p_value:.4f}")
    print("-" * 50)

    # --- 4. Interpretation ---
    print("--- Step 3: Conclusion ---")
    if p_value < alpha_significance:
        print(f"The p-value ({p_value:.4f}) is less than our significance level ({alpha_significance}).")
        print("Result: We reject the null hypothesis. The difference in extinction points is statistically significant.")
    else:
        print(f"The p-value ({p_value:.4f}) is not less than our significance level ({alpha_significance}).")
        print("Result: We fail to reject the null hypothesis. The difference is not statistically significant.")

if __name__ == '__main__':
    compare_extinction_points()