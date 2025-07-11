import numpy as np
from scipy import stats

def compare_extinction_points(data_a, data_b, alpha=0.05):
    """
    Compares the extinction points of two cell types using the appropriate statistical test.

    Args:
        data_a (list or np.array): Extinction time data for cell type A.
        data_b (list or np.array): Extinction time data for cell type B.
        alpha (float): The significance level for the statistical tests.
    """
    print(f"Data for Cell Type A: {data_a}")
    print(f"Data for Cell Type B: {data_b}\n")

    # Step 1: Check for normality using the Shapiro-Wilk test
    shapiro_a_stat, shapiro_a_p = stats.shapiro(data_a)
    shapiro_b_stat, shapiro_b_p = stats.shapiro(data_b)

    print("--- Step 1: Normality Check ---")
    print(f"Cell Type A: Shapiro-Wilk Test p-value = {shapiro_a_p:.4f}")
    print(f"Cell Type B: Shapiro-Wilk Test p-value = {shapiro_b_p:.4f}")

    # The null hypothesis of the Shapiro-Wilk test is that the data is normally distributed.
    # If p > alpha, we do not reject the null hypothesis, and assume normality.
    is_normal_a = shapiro_a_p > alpha
    is_normal_b = shapiro_b_p > alpha

    if is_normal_a:
        print("Data for Cell Type A appears to be normally distributed.")
    else:
        print("Data for Cell Type A does not appear to be normally distributed.")

    if is_normal_b:
        print("Data for Cell Type B appears to be normally distributed.")
    else:
        print("Data for Cell Type B does not appear to be normally distributed.")
    print("-" * 35 + "\n")


    # Step 2: Choose and perform the appropriate significance test
    print("--- Step 2: Significance Testing ---")
    if is_normal_a and is_normal_b:
        print("Both datasets appear normal. Performing an Unpaired T-test.")
        # Perform unpaired t-test
        t_stat, p_value = stats.ttest_ind(data_a, data_b, equal_var=True) # Assuming equal variances for simplicity
        print(f"T-statistic = {t_stat:.4f}")
        print(f"P-value = {p_value:.4f}")

    else:
        print("At least one dataset is not normal. Performing a Wilcoxon rank-sum (Mann-Whitney U) test.")
        # Perform Wilcoxon rank-sum test
        u_stat, p_value = stats.mannwhitneyu(data_a, data_b, alternative='two-sided')
        print(f"U-statistic = {u_stat:.4f}")
        print(f"P-value = {p_value:.4f}")
    print("-" * 35 + "\n")

    # Step 3: Interpret the result
    print("--- Step 3: Conclusion ---")
    if p_value < alpha:
        print(f"Since the p-value ({p_value:.4f}) is less than {alpha}, we reject the null hypothesis.")
        print("Conclusion: The extinction points of the two cell types are significantly different.")
    else:
        print(f"Since the p-value ({p_value:.4f}) is greater than or equal to {alpha}, we fail to reject the null hypothesis.")
        print("Conclusion: There is not enough evidence to say the extinction points are significantly different.")


# --- Example Usage ---
# Case 1: Normally distributed data
print("========== Example 1: Normally Distributed Data ==========")
# Simulated extinction times (in hours) for several replicates
cell_type_A_times_normal = [10.5, 11.2, 10.8, 12.1, 11.5, 10.9, 11.7]
cell_type_B_times_normal = [12.8, 13.5, 13.1, 12.9, 14.0, 13.3, 13.6]
compare_extinction_points(cell_type_A_times_normal, cell_type_B_times_normal)

print("\n\n========== Example 2: Non-Normally Distributed Data ==========")
# Case 2: Non-normally distributed data (e.g., with outliers or skew)
cell_type_A_times_non_normal = [10.1, 10.2, 10.3, 10.4, 10.5, 15.0] # Skewed by an outlier
cell_type_B_times_non_normal = [12.5, 12.6, 12.7, 12.8, 12.9, 13.0]
compare_extinction_points(cell_type_A_times_non_normal, cell_type_B_times_non_normal)