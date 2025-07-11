import numpy as np
from scipy import stats

# This script demonstrates the correct statistical procedure to compare the
# extinction points of two microbial cell types.

# --- Step 1: Define the Data ---
# These lists represent the measured extinction times (e.g., in hours) for
# two different cell types from multiple replicate experiments.
# We will use sample data that is not normally distributed to show the full workflow.
cell_type_A_extinction_times = [20.5, 21.0, 22.3, 19.8, 20.9, 23.5, 18.2, 35.1]
cell_type_B_extinction_times = [25.1, 26.3, 25.8, 27.2, 26.5, 24.9, 25.5, 27.8]

print("--- Step 1: Input Data ---")
print(f"Cell Type A Extinction Times (hours): {cell_type_A_extinction_times}")
print(f"Cell Type B Extinction Times (hours): {cell_type_B_extinction_times}\n")

# --- Step 2: Check for Normality ---
# We use the Shapiro-Wilk test to determine if the data is normally distributed.
# The null hypothesis of this test is that the data is normal.
alpha = 0.05
print(f"--- Step 2: Normality Check (alpha = {alpha}) ---")

# Test Cell Type A
shapiro_stat_A, shapiro_p_A = stats.shapiro(cell_type_A_extinction_times)
print(f"Shapiro-Wilk test for Cell Type A: p-value = {shapiro_p_A:.4f}")

# Test Cell Type B
shapiro_stat_B, shapiro_p_B = stats.shapiro(cell_type_B_extinction_times)
print(f"Shapiro-Wilk test for Cell Type B: p-value = {shapiro_p_B:.4f}\n")

# --- Step 3: Choose and Perform the Correct Significance Test ---
# If either p-value is less than alpha, the data is not normal,
# and we must use a non-parametric test.

data_is_normal = shapiro_p_A > alpha and shapiro_p_B > alpha

print("--- Step 3: Performing Significance Test ---")
if data_is_normal:
    # If data is normal, use the independent (unpaired) t-test.
    # Welch's t-test (equal_var=False) is used as it's more robust.
    print("Decision: Data appears normal. Using Independent t-test.")
    t_stat, p_value = stats.ttest_ind(cell_type_A_extinction_times, cell_type_B_extinction_times, equal_var=False)
    print("Test Results:")
    print(f"  T-statistic = {t_stat:.4f}")
    print(f"  P-value = {p_value:.4f}")
else:
    # If data is not normal, use the Wilcoxon rank-sum test (Mann-Whitney U).
    print("Decision: Data is not normal. Using Wilcoxon rank-sum test.")
    # The 'alternative' is 'two-sided' because we are checking for any difference.
    u_stat, p_value = stats.mannwhitneyu(cell_type_A_extinction_times, cell_type_B_extinction_times, alternative='two-sided')
    print("Test Results:")
    print(f"  U-statistic = {u_stat:.4f}")
    print(f"  P-value = {p_value:.4f}")

# --- Step 4: Conclusion ---
print("\n--- Step 4: Final Conclusion ---")
if p_value < alpha:
    print(f"The p-value ({p_value:.4f}) is less than our significance level ({alpha}).")
    print("Result: We conclude that there IS a statistically significant difference between the extinction points.")
else:
    print(f"The p-value ({p_value:.4f}) is greater than our significance level ({alpha}).")
    print("Result: We CANNOT conclude that there is a statistically significant difference between the extinction points.")
