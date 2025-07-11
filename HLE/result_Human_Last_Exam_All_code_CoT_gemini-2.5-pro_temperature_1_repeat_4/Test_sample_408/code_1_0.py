import numpy as np
from scipy import stats

# Step 1: Define sample data.
# These lists represent the collected extinction times (e.g., in hours) from multiple
# experiments for two different microbial cell types under a stressor.
cell_type_A_extinction_times = [48.5, 51.2, 49.8, 52.1, 50.5, 47.9, 53.0, 51.5, 50.1, 49.2]
cell_type_B_extinction_times = [42.3, 45.1, 41.8, 44.5, 43.2, 46.0, 40.9, 43.8, 44.1, 42.9]
# We'll use a standard significance level (alpha).
alpha = 0.05

print("This script determines if the extinction points of two cell types are significantly different.")
print("-" * 70)
print("Step 1: The data represents the extinction time (in hours) for each replicate.")
print(f"Extinction times for Cell Type A: {cell_type_A_extinction_times}")
print(f"Extinction times for Cell Type B: {cell_type_B_extinction_times}\n")


# Step 2: Check if the data for each group is normally distributed.
# We use the Shapiro-Wilk test. The null hypothesis is that the data is normal.
print("Step 2: Checking for normality using the Shapiro-Wilk test.")
shapiro_A_stat, shapiro_A_p = stats.shapiro(cell_type_A_extinction_times)
shapiro_B_stat, shapiro_B_p = stats.shapiro(cell_type_B_extinction_times)

print(f"Cell Type A: p-value = {shapiro_A_p:.4f}")
print(f"Cell Type B: p-value = {shapiro_B_p:.4f}\n")


# Step 3: Choose and perform the appropriate significance test.
is_A_normal = shapiro_A_p > alpha
is_B_normal = shapiro_B_p > alpha

# If both p-values are > alpha, we assume normality and use a t-test.
if is_A_normal and is_B_normal:
    print(f"Analysis: Both distributions appear normal (p > {alpha}).")
    print("Performing an unpaired t-test to compare the means.\n")

    # Perform the unpaired t-test
    ttest_stat, ttest_p = stats.ttest_ind(cell_type_A_extinction_times, cell_type_B_extinction_times)

    print("--- T-test Results ---")
    print(f"The two data sets being compared are:")
    print(f"Set 1: {cell_type_A_extinction_times}")
    print(f"Set 2: {cell_type_B_extinction_times}")
    print(f"T-statistic: {ttest_stat:.4f}")
    print(f"P-value: {ttest_p:.4f}")

    if ttest_p < alpha:
        print("\nConclusion: The result is statistically significant.")
    else:
        print("\nConclusion: The result is not statistically significant.")

# Otherwise, we use the non-parametric Wilcoxon rank-sum test.
else:
    print(f"Analysis: At least one distribution appears non-normal (p <= {alpha}).")
    print("Performing a Wilcoxon rank-sum test to compare distributions.\n")

    # Perform the Wilcoxon rank-sum test
    mwu_stat, mwu_p = stats.mannwhitneyu(cell_type_A_extinction_times, cell_type_B_extinction_times)

    print("--- Wilcoxon Rank-Sum Test Results ---")
    print(f"The two data sets being compared are:")
    print(f"Set 1: {cell_type_A_extinction_times}")
    print(f"Set 2: {cell_type_B_extinction_times}")
    print(f"U-statistic: {mwu_stat:.4f}")
    print(f"P-value: {mwu_p:.4f}")

    if mwu_p < alpha:
        print("\nConclusion: The result is statistically significant.")
    else:
        print("\nConclusion: The result is not statistically significant.")
print("-" * 70)
