import numpy as np
from scipy import stats

# This script demonstrates the correct statistical procedure to compare the
# extinction points of two cell types.

# --- Step 1: Provide the extinction time data ---
# This data represents the time (e.g., in hours) it took for the cell
# population to reach extinction in several replicate experiments.
# We will use sample data that is not normally distributed to show the full workflow.

cell_type_A_extinction_times = np.array([20.1, 21.5, 22.0, 23.3, 30.8, 31.1])
cell_type_B_extinction_times = np.array([25.4, 26.1, 27.9, 28.3, 40.2, 42.5])

print("--- Step 1: Data ---")
print(f"Cell Type A Extinction Times (hours): {cell_type_A_extinction_times}")
print(f"Cell Type B Extinction Times (hours): {cell_type_B_extinction_times}")
print("-" * 30)


# --- Step 2: Check for normality using the Shapiro-Wilk test ---
# The null hypothesis (H0) is that the data is normally distributed.
# We will use a significance level (alpha) of 0.05.
alpha = 0.05
shapiro_A_stat, shapiro_A_p = stats.shapiro(cell_type_A_extinction_times)
shapiro_B_stat, shapiro_B_p = stats.shapiro(cell_type_B_extinction_times)

print("\n--- Step 2: Normality Check ---")
print(f"Shapiro-Wilk test for Cell Type A: p-value = {shapiro_A_p:.4f}")
print(f"Shapiro-Wilk test for Cell Type B: p-value = {shapiro_B_p:.4f}")


# --- Step 3: Choose and perform the appropriate significance test ---
print("\n--- Step 3: Significance Testing ---")
# The choice of test depends on the normality check results.
if shapiro_A_p > alpha and shapiro_B_p > alpha:
    # This block executes if both datasets appear normally distributed.
    print(f"Decision: Both p-values are > {alpha}, so data is treated as normal.")
    print("Performing an Unpaired t-test.")

    t_stat, p_value = stats.ttest_ind(cell_type_A_extinction_times, cell_type_B_extinction_times)

    print("\nUnpaired t-test Results:")
    print(f"Comparing dataset 1: {cell_type_A_extinction_times}")
    print(f"With dataset 2:    {cell_type_B_extinction_times}")
    print(f"t-statistic: {t_stat:.4f}")
    print(f"p-value: {p_value:.4f}")

else:
    # This block executes if at least one dataset does not appear normal.
    print(f"Decision: At least one p-value is <= {alpha}, so data is treated as non-normal.")
    print("Performing a Wilcoxon rank-sum test.")

    ranksum_stat, p_value = stats.ranksums(cell_type_A_extinction_times, cell_type_B_extinction_times)

    print("\nWilcoxon Rank-Sum Test Results:")
    print(f"Comparing dataset 1: {cell_type_A_extinction_times}")
    print(f"With dataset 2:    {cell_type_B_extinction_times}")
    print(f"Statistic: {ranksum_stat:.4f}")
    print(f"p-value: {p_value:.4f}")


# --- Step 4: Interpret the final result ---
print("\n--- Step 4: Final Conclusion ---")
if p_value < alpha:
    print(f"The p-value ({p_value:.4f}) is less than {alpha}.")
    print("Result: The difference in extinction points between the two cell types is statistically significant.")
else:
    print(f"The p-value ({p_value:.4f}) is not less than {alpha}.")
    print("Result: The difference in extinction points between the two cell types is not statistically significant.")
