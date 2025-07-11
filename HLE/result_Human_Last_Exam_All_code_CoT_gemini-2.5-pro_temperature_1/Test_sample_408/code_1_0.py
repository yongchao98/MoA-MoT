import numpy as np
from scipy import stats

# This plan follows the logic of first checking for normality and then selecting the appropriate test
# to compare the extinction points of two independent cell types.

# Step 1: Define the data.
# These are the measured extinction times (e.g., in hours) for each replicate.
# We have two independent groups: Cell Type A and Cell Type B.
cell_type_A = [20.5, 21.3, 19.8, 20.9, 21.7, 20.2, 19.5]
# Let's make Cell Type B's data skewed to demonstrate the non-parametric path.
# The presence of outliers makes this sample not normally distributed.
cell_type_B = [22.1, 22.5, 21.9, 23.0, 22.3, 28.5, 29.1]

print("This script determines if the extinction points of two cell types are significantly different.")
print("It follows the correct statistical procedure: check normality, then perform the appropriate test.")
print("-" * 60)
print("Input Data:")
print(f"Extinction Times for Cell Type A: {cell_type_A}")
print(f"Extinction Times for Cell Type B: {cell_type_B}")
print("-" * 60)

# Step 2: Check for normality using the Shapiro-Wilk test.
# The null hypothesis of the Shapiro-Wilk test is that the data is normally distributed.
# We use a significance level (alpha) of 0.05.
alpha = 0.05
shapiro_A_stat, shapiro_A_p = stats.shapiro(cell_type_A)
shapiro_B_stat, shapiro_B_p = stats.shapiro(cell_type_B)

print("Step 2: Checking for Normality")
print(f"Shapiro-Wilk test for Cell Type A: p-value = {shapiro_A_p:.4f}")
print(f"Shapiro-Wilk test for Cell Type B: p-value = {shapiro_B_p:.4f}")

# Determine if the data is normal based on the p-values
is_A_normal = shapiro_A_p > alpha
is_B_normal = shapiro_B_p > alpha

if is_A_normal:
    print("-> Data for Cell Type A appears to be normally distributed (p > 0.05).")
else:
    print("-> Data for Cell Type A does NOT appear to be normally distributed (p <= 0.05).")

if is_B_normal:
    print("-> Data for Cell Type B appears to be normally distributed (p > 0.05).")
else:
    print("-> Data for Cell Type B does NOT appear to be normally distributed (p <= 0.05).")
print("-" * 60)

# Step 3: Choose and perform the appropriate significance test.
print("Step 3: Performing Significance Test")
if is_A_normal and is_B_normal:
    # If both are normal, use the unpaired t-test.
    print("Action: Both datasets appear normal. Performing an unpaired t-test.")
    # The t-test compares the means of the two groups.
    t_stat, p_value = stats.ttest_ind(cell_type_A, cell_type_B)
    print(f"T-test Results: T-statistic = {t_stat:.4f}, p-value = {p_value:.4f}")

else:
    # If at least one is not normal, use the Wilcoxon rank-sum test.
    print("Action: At least one dataset does not appear normal. Performing a Wilcoxon rank-sum test.")
    # The Wilcoxon test is a non-parametric test that compares the distributions of the two groups.
    u_stat, p_value = stats.mannwhitneyu(cell_type_A, cell_type_B, alternative='two-sided')
    print(f"Wilcoxon Rank-Sum Test Results: U-statistic = {u_stat:.1f}, p-value = {p_value:.4f}")

print("-" * 60)

# Step 4: Interpret the result.
print("Step 4: Conclusion")
if p_value < alpha:
    print(f"Since the p-value ({p_value:.4f}) is less than {alpha}, we reject the null hypothesis.")
    print("Conclusion: The difference in extinction times between the two cell types is statistically significant.")
else:
    print(f"Since the p-value ({p_value:.4f}) is not less than {alpha}, we fail to reject the null hypothesis.")
    print("Conclusion: The difference in extinction times between the two cell types is NOT statistically significant.")
