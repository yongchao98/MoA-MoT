import numpy as np
from scipy import stats

def analyze_extinction_point_difference():
    """
    Demonstrates the statistical procedure for comparing the extinction points 
    of two cell types.

    This code follows these steps:
    1. Defines the lists of extinction points for two cell types (A and B).
    2. Checks if each list of data is normally distributed using the Shapiro-Wilk test.
    3. Based on the normality check, chooses the correct statistical test:
       - Unpaired t-test if data is normal.
       - Wilcoxon rank-sum (Mann-Whitney U) test if not normal.
    4. Performs the test and prints a conclusion about the significance.
    """
    print("--- Step 1: Define Data ---")
    # In a real experiment, these values would be the measured extinction times
    # from multiple replicates for each cell type.
    extinction_points_A = np.array([10.2, 11.1, 9.8, 10.5, 11.5, 9.5, 10.8, 11.2])
    extinction_points_B = np.array([12.5, 13.1, 12.0, 13.5, 11.9, 12.8, 13.3, 12.9])

    print("Data consists of the distribution of extinction points from replicates.")
    print(f"Extinction Points (hours) for Cell Type A: {extinction_points_A}")
    print(f"Extinction Points (hours) for Cell Type B: {extinction_points_B}\n")
    
    print("--- Summary Statistics ---")
    mean_A = np.mean(extinction_points_A)
    std_A = np.std(extinction_points_A, ddof=1)
    mean_B = np.mean(extinction_points_B)
    std_B = np.std(extinction_points_B, ddof=1)
    
    print(f"Cell Type A: Mean = {mean_A:.2f}, Std Dev = {std_A:.2f}")
    print(f"Cell Type B: Mean = {mean_B:.2f}, Std Dev = {std_B:.2f}\n")

    print("--- Step 2: Check for Normality ---")
    # We use the Shapiro-Wilk test. H0: The data is normally distributed.
    # We assume normality if p > 0.05.
    alpha = 0.05
    _, p_value_A = stats.shapiro(extinction_points_A)
    _, p_value_B = stats.shapiro(extinction_points_B)

    is_A_normal = p_value_A > alpha
    is_B_normal = p_value_B > alpha

    print(f"Shapiro-Wilk test for Cell Type A p-value: {p_value_A:.4f}")
    print(f"Shapiro-Wilk test for Cell Type B p-value: {p_value_B:.4f}")

    if is_A_normal and is_B_normal:
        print("Both distributions can be assumed to be normal.\n")
    else:
        print("At least one distribution does not appear to be normal.\n")

    print("--- Step 3: Perform Significance Test ---")
    if is_A_normal and is_B_normal:
        print("Since data is normal, performing an Unpaired t-test.")
        # The numbers in the final calculation are the means and standard deviations shown above.
        test_statistic, p_value_final = stats.ttest_ind(extinction_points_A, extinction_points_B)
        print(f"Test Result: T-statistic = {test_statistic:.4f}, p-value = {p_value_final:.4f}")
    else:
        print("Since data is not normal, performing a Wilcoxon rank-sum test.")
        test_statistic, p_value_final = stats.mannwhitneyu(extinction_points_A, extinction_points_B)
        print(f"Test Result: U-statistic = {test_statistic}, p-value = {p_value_final:.4f}")

    print("\n--- Step 4: Conclusion ---")
    if p_value_final < alpha:
        print(f"The p-value ({p_value_final:.4f}) is less than {alpha}.")
        print("Conclusion: The extinction points of the two cell types are significantly different.")
    else:
        print(f"The p-value ({p_value_final:.4f}) is not less than {alpha}.")
        print("Conclusion: There is no significant difference between the extinction points.")

# Run the full analysis
analyze_extinction_point_difference()