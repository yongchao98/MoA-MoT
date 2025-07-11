import numpy as np
from scipy import stats

def compare_extinction_points(times_a, times_b, group_a_name="Cell Type A", group_b_name="Cell Type B"):
    """
    Compares the extinction times of two groups to determine if they are significantly different.
    
    This function follows the correct statistical procedure:
    1. Checks if the data for each group is normally distributed using the Shapiro-Wilk test.
    2. If both are normal, it performs an independent t-test.
    3. If at least one is not normal, it performs a Mann-Whitney U test.
    
    Args:
        times_a (list or np.array): A list of extinction times for the first group.
        times_b (list or np.array): A list of extinction times for the second group.
        group_a_name (str): The name of the first group for printing.
        group_b_name (str): The name of the second group for printing.
    """
    print(f"--- Comparing {group_a_name} vs. {group_b_name} ---")
    
    # Ensure data are numpy arrays for calculations
    times_a = np.array(times_a)
    times_b = np.array(times_b)
    
    print(f"Data for {group_a_name} (n={len(times_a)}): {np.round(times_a, 2)}")
    print(f"Data for {group_b_name} (n={len(times_b)}): {np.round(times_b, 2)}\n")

    # Step 1: Check for normality using the Shapiro-Wilk test
    print("Step 1: Checking for Normality")
    shapiro_a_stat, shapiro_a_p = stats.shapiro(times_a)
    shapiro_b_stat, shapiro_b_p = stats.shapiro(times_b)
    
    print(f"Shapiro-Wilk test for {group_a_name}: Statistic={shapiro_a_stat:.4f}, p-value={shapiro_a_p:.4f}")
    print(f"Shapiro-Wilk test for {group_b_name}: Statistic={shapiro_b_stat:.4f}, p-value={shapiro_b_p:.4f}")

    # A p-value > 0.05 suggests the data is normally distributed.
    is_a_normal = shapiro_a_p > 0.05
    is_b_normal = shapiro_b_p > 0.05

    print("\nStep 2: Selecting and Performing the Appropriate Significance Test")
    # Step 2: Choose and perform the test
    if is_a_normal and is_b_normal:
        print("Result: Both datasets appear normally distributed (p > 0.05).")
        print("Action: Performing an Independent t-test.\n")
        
        # Perform the independent t-test
        t_stat, p_value = stats.ttest_ind(times_a, times_b, equal_var=True) # Assuming equal variances for simplicity, can be changed
        
        print("--- Independent t-test Results ---")
        print(f"t-statistic: {t_stat:.4f}")
        print(f"p-value: {p_value:.4f}")

    else:
        if not is_a_normal and not is_b_normal:
            reason = "Both datasets"
        elif not is_a_normal:
            reason = f"The {group_a_name} dataset"
        else:
            reason = f"The {group_b_name} dataset"
        
        print(f"Result: {reason} does not appear normally distributed (p <= 0.05).")
        print("Action: Performing a Mann-Whitney U test (Wilcoxon rank-sum test).\n")

        # Perform the Mann-Whitney U test
        u_stat, p_value = stats.mannwhitneyu(times_a, times_b, alternative='two-sided')

        print("--- Mann-Whitney U Test Results ---")
        print(f"U-statistic: {u_stat:.4f}")
        print(f"p-value: {p_value:.4f}")

    # Step 3: Interpret the result
    print("\n--- Conclusion ---")
    if p_value < 0.05:
        print(f"The result is statistically significant (p < 0.05).")
        print(f"The extinction points of {group_a_name} and {group_b_name} are significantly different.")
    else:
        print(f"The result is not statistically significant (p >= 0.05).")
        print(f"There is no significant difference between the extinction points of {group_a_name} and {group_b_name}.")
    print("-" * 50 + "\n")


if __name__ == '__main__':
    # --- EXAMPLE 1: Data appears Normal ---
    # Replace these lists with your actual extinction time data (e.g., in hours or days)
    # In this example, the extinction times for Cell Type A are generally higher than for B.
    extinction_times_A_normal = [20.5, 22.1, 19.4, 21.8, 23.0, 20.9, 18.7]
    extinction_times_B_normal = [16.2, 17.5, 15.8, 18.1, 17.0, 16.6, 15.5]
    compare_extinction_points(extinction_times_A_normal, extinction_times_B_normal, "Cell Type A (Normal)", "Cell Type B (Normal)")

    # --- EXAMPLE 2: Data does NOT appear Normal ---
    # Here, we use data that will likely fail the normality test to show the other path.
    # Group C has a skewed distribution.
    extinction_times_C_skewed = [12, 13, 13.5, 14, 15, 25, 28] 
    extinction_times_D_normal = [18, 19, 20, 21, 22, 23, 24]
    compare_extinction_points(extinction_times_C_skewed, extinction_times_D_normal, "Cell Type C (Skewed)", "Cell Type D (Normal)")
