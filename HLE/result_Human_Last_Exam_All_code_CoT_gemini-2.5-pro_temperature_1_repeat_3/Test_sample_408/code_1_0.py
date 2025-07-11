import numpy as np
from scipy import stats

def analyze_extinction_points():
    """
    This function demonstrates the correct statistical procedure to compare
    extinction points between two cell types.
    """
    # --- Step 1: Define the Data ---
    # In a real experiment, you would measure these values.
    # Here, we simulate the extinction times (e.g., in hours) for 10 replicates
    # of two different cell types, A and B.
    # We are testing if the two groups have significantly different extinction times.

    # Data for Cell Type A (e.g., collected from 10 experiments)
    extinction_times_A = np.array([10.5, 9.8, 11.2, 8.9, 10.1, 12.0, 9.5, 9.9, 10.8, 11.5])

    # Data for Cell Type B (e.g., collected from 10 experiments)
    extinction_times_B = np.array([12.5, 13.1, 11.9, 12.8, 13.5, 12.2, 13.9, 11.8, 12.7, 13.0])

    print("--- Step 1: The Data ---")
    print(f"Extinction times for Cell Type A: {extinction_times_A}")
    print(f"Extinction times for Cell Type B: {extinction_times_B}\n")

    # --- Step 2: Check for Normality ---
    # We use the Shapiro-Wilk test.
    # The null hypothesis (H0) is that the data is normally distributed.
    # We use a significance level (alpha) of 0.05.
    alpha = 0.05
    shapiro_stat_A, shapiro_p_A = stats.shapiro(extinction_times_A)
    shapiro_stat_B, shapiro_p_B = stats.shapiro(extinction_times_B)

    print("--- Step 2: Check for Normality (Shapiro-Wilk Test) ---")
    print(f"Cell Type A: p-value = {shapiro_p_A:.4f}")
    print(f"Cell Type B: p-value = {shapiro_p_B:.4f}")

    is_A_normal = shapiro_p_A > alpha
    is_B_normal = shapiro_p_B > alpha

    # --- Step 3: Perform Significance Test Based on Normality ---
    print("\n--- Step 3: Perform Significance Test ---")
    if is_A_normal and is_B_normal:
        print("Conclusion: Both datasets appear normal (p > 0.05). Performing an unpaired t-test.\n")
        # Perform unpaired t-test
        t_statistic, p_value = stats.ttest_ind(extinction_times_A, extinction_times_B)
        print(f"Unpaired t-test results:")
        print(f"T-statistic = {t_statistic:.4f}")
        print(f"P-value = {p_value:.4f}")

    else:
        print("Conclusion: At least one dataset is not normal (p <= 0.05). Performing a Wilcoxon rank-sum test.\n")
        # Perform Wilcoxon rank-sum test (Mann-Whitney U)
        u_statistic, p_value = stats.mannwhitneyu(extinction_times_A, extinction_times_B, alternative='two-sided')
        print(f"Wilcoxon rank-sum test results:")
        print(f"U-statistic = {u_statistic:.4f}")
        print(f"P-value = {p_value:.4f}")

    # --- Step 4: Final Conclusion ---
    print("\n--- Step 4: Final Conclusion ---")
    if p_value < alpha:
        print(f"The p-value ({p_value:.4f}) is less than {alpha}.")
        print("Result: The difference in extinction points between the two cell types is statistically significant.")
    else:
        print(f"The p-value ({p_value:.4f}) is not less than {alpha}.")
        print("Result: We cannot conclude there is a statistically significant difference in extinction points.")

# Run the analysis
analyze_extinction_points()