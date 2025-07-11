import numpy as np
from scipy import stats

def compare_extinction_points():
    """
    Demonstrates the statistical procedure to compare extinction points 
    of two microbial cell types.
    """
    # --- Step 1: Gather Data ---
    # These are example datasets representing the time (e.g., in hours) to extinction
    # for multiple replicates of two different cell types.
    # In a real scenario, you would replace this with your actual experimental data.
    
    # Let's simulate data where Cell Type A tends to go extinct earlier than B.
    # We will use a log-normal distribution to simulate data that is not perfectly normal,
    # which is common for biological time-to-event data.
    np.random.seed(42) # for reproducible results
    extinction_times_A = np.random.lognormal(mean=2.3, sigma=0.2, size=15)
    extinction_times_B = np.random.lognormal(mean=2.5, sigma=0.25, size=18)
    
    print("--- Data ---")
    print(f"Cell Type A Extinction Times (hours): {[round(t, 1) for t in extinction_times_A]}")
    print(f"Cell Type B Extinction Times (hours): {[round(t, 1) for t in extinction_times_B]}\n")
    
    # Significance level
    alpha = 0.05
    
    # --- Step 2: Check for Normality using Shapiro-Wilk test ---
    print("--- Normality Check ---")
    shapiro_p_A = stats.shapiro(extinction_times_A).pvalue
    shapiro_p_B = stats.shapiro(extinction_times_B).pvalue
    
    print(f"Shapiro-Wilk test p-value for Cell Type A: {shapiro_p_A:.4f}")
    print(f"Shapiro-Wilk test p-value for Cell Type B: {shapiro_p_B:.4f}\n")
    
    is_normal_A = shapiro_p_A > alpha
    is_normal_B = shapiro_p_B > alpha

    # --- Step 3: Select and Perform the Appropriate Test ---
    print("--- Significance Test ---")
    if is_normal_A and is_normal_B:
        print("Both datasets appear normally distributed. Performing Unpaired T-test.")
        # Perform independent t-test
        t_stat, p_value = stats.ttest_ind(extinction_times_A, extinction_times_B, equal_var=False) # Welch's t-test
        print(f"T-statistic: {t_stat:.4f}")
        print(f"P-value: {p_value:.4f}")
    else:
        print("At least one dataset does not appear normally distributed. Performing Wilcoxon rank-sum test.")
        # Perform Wilcoxon rank-sum test (Mann-Whitney U)
        u_stat, p_value = stats.mannwhitneyu(extinction_times_A, extinction_times_B, alternative='two-sided')
        print(f"U-statistic: {u_stat:.4f}")
        print(f"P-value: {p_value:.4f}")

    # --- Step 4: Interpret the Results ---
    print("\n--- Conclusion ---")
    if p_value < alpha:
        print(f"The p-value ({p_value:.4f}) is less than {alpha}, so we reject the null hypothesis.")
        print("The extinction points of the two cell types are significantly different.")
    else:
        print(f"The p-value ({p_value:.4f}) is not less than {alpha}, so we fail to reject the null hypothesis.")
        print("There is not enough evidence to claim a significant difference in extinction points.")

# Run the analysis
compare_extinction_points()