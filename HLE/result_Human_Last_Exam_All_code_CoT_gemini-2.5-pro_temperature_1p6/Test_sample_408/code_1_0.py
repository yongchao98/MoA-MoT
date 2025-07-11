import numpy as np
from scipy import stats

def compare_extinction_points(data_a, data_b, alpha=0.05):
    """
    Checks for normality and then performs the appropriate statistical test
    to compare two independent groups of extinction time data.

    Args:
        data_a (array-like): Extinction times for cell type A.
        data_b (array-like): Extinction times for cell type B.
        alpha (float): The significance level for the tests.
    """
    print("--- Step 1: Data Summary ---")
    mean_a, std_a, n_a = np.mean(data_a), np.std(data_a, ddof=1), len(data_a)
    mean_b, std_b, n_b = np.mean(data_b), np.std(data_b, ddof=1), len(data_b)
    print(f"Cell Type A: n={n_a}, mean={mean_a:.2f}, std_dev={std_a:.2f}")
    print(f"Cell Type B: n={n_b}, mean={mean_b:.2f}, std_dev={std_b:.2f}\n")

    print("--- Step 2: Normality Check (Shapiro-Wilk Test) ---")
    shapiro_a_stat, shapiro_a_p = stats.shapiro(data_a)
    shapiro_b_stat, shapiro_b_p = stats.shapiro(data_b)
    print(f"Cell Type A Normality p-value: {shapiro_a_p:.4f}")
    print(f"Cell Type B Normality p-value: {shapiro_b_p:.4f}")

    is_normal_a = shapiro_a_p > alpha
    is_normal_b = shapiro_b_p > alpha

    if is_normal_a and is_normal_b:
        print(f"\nConclusion: Both datasets appear normally distributed (p > {alpha}).")
        print("--- Step 3: Performing Unpaired T-test ---\n")
        
        # Perform unpaired t-test
        t_stat, p_value = stats.ttest_ind(data_a, data_b, equal_var=False) # Welch's t-test
        
        print("Unpaired T-test Results:")
        print(f"T-statistic = {t_stat:.4f}")
        print(f"P-value = {p_value:.4f}")
        
    else:
        print(f"\nConclusion: At least one dataset is not normally distributed (p <= {alpha}).")
        print("--- Step 3: Performing Wilcoxon Rank-Sum Test ---\n")

        # Perform Wilcoxon rank-sum test
        statistic, p_value = stats.ranksums(data_a, data_b)
        
        print("Wilcoxon Rank-Sum Test Results:")
        print(f"Test Statistic = {statistic:.4f}")
        print(f"P-value = {p_value:.4f}")

    print("\n--- Step 4: Final Conclusion ---")
    if p_value < alpha:
        print(f"Since the p-value ({p_value:.4f}) is less than {alpha}, we reject the null hypothesis.")
        print("The difference in extinction points between the two cell types is statistically significant.")
    else:
        print(f"Since the p-value ({p_value:.4f}) is greater than {alpha}, we fail to reject the null hypothesis.")
        print("There is not enough evidence to claim a significant difference in extinction points.")

# --- Example Usage ---
# Let's create some sample data representing the time (in hours) to extinction.
# We will use a log-normal distribution to simulate skewed, non-normal data.
np.random.seed(42) # for reproducible results
extinction_times_A = np.random.lognormal(mean=3.5, sigma=0.4, size=30)
extinction_times_B = np.random.lognormal(mean=4.0, sigma=0.5, size=35)

# Run the analysis
compare_extinction_points(extinction_times_A, extinction_times_B)