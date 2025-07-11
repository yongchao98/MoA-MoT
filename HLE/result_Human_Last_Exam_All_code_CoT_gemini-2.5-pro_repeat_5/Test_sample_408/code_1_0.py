import numpy as np
from scipy import stats

def compare_extinction_points(data_a, data_b, alpha=0.05):
    """
    Compares the extinction points of two cell types (A and B).

    This function first checks if the data for each cell type is normally distributed
    using the Shapiro-Wilk test. Based on the result, it performs either an
    unpaired t-test (if both are normal) or a Wilcoxon rank-sum test (if not)
    to determine if there is a significant difference between the two groups.

    Args:
        data_a (list or np.array): A list of extinction points for cell type A.
        data_b (list or np.array): A list of extinction points for cell type B.
        alpha (float): The significance level for the statistical tests.
    """
    print("--- Step 1: Simulating Extinction Point Data ---")
    print(f"Cell Type A Extinction Times (hours): {np.round(data_a, 2)}")
    print(f"Cell Type B Extinction Times (hours): {np.round(data_b, 2)}\n")

    # --- Step 2: Check for Normality ---
    print("--- Step 2: Checking for Normality (Shapiro-Wilk Test) ---")
    shapiro_a_stat, shapiro_a_p = stats.shapiro(data_a)
    shapiro_b_stat, shapiro_b_p = stats.shapiro(data_b)

    print(f"Cell Type A: p-value = {shapiro_a_p:.4f}")
    print(f"Cell Type B: p-value = {shapiro_b_p:.4f}\n")

    # Assume normality if p-value is greater than alpha
    is_a_normal = shapiro_a_p > alpha
    is_b_normal = shapiro_b_p > alpha

    if is_a_normal:
        print("Data for Cell Type A appears to be normally distributed.")
    else:
        print("Data for Cell Type A does not appear to be normally distributed.")

    if is_b_normal:
        print("Data for Cell Type B appears to be normally distributed.\n")
    else:
        print("Data for Cell Type B does not appear to be normally distributed.\n")

    # --- Step 3: Select and Perform Significance Test ---
    print("--- Step 3: Performing Significance Test ---")
    if is_a_normal and is_b_normal:
        print("Both datasets are normal. Performing Unpaired T-test...")
        # Welch's t-test is used by default, which doesn't assume equal variance
        t_stat, p_value = stats.ttest_ind(data_a, data_b, equal_var=False)
        print(f"T-statistic: {t_stat:.4f}")
        print(f"P-value: {p_value:.4f}\n")
    else:
        print("At least one dataset is not normal. Performing Wilcoxon Rank-Sum test...")
        rank_stat, p_value = stats.ranksums(data_a, data_b)
        print(f"Wilcoxon statistic: {rank_stat:.4f}")
        print(f"P-value: {p_value:.4f}\n")

    # --- Step 4: Conclusion ---
    print("--- Step 4: Conclusion ---")
    if p_value < alpha:
        print(f"Since the p-value ({p_value:.4f}) is less than {alpha}, we reject the null hypothesis.")
        print("Conclusion: There is a statistically significant difference between the extinction points of the two cell types.")
    else:
        print(f"Since the p-value ({p_value:.4f}) is greater than or equal to {alpha}, we fail to reject the null hypothesis.")
        print("Conclusion: There is no statistically significant difference between the extinction points of the two cell types.")


if __name__ == '__main__':
    # --- Example Data Generation ---
    # We will generate two sets of data representing extinction times in hours.
    # To demonstrate the logic, we'll create two normally distributed samples
    # with different means to show a case where a t-test is appropriate.
    np.random.seed(42) # for reproducible results
    
    # Cell Type A: normal distribution, mean=50 hours, std=5, 15 replicates
    extinction_times_A = np.random.normal(loc=50, scale=5, size=15)
    
    # Cell Type B: normal distribution, mean=60 hours, std=6, 15 replicates
    extinction_times_B = np.random.normal(loc=60, scale=6, size=15)

    # Run the comparison
    compare_extinction_points(extinction_times_A, extinction_times_B)
    
    # To see the other path, uncomment the lines below to make one dataset non-normal
    # print("\n" + "="*50 + "\n")
    # print("DEMONSTRATING NON-NORMAL DATA CASE\n")
    # extinction_times_C = stats.lognorm.rvs(s=0.9, loc=45, scale=10, size=15, random_state=42)
    # extinction_times_D = stats.lognorm.rvs(s=0.9, loc=60, scale=10, size=15, random_state=1)
    # compare_extinction_points(extinction_times_C, extinction_times_D)