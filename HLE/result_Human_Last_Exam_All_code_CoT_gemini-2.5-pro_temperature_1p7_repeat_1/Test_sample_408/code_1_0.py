import numpy as np
from scipy import stats

def analyze_extinction_times():
    """
    This script determines if the extinction point of two different microbial 
    cell types are significantly different.

    It implements the following logic:
    1. It uses two sample arrays of data representing the time-to-extinction
       for each cell type.
    2. It checks if each data sample is normally distributed using the Shapiro-Wilk test.
    3. If both are normal, it performs an unpaired t-test.
    4. If at least one is not normal, it performs a Wilcoxon rank-sum test.
    5. It prints the results and a final conclusion.
    """
    
    # --- Sample Data ---
    # Imagine you ran 12 replicates for Cell Type A and 12 for Cell Type B.
    # The numbers below are the hours it took for each replicate to go extinct.
    # We will use slightly skewed data to demonstrate the non-parametric path.
    np.random.seed(42) # for reproducibility
    extinction_times_A = np.random.lognormal(mean=4.0, sigma=0.25, size=12)
    extinction_times_B = np.random.lognormal(mean=4.2, sigma=0.25, size=12)

    # These are the data points that will be used in the test
    print("--- Data for Statistical Test ---")
    # Using np.round to make the output cleaner
    print(f"Extinction Times for Cell Type A (hours): {np.round(extinction_times_A, 2)}")
    print(f"Extinction Times for Cell Type B (hours): {np.round(extinction_times_B, 2)}\n")
    
    # --- Step 1: Check for Normality ---
    print("--- Step 1: Checking for Normality (Shapiro-Wilk Test) ---")
    alpha = 0.05
    
    stat_A, p_val_A = stats.shapiro(extinction_times_A)
    stat_B, p_val_B = stats.shapiro(extinction_times_B)
    
    print(f"Cell Type A: p-value = {p_val_A:.4f}")
    print(f"Cell Type B: p-value = {p_val_B:.4f}")

    # --- Step 2: Choose and Perform the Correct Significance Test ---
    is_A_normal = p_val_A > alpha
    is_B_normal = p_val_B > alpha

    if is_A_normal and is_B_normal:
        print(f"\nConclusion: Both samples appear normally distributed (p > {alpha}).")
        print("--- Step 2: Performing Unpaired t-test ---")
        
        # Levene's test for equality of variances is good practice before a t-test
        levene_stat, levene_p = stats.levene(extinction_times_A, extinction_times_B)
        # Use Welch's t-test if variances are not equal
        equal_variances = levene_p > alpha
        if not equal_variances:
            print("(Note: Variances are unequal, using Welch's t-test)")

        t_statistic, p_value_final = stats.ttest_ind(extinction_times_A, extinction_times_B, equal_var=equal_variances)
        
        print(f"Test Statistic (t-value): {t_statistic:.4f}")
        print(f"Final P-value: {p_value_final:.4f}")
    
    else:
        print(f"\nConclusion: At least one sample is not normally distributed (p <= {alpha}).")
        print("--- Step 2: Performing Non-Parametric Wilcoxon Rank-Sum Test ---")
        
        statistic, p_value_final = stats.mannwhitneyu(extinction_times_A, extinction_times_B, alternative='two-sided')

        print(f"Test Statistic (U-value): {statistic:.4f}")
        print(f"Final P-value: {p_value_final:.4f}")

    # --- Step 3: Final Conclusion ---
    print("\n--- Step 3: Interpreting the Result ---")
    if p_value_final < alpha:
        print(f"Since the p-value ({p_value_final:.4f}) is less than our significance level ({alpha}), we reject the null hypothesis.")
        print("Final Answer: The extinction times between Cell Type A and Cell Type B are significantly different.")
    else:
        print(f"Since the p-value ({p_value_final:.4f}) is greater than our significance level ({alpha}), we fail to reject the null hypothesis.")
        print("Final Answer: There is not enough evidence to say the extinction times are significantly different.")

# Run the analysis function
analyze_extinction_times()