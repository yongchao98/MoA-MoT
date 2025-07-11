import numpy as np
from scipy import stats

def analyze_extinction_points():
    """
    Analyzes the statistical significance of the difference between extinction points
    of two microbial cell types.
    """
    # --- Step 1: Define the Data ---
    # These are example datasets representing the measured extinction points (e.g., time in hours)
    # for two different cell types under a stressor.
    # We use a fixed seed for reproducibility of the random data.
    np.random.seed(0)
    
    # Let's assume Cell Type B is more resistant.
    # Data is generated from a normal distribution for demonstration purposes.
    cell_type_A_extinction = np.random.normal(loc=20, scale=2.5, size=15)
    cell_type_B_extinction = np.random.normal(loc=24, scale=2.8, size=15)

    print("--- Step 1: Extinction Point Data ---")
    print(f"Cell Type A Extinction Points (hours): {[round(x, 2) for x in cell_type_A_extinction]}")
    print(f"Cell Type B Extinction Points (hours): {[round(x, 2) for x in cell_type_B_extinction]}")
    print("-" * 35)

    # --- Step 2: Check for Normality ---
    print("\n--- Step 2: Checking for Normality (Shapiro-Wilk Test) ---")
    alpha = 0.05
    
    shapiro_A = stats.shapiro(cell_type_A_extinction)
    shapiro_B = stats.shapiro(cell_type_B_extinction)

    print(f"Cell Type A: p-value = {shapiro_A.pvalue:.4f}")
    print(f"Cell Type B: p-value = {shapiro_B.pvalue:.4f}")

    # Check if we should assume normality
    is_normal_A = shapiro_A.pvalue > alpha
    is_normal_B = shapiro_B.pvalue > alpha

    if is_normal_A and is_normal_B:
        print(f"\nBoth p-values are > {alpha}, so we assume the data is normally distributed.")
        is_normal = True
    else:
        print(f"\nAt least one p-value is <= {alpha}, so we cannot assume normality.")
        is_normal = False
    print("-" * 35)

    # --- Step 3: Perform Significance Test ---
    print("\n--- Step 3: Performing Significance Test ---")

    if is_normal:
        print("Data is normal. Performing an Independent (Unpaired) T-test.")

        # Calculate components for the t-test equation
        mean_A = np.mean(cell_type_A_extinction)
        mean_B = np.mean(cell_type_B_extinction)
        std_A = np.std(cell_type_A_extinction, ddof=1)
        std_B = np.std(cell_type_B_extinction, ddof=1)
        n_A = len(cell_type_A_extinction)
        n_B = len(cell_type_B_extinction)

        print("\nT-test Equation Components:")
        print(f"Mean (Group A): {mean_A:.4f}")
        print(f"Mean (Group B): {mean_B:.4f}")
        print(f"Standard Deviation (Group A): {std_A:.4f}")
        print(f"Standard Deviation (Group B): {std_B:.4f}")
        print(f"Sample Size (n_A): {n_A}")
        print(f"Sample Size (n_B): {n_B}")
        
        # Perform the t-test
        t_statistic, p_value = stats.ttest_ind(cell_type_A_extinction, cell_type_B_extinction)

        print("\nTest Result:")
        print(f"T-statistic = {t_statistic:.4f}")
        print(f"P-value = {p_value:.4f}")

    else:
        print("Data is not normal. Performing the Wilcoxon Rank-Sum Test.")
        
        # Perform the Wilcoxon rank-sum test
        statistic, p_value = stats.ranksums(cell_type_A_extinction, cell_type_B_extinction)

        print("\nTest Result:")
        print(f"Wilcoxon statistic = {statistic:.4f}")
        print(f"P-value = {p_value:.4f}")
    
    print("-" * 35)

    # --- Step 4: Conclusion ---
    print("\n--- Step 4: Conclusion ---")
    if p_value < alpha:
        print(f"The p-value ({p_value:.4f}) is less than {alpha}.")
        print("We reject the null hypothesis.")
        print("Conclusion: The extinction points of the two cell types are significantly different.")
    else:
        print(f"The p-value ({p_value:.4f}) is greater than or equal to {alpha}.")
        print("We fail to reject the null hypothesis.")
        print("Conclusion: There is not enough evidence to say the extinction points are significantly different.")

if __name__ == '__main__':
    analyze_extinction_points()
