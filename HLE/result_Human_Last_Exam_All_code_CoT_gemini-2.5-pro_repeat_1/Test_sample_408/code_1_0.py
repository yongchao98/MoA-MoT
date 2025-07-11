import numpy as np
from scipy import stats

def analyze_extinction_points():
    """
    Analyzes the difference between extinction points of two cell types.
    
    This function demonstrates the correct statistical workflow:
    1. Defines two sample datasets of extinction times.
    2. Checks both datasets for normality using the Shapiro-Wilk test.
    3. Based on the normality check, performs either:
        a. An independent t-test (if data is normal).
        b. A Wilcoxon rank-sum test (if data is not normal).
    4. Prints the results and interpretation.
    """
    # Sample data: Extinction times in hours for two cell types.
    # In a real scenario, you would replace this with your actual experimental data.
    cell_type_1_extinction_times = [25.1, 26.3, 24.8, 25.5, 27.1, 26.8, 25.9]
    cell_type_2_extinction_times = [28.2, 29.1, 27.8, 30.0, 29.5, 28.6, 28.9]

    print("--- Step 1: Data Definition ---")
    print(f"Cell Type 1 Extinction Times: {cell_type_1_extinction_times}")
    print(f"Cell Type 2 Extinction Times: {cell_type_2_extinction_times}\n")
    
    # --- Step 2: Check for Normality ---
    print("--- Step 2: Normality Check (Shapiro-Wilk Test) ---")
    alpha = 0.05
    shapiro_stat_1, shapiro_p_1 = stats.shapiro(cell_type_1_extinction_times)
    shapiro_stat_2, shapiro_p_2 = stats.shapiro(cell_type_2_extinction_times)

    print(f"Cell Type 1: P-value = {shapiro_p_1:.4f}")
    print(f"Cell Type 2: P-value = {shapiro_p_2:.4f}")

    is_normal_1 = shapiro_p_1 > alpha
    is_normal_2 = shapiro_p_2 > alpha

    if is_normal_1 and is_normal_2:
        print(f"\nBoth datasets appear to be normally distributed (p > {alpha}).")
        # --- Step 3a: Perform Unpaired T-test ---
        print("\n--- Step 3: Performing an Independent T-test ---\n")
        
        # Calculate components for the equation
        mean1 = np.mean(cell_type_1_extinction_times)
        mean2 = np.mean(cell_type_2_extinction_times)
        std1 = np.std(cell_type_1_extinction_times, ddof=1) # ddof=1 for sample std dev
        std2 = np.std(cell_type_2_extinction_times, ddof=1)
        n1 = len(cell_type_1_extinction_times)
        n2 = len(cell_type_2_extinction_times)
        
        # Perform the t-test using scipy
        t_stat, p_value = stats.ttest_ind(cell_type_1_extinction_times, cell_type_2_extinction_times, equal_var=False) # Welch's t-test

        print("T-test Equation: t = (mean1 - mean2) / sqrt(s1^2/n1 + s2^2/n2)")
        print("---------------------------------------------")
        print(f"Mean 1 (mean1): {mean1:.4f}")
        print(f"Mean 2 (mean2): {mean2:.4f}")
        print(f"Std Dev 1 (s1): {std1:.4f}")
        print(f"Std Dev 2 (s2): {std2:.4f}")
        print(f"Sample Size 1 (n1): {n1}")
        print(f"Sample Size 2 (n2): {n2}")
        print("---------------------------------------------")
        print(f"Plugging in the numbers:")
        print(f"t = ({mean1:.4f} - {mean2:.4f}) / sqrt({std1**2:.4f}/{n1} + {std2**2:.4f}/{n2})")
        
        print("\n--- Results ---")
        print(f"Calculated T-statistic: {t_stat:.4f}")
        print(f"P-value: {p_value:.4f}")

    else:
        print(f"\nAt least one dataset does not appear normally distributed (p <= {alpha}).")
        # --- Step 3b: Perform Wilcoxon Rank-Sum Test ---
        print("\n--- Step 3: Performing a Wilcoxon Rank-Sum Test ---\n")
        
        # The equation for Wilcoxon test is based on ranks and is not shown here.
        stat, p_value = stats.ranksums(cell_type_1_extinction_times, cell_type_2_extinction_times)
        print("--- Results ---")
        print(f"Wilcoxon Statistic: {stat:.4f}")
        print(f"P-value: {p_value:.4f}")

    # --- Step 4: Interpretation ---
    print("\n--- Step 4: Conclusion ---")
    if p_value < alpha:
        print(f"Since the p-value ({p_value:.4f}) is less than {alpha}, we reject the null hypothesis.")
        print("Conclusion: There is a statistically significant difference between the extinction points of the two cell types.")
    else:
        print(f"Since the p-value ({p_value:.4f}) is greater than {alpha}, we fail to reject the null hypothesis.")
        print("Conclusion: There is no statistically significant difference between the extinction points of the two cell types.")

if __name__ == '__main__':
    analyze_extinction_points()