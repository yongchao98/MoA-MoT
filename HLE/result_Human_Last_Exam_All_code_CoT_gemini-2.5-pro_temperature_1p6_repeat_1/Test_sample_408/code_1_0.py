import numpy as np
from scipy import stats

def compare_extinction_points(cell_type_1_data, cell_type_2_data, group1_name="Cell Type 1", group2_name="Cell Type 2"):
    """
    Performs a statistical analysis to compare the extinction points of two microbial cell types.
    
    This function:
    1. Checks if the data for each cell type is normally distributed using the Shapiro-Wilk test.
    2. If both are normal, it performs an independent (unpaired) t-test.
    3. If at least one is not normal, it performs a Wilcoxon rank-sum (Mann-Whitney U) test.
    """
    
    print(f"--- Statistical Analysis of Extinction Points: {group1_name} vs {group2_name} ---")

    # Step 1: Check for normality on both datasets
    print("\nStep 1: Checking for Normality (Shapiro-Wilk Test)")
    shapiro_stat_1, shapiro_p_1 = stats.shapiro(cell_type_1_data)
    shapiro_stat_2, shapiro_p_2 = stats.shapiro(cell_type_2_data)

    print(f"  - {group1_name}: Statistic={shapiro_stat_1:.4f}, p-value={shapiro_p_1:.4f}")
    print(f"  - {group2_name}: Statistic={shapiro_stat_2:.4f}, p-value={shapiro_p_2:.4f}")

    # A common significance level (alpha) is 0.05.
    # If p > alpha, we fail to reject the null hypothesis that the data is normal.
    is_normal_1 = shapiro_p_1 > 0.05
    is_normal_2 = shapiro_p_2 > 0.05

    # Step 2 & 3: Choose and perform the appropriate significance test
    print("\nStep 2: Performing Significance Test")
    if is_normal_1 and is_normal_2:
        print("Both datasets appear normally distributed. Using Independent t-test.")
        
        # Levene's test to check for equal variances (optional but good practice)
        levene_stat, levene_p = stats.levene(cell_type_1_data, cell_type_2_data)
        equal_variances = levene_p > 0.05
        
        # Perform independent t-test. equal_var is set based on Levene's test.
        # Welch's t-test (equal_var=False) is generally robust.
        ttest_stat, ttest_p = stats.ttest_ind(cell_type_1_data, cell_type_2_data, equal_var=equal_variances)
        
        print("\n--- T-test Results ---")
        print(f"T-statistic = {ttest_stat:.4f}")
        print(f"P-value = {ttest_p:.4f}")
        
    else:
        print("At least one dataset is not normally distributed. Using Wilcoxon rank-sum test.")
        
        # Perform Wilcoxon rank-sum test (Mann-Whitney U test)
        mwu_stat, mwu_p = stats.mannwhitneyu(cell_type_1_data, cell_type_2_data, alternative='two-sided')
        
        print("\n--- Wilcoxon Rank-Sum Test Results ---")
        print(f"U-statistic = {mwu_stat:.4f}")
        print(f"P-value = {mwu_p:.4f}")
    
    # Final Conclusion
    final_p_value = ttest_p if (is_normal_1 and is_normal_2) else mwu_p
    print("\n--- Conclusion ---")
    if final_p_value < 0.05:
        print(f"The difference in extinction points between {group1_name} and {group2_name} is statistically significant (p < 0.05).")
    else:
        print(f"The difference in extinction points between {group1_name} and {group2_name} is not statistically significant (p >= 0.05).")


# --- EXAMPLE USAGE ---
# Replace these lists with your actual extinction time data (e.g., in hours)

# Scenario 1: Assume your collected data looks like this (and happens to be normal)
np.random.seed(42) # for reproducible results
extinction_data_A = np.random.normal(loc=120, scale=15, size=25) # Mean 120h, SD 15h, 25 replicates
extinction_data_B = np.random.normal(loc=150, scale=20, size=28) # Mean 150h, SD 20h, 28 replicates

compare_extinction_points(extinction_data_A, extinction_data_B, group1_name="Resistant Strain", group2_name="Wild Type")

print("\n" + "="*80 + "\n")

# Scenario 2: Assume your data is skewed (not normal)
extinction_data_C = stats.lognorm.rvs(s=0.5, loc=80, scale=20, size=30, random_state=1)
extinction_data_D = stats.lognorm.rvs(s=0.5, loc=100, scale=25, size=32, random_state=2)

compare_extinction_points(extinction_data_C, extinction_data_D, group1_name="Mutant C", group2_name="Mutant D")
