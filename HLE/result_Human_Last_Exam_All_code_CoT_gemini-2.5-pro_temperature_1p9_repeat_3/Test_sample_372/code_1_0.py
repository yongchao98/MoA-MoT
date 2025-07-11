import math

def solve_epigenetics_problem():
    """
    Solves a two-part problem related to histone methylation and gene expression.
    """
    
    # --- Part 1: Calculate the initial percentage of H3K4me3 sites ---
    
    print("--- Part 1: Initial H3K4me3 Percentage Calculation ---")
    
    # Given parameters
    rate_per_hour = 0.10  # 10%
    time_hours = 10
    final_proportion_percent = 11.04
    final_proportion = final_proportion_percent / 100
    
    # The turnover of H3K4me3 follows an exponential decay model.
    # The formula is: Final Proportion = Initial Proportion * (1 - Rate)^Time
    # We can rearrange this to solve for the Initial Proportion:
    # Initial Proportion = Final Proportion / (1 - Rate)^Time
    
    initial_proportion = final_proportion / ((1 - rate_per_hour)**time_hours)
    initial_percentage = initial_proportion * 100
    
    print("The model for this process is exponential decay.")
    print("The equation to find the initial percentage (P_initial) is: P_initial = P_final / (1 - rate)^time")
    print("\nCalculating the initial percentage of H3K4me3 sites:")
    print(f"P_initial = {final_proportion_percent}% / (1 - {rate_per_hour})^{time_hours}")
    print(f"The initial percentage of trimethylated sites was: {initial_percentage:.2f}%\n\n")

    # --- Part 2: Determine the impact on gene expression ---

    print("--- Part 2: Gene Expression Impact Calculation ---")
    
    # Given parameters
    initial_expression_rpkm = 200
    proportion_decrease_percentage = 10
    proportion_decrease = proportion_decrease_percentage / 100
    
    # A linear relationship means a % decrease in proportion leads to the same % decrease in expression.
    final_expression_factor = 1 - proportion_decrease
    final_expression_rpkm = initial_expression_rpkm * final_expression_factor
    
    print("The relationship between methylation and gene expression is linear.")
    print(f"A decrease of {proportion_decrease_percentage}% in H3K4me3 proportion will cause a {proportion_decrease_percentage}% decrease in gene expression.")
    print("The equation to find the new expression level (E_final) is: E_final = E_initial * (1 - decrease_rate)")
    print("\nCalculating the new average gene expression:")
    print(f"E_final = {initial_expression_rpkm} RPKM * (1 - {proportion_decrease})")
    print(f"The new average gene expression level would be: {final_expression_rpkm:.2f} RPKM")

solve_epigenetics_problem()