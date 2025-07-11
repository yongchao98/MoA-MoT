import math

def solve_methylation_and_expression():
    """
    Calculates the initial histone methylation percentage and the resulting
    impact on gene expression based on the provided parameters.
    """
    
    # --- Part 1: Calculate the initial percentage of H3K4me3 ---
    
    final_proportion_me3 = 0.1104  # 11.04%
    time_hours = 10
    turnover_rate = 0.10  # 10% per hour

    # The process follows continuous exponential decay: P(t) = P0 * e^(-k*t)
    # We need to find P0 (initial proportion).
    # P0 = P(t) / e^(-k*t)
    
    initial_proportion_me3 = final_proportion_me3 / math.exp(-turnover_rate * time_hours)

    print("--- Part 1: Initial Methylation Calculation ---")
    print(f"The decay of H3K4me3 is modeled by the equation: P(t) = P0 * e^(-k*t)")
    print(f"To find the initial proportion (P0), we use the rearranged formula: P0 = P(t) / e^(-k*t)")
    print("Plugging in the given values:")
    print(f"P0 = {final_proportion_me3 * 100:.2f}% / e^(-{turnover_rate} * {time_hours})")
    print(f"The calculated initial percentage of H3K4me3 sites is: {initial_proportion_me3 * 100:.2f}%\n")

    # --- Part 2: Determine the impact on gene expression ---
    
    initial_expression_rpkm = 200
    relative_decrease_me3 = 0.10  # 10% decrease

    # Based on the linear relationship, a relative change in methylation causes the same relative change in expression.
    # New Expression = Initial Expression * (1 - Relative Decrease)
    
    new_expression_rpkm = initial_expression_rpkm * (1 - relative_decrease_me3)
    
    print("--- Part 2: Gene Expression Impact Calculation ---")
    print("Assuming a linear relationship between H3K4me3 proportion and gene expression.")
    print(f"The initial methylation of {initial_proportion_me3 * 100:.2f}% corresponds to an expression level of {initial_expression_rpkm} RPKM.")
    print(f"When the proportion of H3K4me3 sites decreases by {relative_decrease_me3 * 100:.0f}%, the expression level also decreases by the same relative amount.")
    print("The new expression level is calculated as follows:")
    print(f"New Expression = Initial Expression * (1 - Relative Decrease)")
    print(f"New Expression = {initial_expression_rpkm} * (1 - {relative_decrease_me3})")
    print(f"The new, impacted gene expression level is: {new_expression_rpkm:.1f} RPKM")

# Execute the function to print the results
solve_methylation_and_expression()

<<<180.0>>>