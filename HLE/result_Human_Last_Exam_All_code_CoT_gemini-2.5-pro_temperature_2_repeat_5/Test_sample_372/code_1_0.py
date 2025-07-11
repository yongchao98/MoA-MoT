import math

def solve_methylation_problem():
    """
    Solves a two-part problem related to histone methylation and gene expression.
    """
    
    # Part 1: Calculate the initial percentage of H3K4me3 sites.
    print("--- Part 1: Calculating the initial percentage of H3K4me3 ---")
    
    p_after_10h = 11.04  # Percentage after 10 hours
    rate_per_hour = 0.10 # 10% turnover rate per hour
    time_hours = 10      # 10 hours duration
    
    # The formula for exponential decay is P(t) = P0 * (1 - r)^t
    # We need to solve for P0: P0 = P(t) / (1 - r)^t
    
    decay_factor = (1 - rate_per_hour)**time_hours
    initial_proportion = p_after_10h / decay_factor
    
    print("The formula for the final proportion P(t) is: P(t) = P0 * (1 - r)^t")
    print(f"Given: P({time_hours}) = {p_after_10h}%, r = {rate_per_hour}, t = {time_hours} hours")
    print("To find the initial proportion (P0), we rearrange the formula: P0 = P(t) / (1 - r)^t")
    print(f"P0 = {p_after_10h} / (1 - {rate_per_hour})^{time_hours}")
    print(f"P0 = {p_after_10h} / {decay_factor:.4f}")
    print(f"The initial percentage of trimethylated sites (P0) was: {initial_proportion:.2f}%\n")

    # Part 2: Determine the impact on gene expression.
    print("--- Part 2: Calculating the impact on gene expression ---")

    initial_expression = 200  # RPKM
    decrease_fraction = 0.10  # 10% decrease in the proportion of H3K4me3 sites

    # Assuming a linear (proportional) relationship: Expression = k * Proportion
    # A 10% decrease in proportion leads to a 10% decrease in expression.
    # New Expression = Initial Expression * (1 - Decrease Fraction)
    
    final_expression = initial_expression * (1 - decrease_fraction)
    
    print("Assuming a linear relationship between methylation and gene expression.")
    print(f"The initial expression level is {initial_expression} RPKM.")
    print("If the H3K4me3 proportion decreases by 10%, the expression level will also decrease proportionally.")
    print("Final Expression = Initial Expression * (1 - Decrease Fraction)")
    print(f"Final Expression = {initial_expression} * (1 - {decrease_fraction})")
    print(f"The new gene expression level would be: {final_expression:.2f} RPKM")


solve_methylation_problem()
<<<Initial methylation: 31.66%, New expression: 180.00 RPKM>>>