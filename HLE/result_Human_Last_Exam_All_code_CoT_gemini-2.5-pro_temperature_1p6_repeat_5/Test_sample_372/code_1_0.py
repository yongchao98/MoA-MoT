import math

def solve_epigenetics_problem():
    """
    Solves a two-part problem related to histone methylation and gene expression.
    """
    
    # --- Part 1: Calculate the initial proportion of H3K4me3 sites ---

    # Given values for Part 1
    final_proportion_p1 = 11.04  # in percent
    turnover_rate = 0.10         # 10% per hour
    time_hours = 10              # 10 hours

    # The turnover is an exponential decay process. The formula is:
    # Final_Proportion = Initial_Proportion * (1 - turnover_rate)^time
    # We rearrange it to solve for the Initial_Proportion:
    # Initial_Proportion = Final_Proportion / (1 - turnover_rate)^time
    
    decay_factor_base = 1 - turnover_rate
    decay_factor_total = math.pow(decay_factor_base, time_hours)
    initial_proportion_p1 = final_proportion_p1 / decay_factor_total

    print("--- Part 1: Initial H3K4me3 Proportion ---")
    print(f"The proportion of H3K4me3 sites is {final_proportion_p1:.2f}% after {time_hours} hours, with a turnover rate of {turnover_rate*100:.0f}% per hour.")
    print("To find the initial proportion (P_initial), we use the formula: P_initial = P_final / (1 - rate)^time")
    print(f"P_initial = {final_proportion_p1:.2f} / (1 - {turnover_rate:.2f}) ** {time_hours}")
    print(f"P_initial = {final_proportion_p1:.2f} / {decay_factor_base:.1f} ** {time_hours}")
    print(f"P_initial = {final_proportion_p1:.2f} / {decay_factor_total:.4f}")
    print(f"The initial percentage of trimethylated sites was: {initial_proportion_p1:.2f}%\n")


    # --- Part 2: Calculate the impact on gene expression ---

    # Given values for Part 2
    initial_expression = 200    # in RPKM
    proportion_decrease = 0.10  # 10% or 0.10

    # The relationship between methylation and expression is linear.
    # A 10% decrease in methylation leads to a 10% decrease in expression.
    # New_Expression = Initial_Expression * (1 - proportion_decrease)
    
    new_expression = initial_expression * (1 - proportion_decrease)
    
    print("--- Part 2: Impact on Gene Expression ---")
    print(f"The initial average expression level is {initial_expression} RPKM.")
    print(f"The proportion of H3K4me3 sites decreases by {proportion_decrease*100:.0f}% over 10 hours.")
    print("Assuming a linear relationship, the new expression level is calculated as:")
    print("New Expression = Initial Expression * (1 - Proportion Decrease)")
    print(f"New Expression = {initial_expression} * (1 - {proportion_decrease:.2f})")
    print(f"New Expression = {initial_expression} * {1 - proportion_decrease:.2f}")
    print(f"The new target gene expression level is: {new_expression:.2f} RPKM.")

# Execute the function to print the solution
solve_epigenetics_problem()