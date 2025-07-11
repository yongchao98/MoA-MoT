import math

def solve_epigenetics_problem():
    """
    Solves a two-part problem related to histone methylation and gene expression.
    """
    
    # --- Part 1: Calculate the initial percentage of H3K4me3 sites ---
    
    # Given values for Part 1
    p_t_percent = 11.04  # Percentage of H3K4me3 sites after 10 hours
    p_t = p_t_percent / 100.0  # Convert percentage to a proportion
    decay_rate = 0.10  # 10% turnover rate per hour
    time = 10  # Time in hours

    # The decay is modeled by the equation: P(t) = P0 * (1 - r)^t
    # We solve for the initial proportion, P0: P0 = P(t) / (1 - r)^t
    initial_proportion = p_t / math.pow((1 - decay_rate), time)
    initial_percentage = initial_proportion * 100

    print("--- Part 1: Initial Percentage of Trimethylated Sites ---")
    print("The formula to find the initial proportion (P0) from the final proportion (P(t)) is:")
    print("P0 = P(t) / (1 - r)^t\n")
    print("Substituting the given values into the equation:")
    print(f"P0 = {p_t} / (1 - {decay_rate})^{time}")
    print(f"The calculated initial percentage of H3K4me3 sites is: {initial_percentage:.2f}%\n")

    # --- Part 2: Determine the impact on gene expression ---
    
    # Given values for Part 2
    initial_expression = 200  # RPKM
    proportion_decrease = 0.10  # 10% decrease

    # Due to the linear relationship, a 10% decrease in the methylation proportion
    # results in a 10% decrease in gene expression.
    # New Expression = Initial Expression * (1 - Decrease Percentage)
    new_expression = initial_expression * (1 - proportion_decrease)

    print("--- Part 2: Impact on Gene Expression ---")
    print("The formula to find the new gene expression (E_new) is:")
    print("E_new = E_initial * (1 - proportion_decrease)\n")
    print("Substituting the given values into the equation:")
    print(f"E_new = {initial_expression} * (1 - {proportion_decrease})")
    print(f"The new gene expression level is: {new_expression:.2f} RPKM")

solve_epigenetics_problem()
<<<31.66>>>