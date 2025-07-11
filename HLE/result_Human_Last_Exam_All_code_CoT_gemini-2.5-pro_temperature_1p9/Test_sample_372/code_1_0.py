import math

def solve_epigenetics_problem():
    """
    Solves a two-part problem involving histone methylation decay and its impact on gene expression.
    """

    # --- Part 1: Calculate the initial percentage of H3K4me3 sites ---

    print("--- Part 1: Finding the Initial Trimethylation Percentage ---")

    # Given values for Part 1
    p_final_percent = 11.04  # Percentage of H3K4me3 after 10 hours
    time_hours = 10.0          # Time elapsed in hours
    rate_per_hour = 0.10     # 10% turnover rate per hour

    print("The turnover of H3K4me3 to H3K4me2 is a continuous process modeled by the exponential decay formula:")
    print("P(t) = P(0) * e^(-r*t)\n")
    print("Where:")
    print(f"  P(t) = Percentage after time t = {p_final_percent}%")
    print("  P(0) = Initial percentage (the value we need to find)")
    print("  e   = Euler's number (base of natural logarithm)")
    print(f"  r   = Rate of turnover = {rate_per_hour} per hour")
    print(f"  t   = Time = {time_hours} hours\n")

    # To find the initial percentage P(0), we rearrange the formula:
    # P(0) = P(t) / e^(-r*t)
    p_initial_percent = p_final_percent / math.exp(-rate_per_hour * time_hours)

    print("Calculation:")
    print(f"Step 1: Substitute the known values into the equation.")
    print(f"  {p_final_percent} = P(0) * e^(-{rate_per_hour} * {time_hours})")
    print(f"  {p_final_percent} = P(0) * e^({-rate_per_hour * time_hours})")
    print("\nStep 2: Solve for P(0).")
    print(f"  P(0) = {p_final_percent} / e^({-rate_per_hour * time_hours})")
    print(f"  P(0) = {p_initial_percent:.2f}\n")
    print(f"Therefore, the percentage of sites trimethylated at the start was {p_initial_percent:.2f}%.\n")


    # --- Part 2: Determine the impact on gene expression ---

    print("--- Part 2: Determining the Impact on Gene Expression ---")

    # Given values for Part 2
    initial_expression_rpkm = 200.0  # Average expression level in RPKM
    proportion_decrease = 0.10     # 10% decrease in H3K4me3 sites

    print("A linear relationship is assumed between the proportion of H3K4me3 sites and gene expression.")
    print("This means gene expression (E) is directly proportional to the methylation proportion (P): E = k * P.")
    print("Therefore, a 10% decrease in the methylation proportion will result in a 10% decrease in gene expression.\n")

    # Calculate the new gene expression level
    final_expression_rpkm = initial_expression_rpkm * (1 - proportion_decrease)

    print("Calculation:")
    print(f"Step 1: Define the initial state and the change.")
    print(f"  Initial Expression = {initial_expression_rpkm} RPKM")
    print(f"  Decrease in Methylation = {proportion_decrease * 100}%")
    print("\nStep 2: Apply the proportional decrease to the expression level.")
    print(f"  Final Expression = Initial Expression * (1 - Decrease in Methylation)")
    print(f"  Final Expression = {initial_expression_rpkm} * (1 - {proportion_decrease})")
    print(f"  Final Expression = {final_expression_rpkm:.2f}\n")
    print(f"The impact on the target gene expression is a reduction to {final_expression_rpkm:.2f} RPKM.")

solve_epigenetics_problem()
<<<180.00>>>