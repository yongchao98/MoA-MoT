import math

def solve_epigenetics_problem():
    """
    Calculates the initial methylation percentage and the resulting
    impact on gene expression based on the given parameters.
    """
    print("### Part 1: Calculating the initial percentage of trimethylated sites ###\n")

    # Given values for Part 1
    final_prop_perc = 11.04
    rate_per_hour = 0.10  # 10%
    time_hours = 10

    # Convert percentage to a proportion
    final_prop = final_prop_perc / 100.0

    # Using the exponential decay formula: P(t) = P(0) * e^(-k*t)
    # We solve for P(0): P(0) = P(t) / e^(-k*t)
    initial_prop = final_prop / math.exp(-rate_per_hour * time_hours)
    initial_prop_perc = initial_prop * 100

    print("The formula to find the initial proportion P(0) is: P(t) / e^(-k*t)")
    print(f"Using the given values:")
    print(f"P(0) = {final_prop} / e^(-{rate_per_hour} * {time_hours})")
    print(f"P(0) = {initial_prop:.4f}\n")
    print(f"Therefore, the initial percentage of trimethylated sites was {initial_prop_perc:.2f}%.\n")


    print("### Part 2: Calculating the impact on gene expression ###\n")

    # Given values for Part 2
    initial_expression = 200  # RPKM
    methylation_decrease = 0.10 # 10%

    # The new proportion is a 10% decrease from the initial proportion found in Part 1
    new_prop = initial_prop * (1 - methylation_decrease)

    # Assuming a linear relationship: New_Expression / Initial_Expression = New_Proportion / Initial_Proportion
    # New_Expression = Initial_Expression * (New_Proportion / Initial_Proportion)
    new_expression = initial_expression * (new_prop / initial_prop)

    print("A 10% decrease in the initial methylation of "
          f"{initial_prop_perc:.2f}% results in a new proportion of "
          f"{initial_prop_perc * (1 - methylation_decrease):.2f}%.")
    print("\nBased on a linear relationship between methylation and gene expression, the new expression level is calculated as:")
    print(f"New Expression = Initial Expression * (New Proportion / Initial Proportion)")
    print(f"New Expression = {initial_expression} RPKM * ({new_prop:.4f} / {initial_prop:.4f})")
    print(f"New Expression = {initial_expression} RPKM * {1 - methylation_decrease}")
    print(f"\nThe impact on target gene expression results in a new level of {new_expression:.0f} RPKM.")

solve_epigenetics_problem()