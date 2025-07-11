import math

def solve_methylation_problem():
    """
    Solves a two-part problem involving histone methylation kinetics and its effect on gene expression.
    """

    # --- Part 1: Calculate the initial percentage of trimethylated sites ---

    print("--- Part 1: Initial Trimethylation Percentage ---")
    
    # Given values for Part 1
    final_proportion_percent = 11.04
    rate = 0.10  # 10% per hour
    time = 10  # hours

    # Convert final proportion from percent to decimal for calculation
    final_proportion_decimal = final_proportion_percent / 100

    # The decay is modeled by: Final = Initial * e^(-rate * time)
    # We solve for the initial proportion: Initial = Final / e^(-rate * time)
    initial_proportion_decimal = final_proportion_decimal / math.exp(-rate * time)
    initial_proportion_percent = initial_proportion_decimal * 100

    print("The final proportion of H3K4me3 sites is calculated from the initial proportion using the formula:")
    print(f"Final Proportion = Initial Proportion * e^(-rate * time)")
    print("To find the initial proportion, we rearrange the formula:")
    print(f"Initial Proportion = Final Proportion / e^(-rate * time)")
    print(f"Initial Proportion (%) = {final_proportion_percent} / e^(-{rate} * {time})")
    print(f"Initial Proportion (%) = {final_proportion_percent} / {math.exp(-rate * time):.4f}")
    print(f"Initial Proportion (%) = {initial_proportion_percent:.2f}")
    print(f"\nTherefore, the percentage of sites trimethylated at the start of this period was {initial_proportion_percent:.2f}%.")
    
    print("\n" + "="*60 + "\n")

    # --- Part 2: Determine the impact on gene expression ---
    
    print("--- Part 2: Impact on Gene Expression ---")
    
    # Given values for Part 2
    initial_expression = 200  # RPKM
    proportion_decrease = 0.10  # 10% decrease

    # Due to the linear relationship, a 10% decrease in methylation leads to a 10% decrease in expression.
    final_expression = initial_expression * (1 - proportion_decrease)

    print("With a linear relationship between methylation and gene expression, a 10% decrease in H3K4me3 sites will cause a proportional decrease in expression.")
    print(f"New Expression = Initial Expression * (1 - Proportion Decrease)")
    print(f"New Expression = {initial_expression} RPKM * (1 - {proportion_decrease})")
    print(f"New Expression = {initial_expression} RPKM * {1 - proportion_decrease}")
    print(f"New Expression = {final_expression:.2f} RPKM")
    print(f"\nTherefore, the impact on target gene expression is a decrease to {final_expression:.2f} RPKM.")

solve_methylation_problem()

<<<30.01>>>