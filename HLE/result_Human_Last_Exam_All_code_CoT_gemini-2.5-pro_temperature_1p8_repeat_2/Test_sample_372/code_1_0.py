import math

def solve_methylation_and_expression():
    """
    Calculates the initial histone methylation percentage and the resulting
    change in gene expression based on the provided parameters.
    """

    # --- Part 1: Calculate the initial percentage of H3K4me3 ---
    print("--- Part 1: Initial H3K4me3 Percentage Calculation ---")
    print("The turnover of H3K4me3 is modeled as an exponential decay process.")
    print("The formula is: Final_Prop = Initial_Prop * (1 - Decay_Rate)^Time\n")

    # Given parameters for Part 1
    final_prop = 0.1104  # 11.04%
    decay_rate = 0.10    # 10% per hour
    time = 10.0          # 10 hours

    # We need to solve for Initial_Prop
    # Rearranged formula: Initial_Prop = Final_Prop / (1 - Decay_Rate)^Time
    
    # Perform the calculation
    initial_prop = final_prop / ((1 - decay_rate) ** time)
    initial_prop_percent = initial_prop * 100

    print("To find the initial percentage, we use the following values:")
    print(f"  Final Proportion (after 10 hours): {final_prop * 100:.2f}%")
    print(f"  Decay Rate per hour: {decay_rate * 100:.0f}%")
    print(f"  Time: {time:.0f} hours\n")

    print("Calculation:")
    print(f"Initial Proportion = {final_prop} / (1 - {decay_rate})^{int(time)}")
    print(f"Initial Proportion = {final_prop} / {math.pow(1 - decay_rate, time):.6f}")
    print(f"Initial Proportion = {initial_prop:.6f}\n")
    print(f"Result: The percentage of sites trimethylated at the start of this period was {initial_prop_percent:.2f}%.\n")


    # --- Part 2: Determine the impact on gene expression ---
    print("--- Part 2: Gene Expression Impact Calculation ---")
    print("The problem assumes a linear relationship between H3K4me3 proportion and gene expression.")
    print("Therefore, a percentage decrease in methylation will result in the same percentage decrease in expression.\n")

    # Given parameters for Part 2
    initial_expression = 200.0  # RPKM
    proportion_decrease = 0.10  # 10%

    # Calculate the new expression level
    new_expression = initial_expression * (1 - proportion_decrease)
    change_in_expression = initial_expression - new_expression

    print("To find the impact on gene expression, we use the following values:")
    print(f"  Initial Gene Expression: {initial_expression:.0f} RPKM")
    print(f"  Decrease in H3K4me3 Proportion: {proportion_decrease * 100:.0f}%\n")

    print("Calculation:")
    print(f"New Expression = Initial Expression * (1 - Proportion Decrease)")
    print(f"New Expression = {initial_expression:.0f} * (1 - {proportion_decrease})")
    print(f"New Expression = {initial_expression:.0f} * {1 - proportion_decrease}\n")
    print(f"Result: The new target gene expression level is {new_expression:.2f} RPKM.")
    print(f"This is a decrease of {change_in_expression:.2f} RPKM from the initial level.")


# Execute the function to print the solutions
solve_methylation_and_expression()
<<<Initial H3K4me3 was 31.66%, and the new expression level is 180.00 RPKM.>>>