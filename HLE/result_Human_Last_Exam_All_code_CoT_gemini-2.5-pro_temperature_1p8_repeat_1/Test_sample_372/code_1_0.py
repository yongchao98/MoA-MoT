import math

def solve_methylation_problem():
    """
    Solves a two-part problem regarding histone methylation and gene expression.

    Part 1: Calculates the initial percentage of H3K4me3 sites based on a decay model.
    Part 2: Calculates the new gene expression level after a decrease in H3K4me3 sites.
    """

    # --- Part 1: Calculate the initial percentage of trimethylated sites ---

    print("--- Part 1: Initial Percentage of Trimethylated Sites ---")
    
    # Given values
    final_proportion_val = 11.04  # in percent
    decay_rate = 0.10          # 10% per hour
    time = 10                  # in hours

    # The decay model is P_t = P_0 * (1 - rate)^time
    # We solve for the initial proportion, P_0 = P_t / (1 - rate)^time
    
    final_proportion_dec = final_proportion_val / 100.0
    
    # Calculate the denominator (1 - rate)^time
    decay_factor = (1 - decay_rate) ** time
    
    # Calculate the initial proportion as a decimal
    initial_proportion_dec = final_proportion_dec / decay_factor
    
    # Convert to percentage
    initial_proportion_pct = initial_proportion_dec * 100

    print("The initial proportion is calculated using the decay formula: Final = Initial * (1 - Rate)^Time")
    print("Rearranged for the initial value: Initial = Final / (1 - Rate)^Time")
    print("\nSubstituting the given values into the equation:")
    print(f"Initial Proportion = {final_proportion_dec} / (1 - {decay_rate})^{time}")
    print(f"Initial Proportion = {final_proportion_dec} / {decay_factor:.4f}")
    print(f"The initial proportion was {initial_proportion_dec:.4f}, which is {initial_proportion_pct:.2f}%.")

    # --- Part 2: Impact on Gene Expression ---

    print("\n--- Part 2: Impact on Gene Expression ---")

    # Given values
    initial_expression = 200     # in RPKM
    proportion_decrease = 0.10   # 10% decrease

    # With a linear relationship, a 10% decrease in methylation leads to a 10% decrease in expression.
    final_expression = initial_expression * (1 - proportion_decrease)

    print("Assuming a linear relationship, a 10% decrease in methylation proportion leads to a proportional decrease in gene expression.")
    print("The new expression is calculated as: New Expression = Initial Expression * (1 - Proportional Decrease)")
    print("\nSubstituting the given values into the equation:")
    print(f"New Expression = {initial_expression} * (1 - {proportion_decrease})")
    print(f"New Expression = {initial_expression} * {1 - proportion_decrease}")
    print(f"The new gene expression level is {final_expression:.2f} RPKM.")

# Execute the function to print the results
solve_methylation_problem()