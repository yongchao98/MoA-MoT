import math

def solve_methylation_expression():
    """
    Calculates the initial methylation percentage and the impact on gene expression.
    """
    # Part 1: Calculate the initial percentage of H3K4me3
    
    # Given values
    final_proportion_percent = 11.04
    rate_per_hour = 0.10  # 10%
    time_hours = 10

    # Convert percentage to decimal for calculation
    final_proportion = final_proportion_percent / 100

    # The formula for exponential decay is P(t) = P(0) * (1 - r)^t
    # We solve for the initial proportion, P(0) = P(t) / (1 - r)^t
    initial_proportion = final_proportion / ((1 - rate_per_hour)**time_hours)
    initial_proportion_percent = initial_proportion * 100

    print("--- Part 1: Initial Percentage of Trimethylated Sites ---")
    print("The final proportion is calculated from the initial proportion using the exponential decay formula:")
    print(f"Final Proportion = Initial Proportion * (1 - Decay Rate)^Time")
    print("To find the initial proportion, we rearrange the formula:")
    print(f"Initial Proportion = Final Proportion / (1 - Decay Rate)^Time")
    print("\nPlugging in the numbers:")
    print(f"Initial Proportion (%) = {final_proportion_percent} / (1 - {rate_per_hour})^{time_hours}")
    print(f"Initial Proportion (%) = {initial_proportion_percent:.2f}%")
    
    print("\n" + "="*60 + "\n")

    # Part 2: Determine the impact on gene expression
    
    # Given values
    initial_expression = 200  # RPKM
    proportion_decrease_percent = 10 # 10%
    proportion_decrease = proportion_decrease_percent / 100

    # Due to the linear relationship, a 10% decrease in methylation proportion
    # leads to a 10% decrease in gene expression.
    new_expression = initial_expression * (1 - proportion_decrease)
    expression_change = new_expression - initial_expression

    print("--- Part 2: Impact on Gene Expression ---")
    print("Given a linear relationship between methylation and expression, a change in methylation proportion directly impacts expression by the same percentage.")
    print("\nCalculating the new expression level after a 10% decrease:")
    print(f"New Expression = Initial Expression * (1 - {proportion_decrease})")
    print(f"New Expression = {initial_expression} RPKM * (1 - {proportion_decrease}) = {new_expression:.2f} RPKM")

    print("\nThe impact is the change in expression:")
    print(f"Change in Expression = New Expression - Initial Expression")
    print(f"Change in Expression = {new_expression:.2f} RPKM - {initial_expression} RPKM = {expression_change:.2f} RPKM")

# Execute the function to print the results
solve_methylation_expression()

<<< -20.00 >>>