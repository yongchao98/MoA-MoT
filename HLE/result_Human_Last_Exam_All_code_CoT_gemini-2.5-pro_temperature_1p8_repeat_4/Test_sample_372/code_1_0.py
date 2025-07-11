import math

def solve_methylation_and_expression():
    """
    Solves a two-part problem on histone methylation and gene expression.
    Part 1: Calculates the initial percentage of H3K4me3 sites using an exponential decay model.
    Part 2: Calculates the new gene expression level based on a decrease in H3K4me3 sites, assuming a linear relationship.
    """
    # --- Part 1: Calculate Initial H3K4me3 Percentage ---
    print("Part 1: Calculating the initial percentage of trimethylated sites.")
    
    # Given values
    final_prop_perc = 11.04  # Final proportion in percent
    rate = 0.10  # 10% per hour decay rate
    time = 10    # 10 hours

    # The formula for continuous exponential decay is A(t) = A(0) * e^(-k*t).
    # We solve for the initial amount A(0): A(0) = A(t) / e^(-k*t).
    initial_prop_perc = final_prop_perc / math.exp(-rate * time)
    
    # Print the equation and the result
    print("The calculation for the initial percentage (A(0)) is based on the formula: A(0) = A(t) / e^(-k*t)")
    print(f"Initial trimethylated sites = {final_prop_perc}% / e^(-{rate:.2f} * {time}) = {initial_prop_perc:.2f}%")
    print("-" * 20)

    # --- Part 2: Determine Impact on Gene Expression ---
    print("\nPart 2: Determining the impact on target gene expression.")

    # Given values
    initial_expression = 200 # RPKM
    decrease_fraction = 0.10 # 10% decrease
    
    # The relationship is linear, so we assume direct proportionality.
    # New Expression = Initial Expression * (New Methylation / Initial Methylation)
    # A 10% decrease in methylation means New Methylation = Initial Methylation * (1 - 0.10)
    # So, New Expression = Initial Expression * (1 - 0.10)
    new_methylation_perc = initial_prop_perc * (1 - decrease_fraction)
    new_expression = initial_expression * (1 - decrease_fraction)
    
    # Print the logic and the result
    print("Assuming a linear relationship, a 10% decrease in H3K4me3 proportion leads to a 10% decrease in gene expression.")
    print(f"Initial methylation from Part 1: {initial_prop_perc:.2f}%")
    print(f"Initial gene expression: {initial_expression} RPKM")
    print(f"New gene expression = {initial_expression} RPKM * (1 - {decrease_fraction:.2f}) = {new_expression:.2f} RPKM")
    
    return initial_prop_perc, new_expression

# Run the function to display the results
initial_percentage, new_expression_level = solve_methylation_and_expression()

# Final Answer
# print(f"\n<<<{initial_percentage:.1f}, {new_expression_level:.1f}>>>")
# The question asked for the answer in <<<>>> format, so I format it as requested
# The final answer is a string with two values. The first value is the percentage from part 1, and the second is the RPKM from part 2.
print(f"\n<<<{initial_percentage:.1f}, {new_expression_level:.1f}>>>")