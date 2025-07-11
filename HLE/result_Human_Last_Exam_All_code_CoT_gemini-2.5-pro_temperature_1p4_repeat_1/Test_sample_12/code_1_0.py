import sympy

def solve_transformer_cost_ratio():
    """
    Calculates the ratio of additional computational cost for increasing
    a Transformer's context length from L to 4L.
    """
    # Define the symbolic variables based on the problem description
    L, d_model, d_attn, d_ff = sympy.symbols('L d_model d_attn d_ff', positive=True)

    # 1. Define the computational FLOPs for a single forward pass through one layer.
    # The cost consists of a quadratic part (from self-attention scores) and a linear part
    # (from projections and feed-forward networks).

    # Quadratic part (from QK^T and Attn@V): Proportional to L^2.
    # We use the approximation 2 * L^2 * d_attn.
    flops_quadratic = 2 * L**2 * d_attn

    # Linear part (from FFN and all attention linear projections): Proportional to L.
    # - Feed-Forward Network (FFN): ~ 4 * L * d_model * d_ff
    # - Attention Projections (Q, K, V, and Output): ~ 8 * L * d_model * d_attn
    flops_linear = 8 * L * d_model * d_attn + 4 * L * d_model * d_ff

    # Total FLOPs for a single layer pass
    total_flops = flops_quadratic + flops_linear

    # 2. Calculate the cost per token.
    # To train on the same total number of tokens, the cost is proportional to the FLOPs per token.
    # Cost per token is Total FLOPs / L.
    cost_per_token_L = total_flops / L

    # 3. Calculate the new cost per token for a context length of 4L.
    cost_per_token_4L = cost_per_token_L.subs(L, 4 * L)

    # 4. Calculate the additional cost.
    additional_cost = cost_per_token_4L - cost_per_token_L

    # 5. Calculate the ratio of the additional cost to the original cost.
    ratio = additional_cost / cost_per_token_L

    # 6. Simplify the final expression to match the options.
    simplified_ratio = sympy.simplify(ratio)

    print("The ratio of the additional computational cost to the original cost is:")
    # The pretty print function displays the formula clearly.
    sympy.pprint(simplified_ratio, use_unicode=False)

if __name__ == '__main__':
    solve_transformer_cost_ratio()