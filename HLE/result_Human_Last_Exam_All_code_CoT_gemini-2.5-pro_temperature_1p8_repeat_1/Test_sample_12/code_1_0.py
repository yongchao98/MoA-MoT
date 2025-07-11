def solve():
    """
    This function prints the derived formula for the ratio of additional computational cost.
    """
    # Define symbolic variable names as strings
    L = "L"
    d_attn = "d_attn"
    d_model = "d_model"
    d_ff = "d_ff"

    # Construct the numerator and denominator strings based on the derivation
    numerator_str = f"3 * {L} * {d_attn}"
    denominator_str = f"2 * {d_model} * (2 * {d_attn} + {d_ff}) + {L} * {d_attn}"

    # Print the full expression for the ratio
    print("The best estimate for the ratio of the additional computational cost is:")
    print(f"Ratio = ({numerator_str}) / ({denominator_str})")

solve()