def solve_linear_logic_encoding():
    """
    This script provides the formula F(r, z, nz) for encoding the JZDEC
    instruction of a Minsky machine in Linear Logic.
    
    It prints the final formula using common text representations for the
    logical operators.
    """

    # Define text representations for Linear Logic symbols
    # We use 'r', 'z', and 'nz' as symbolic placeholders for the
    # register index and state labels.
    r = 'r'
    z = 'z'
    nz = 'nz'

    # 'otimes' represents the multiplicative conjunction (tensor, ⊗)
    # '&' represents the additive conjunction (with, &)
    # 'multimap' represents linear implication (⊸)
    tensor = "otimes"
    with_op = "&"
    implication = "multimap"

    # Define literals
    state_literal = "S"
    counter_literal = "C"
    zero_test_literal = "Z"

    # Construct the formula for the zero case: Z_r ⊗ S_z
    # This branch is only viable if counter r is zero.
    zero_case = f"({zero_test_literal}_{r} {tensor} {state_literal}_{z})"

    # Construct the formula for the non-zero case: C_r ⊸ S_nz
    # This branch is only viable if counter r is greater than zero.
    nonzero_case = f"({counter_literal}_{r} {implication} {state_literal}_{nz})"

    # Combine the two cases with the '&' connective
    # F(r, z, nz) = (zero_case) & (non_zero_case)
    final_formula = f"{zero_case} {with_op} {nonzero_case}"

    print("The appropriate formula for F(r,z,nz) is:")
    print(final_formula)

solve_linear_logic_encoding()