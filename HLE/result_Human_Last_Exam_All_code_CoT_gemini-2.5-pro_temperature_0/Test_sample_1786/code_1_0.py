def solve():
    """
    This function constructs and prints the Linear Logic formula for F(r,z,nz).
    """
    # Define the components of the formula as strings
    # S_z: Literal for the 'zero' target state
    # S_nz: Literal for the 'non-zero' target state
    # C_r: Literal for a unit of counter 'r'
    # Z_r: Literal to initiate the zero-test on counter 'r'
    # ⊗: Tensor (multiplicative conjunction)
    # &: With (additive conjunction)
    # ⊸: Linear Implication
    
    # The case where counter r is non-zero (x_r > 0).
    # It consumes one C_r and produces S_nz.
    f_nonzero = "(C_r ⊸ S_{nz})"
    
    # The case where counter r is zero (x_r = 0).
    # It produces the new state S_z and the zero-test literal Z_r.
    f_zero = "(S_z ⊗ Z_r)"
    
    # The full formula F(r,z,nz) is an additive conjunction (a choice)
    # between the zero case and the non-zero case.
    final_formula = f"{f_zero} & {f_nonzero}"
    
    print("The appropriate formula F(r,z,nz) is:")
    print(final_formula)

solve()